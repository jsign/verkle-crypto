const std = @import("std");
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Banderwagon = banderwagon.Banderwagon;
const Fr = banderwagon.Fr;
const crs = @import("../crs/crs.zig");
const Transcript = @import("transcript.zig");
const common = @import("common.zig");
const assert = std.debug.assert;

pub fn IPA(comptime VECTOR_LENGTH: comptime_int) type {
    assert(std.math.isPowerOfTwo(VECTOR_LENGTH));

    return struct {
        const Self = @This();

        const NUM_STEPS = std.math.log2(VECTOR_LENGTH);

        pub const ProverQuery = struct {
            // TODO: rename to A
            polynomial: [VECTOR_LENGTH]Fr,
            commitment: Banderwagon,
            // Input point
            point: Fr,
            // TODO: Rename to B
            point_evaluations: [VECTOR_LENGTH]Fr,
        };

        pub const VerifierQuery = struct {
            commitment: Banderwagon,
            point: Fr,
            // TODO: rename to B
            point_evaluations: [VECTOR_LENGTH]Fr,
            // TODO: rename
            output_point: Fr,
            proof: IPAProof,
        };

        pub const IPAProof = struct {
            L: [NUM_STEPS]Banderwagon,
            R: [NUM_STEPS]Banderwagon,
            a: Fr,
        };

        const ProofResult = struct {
            result: Fr,
            proof: IPAProof,
        };

        pub fn createProof(xcrs: crs.CRS, transcript: *Transcript, query: ProverQuery) ProofResult {
            transcript.domainSep("ipa");
            transcript.appendPoint(query.commitment, "C");
            transcript.appendScalar(query.point, "input point");

            var _A = query.polynomial;
            var A: []Fr = _A[0..];
            var _B = query.point_evaluations;
            var B: []Fr = _B[0..];
            const y = common.innerProduct(A, B);
            transcript.appendScalar(y, "output point");

            // Rescale Q.
            const w = transcript.challengeScalar("w");
            const q = xcrs.Q.scalarMul(w);

            var L: [NUM_STEPS]Banderwagon = undefined;
            var R: [NUM_STEPS]Banderwagon = undefined;
            var _basis = xcrs.Gs;
            var basis: []Banderwagon = _basis[0..];

            var step: usize = 0;
            var n: usize = VECTOR_LENGTH;
            while (n > 1) {
                const m = n / 2;

                const a_L = A[0..m];
                const a_R = A[m..];
                const b_L = B[0..m];
                const b_R = B[m..];
                const z_L = common.innerProduct(a_R, b_L);
                const z_R = common.innerProduct(a_L, b_R);

                var C_L: Banderwagon = undefined;
                C_L.add(Banderwagon.msm(basis[0..m], a_R), q.scalarMul(z_L));
                var C_R: Banderwagon = undefined;
                C_R.add(Banderwagon.msm(basis[m..], a_L), q.scalarMul(z_R));

                L[step] = C_L;
                R[step] = C_R;

                transcript.appendPoint(C_L, "L");
                transcript.appendPoint(C_R, "R");
                const x = transcript.challengeScalar("x");

                var xinv: Fr = x.inv().?;

                // Compute updates for next round
                for (0..m) |i| {
                    A[i] = a_L[i].add(x.mul(a_R[i]));
                }
                A = A[0..m];

                for (b_L, b_R, 0..) |v, z, i| {
                    B[i] = v.add(xinv.mul(z));
                }
                B = B[0..m];

                for (basis[0..m], basis[m..], 0..) |bl, br, i| {
                    basis[i].add(bl, Banderwagon.scalarMul(br, xinv));
                }
                basis = basis[0..m];

                n = m;
                step += 1;
            }

            return .{ .result = y, .proof = IPAProof{ .a = A[0], .L = L, .R = R } };
        }

        pub fn verifyProof(xcrs: crs.CRS, transcript: *Transcript, query: VerifierQuery) bool {
            transcript.domainSep("ipa");

            const C = query.commitment;
            transcript.appendPoint(C, "C");

            const z = query.point;
            transcript.appendScalar(z, "input point");

            const proof = query.proof;
            const y = query.output_point;
            transcript.appendScalar(y, "output point");

            const w = transcript.challengeScalar("w");
            const q = xcrs.Q.scalarMul(w);

            var current_commitment: Banderwagon = undefined;
            current_commitment.add(C, q.scalarMul(y));

            var xs: [NUM_STEPS]Fr = undefined;
            var xinvs: [NUM_STEPS]Fr = undefined;

            var step: usize = 0;
            var n: usize = VECTOR_LENGTH;
            while (n > 1) {
                const m = n / 2;

                const C_L = proof.L[step];
                const C_R = proof.R[step];
                transcript.appendPoint(C_L, "L");
                transcript.appendPoint(C_R, "R");

                const x = transcript.challengeScalar("x");
                const x_inv = x.inv().?;

                xs[step] = x;
                xinvs[step] = x_inv;

                var tmp: Banderwagon = undefined;
                tmp.add(Banderwagon.scalarMul(C_L, x), Banderwagon.scalarMul(C_R, x_inv));
                current_commitment.add(current_commitment, tmp);

                n = m;
                step = step + 1;
            }

            // Do it the inefficient way (TODO: optimize)
            var _Gs = xcrs.Gs;
            var current_basis: []Banderwagon = _Gs[0..];

            var _b = query.point_evaluations;
            var b: []Fr = _b[0..];
            for (0..xs.len) |j| {
                const current_basis_split = splitInHalf(Banderwagon, current_basis);
                const b_split = splitInHalf(Fr, b);
                const x_inv = xinvs[j];
                foldScalars(b_split.L, b_split.R, x_inv);
                b = b_split.L;
                foldPoints(current_basis_split.L, current_basis_split.R, x_inv);
                current_basis = current_basis_split.L;
            }

            assert(b.len == current_basis.len);
            assert(b.len == 1);

            const b_0 = b[0];
            const G_0 = current_basis[0];

            // G[0] * a + (a * b) * Q;
            var got_commitment: Banderwagon = undefined;
            got_commitment.add(Banderwagon.scalarMul(G_0, proof.a), Banderwagon.scalarMul(q, Fr.mul(proof.a, b_0)));

            return current_commitment.equal(got_commitment);
        }

        // foldScalars computes a[i] = a[i] + b[i] * challenge
        fn foldScalars(a: []Fr, b: []const Fr, folding_challenge: Fr) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a[i] = Fr.add(a[i], Fr.mul(b[i], folding_challenge));
            }
        }

        // foldPoints computes a[i] = a[i] + b[i] * challenge
        fn foldPoints(a: []Banderwagon, b: []const Banderwagon, folding_challenge: Fr) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a[i].add(a[i], Banderwagon.scalarMul(b[i], folding_challenge));
            }
        }

        // Splits a list into two lists of equal length
        // Eg[S1, S2, S3, S4] becomes[S1, S2], [S3, S4]
        fn splitInHalf(comptime T: type, x: []T) struct { L: []T, R: []T } {
            assert(x.len % 2 == 0);
            const mid = x.len / 2;
            return .{ .L = x[0..mid], .R = x[mid..] };
        }
    };
}

test "basic proof" {
    const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights(crs.DomainSize, crs.Domain);
    const ipa = IPA(crs.DomainSize);

    // Test a simple IPA proof
    var weights = PrecomputedWeights.init();

    // Polynomial in lagrange basis
    var lagrange_poly: [crs.DomainSize]Fr = undefined;
    for (0..lagrange_poly.len) |i| {
        lagrange_poly[i] = Fr.fromInteger(i % 32 + 1);
    }

    // Commit to the polynomial in lagrange basis
    const xcrs = crs.CRS.init();
    const commitment = xcrs.commit(lagrange_poly);

    const expected_comm = std.fmt.bytesToHex(commitment.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("1b9dff8f5ebbac250d291dfe90e36283a227c64b113c37f1bfb9e7a743cdb128", &expected_comm);

    var prover_transcript = Transcript.init("test");

    // create a opening proof for a point outside of the domain
    const input_point = Fr.fromInteger(2101);
    const b = weights.barycentricFormulaConstants(input_point);
    const output_point_check = common.innerProduct(&lagrange_poly, &b);
    const output_point_check_hex = std.fmt.bytesToHex(output_point_check.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("4a353e70b03c89f161de002e8713beec0d740a5e20722fd5bd68b30540a33208", &output_point_check_hex);

    var query = ipa.ProverQuery{
        .polynomial = lagrange_poly,
        .commitment = commitment,
        .point = input_point,
        .point_evaluations = b,
    };

    var ipa_proof = ipa.createProof(xcrs, &prover_transcript, query);

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");
    const p_challenge_hex = std.fmt.bytesToHex(p_challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("0a81881cbfd7d7197a54ebd67ed6a68b5867f3c783706675b34ece43e85e7306", &p_challenge_hex);

    // Verify the proof.
    var verifier_transcript = Transcript.init("test");
    const b_verifier = weights.barycentricFormulaConstants(input_point);
    const verifier_query = ipa.VerifierQuery{
        .commitment = commitment,
        .point = input_point,
        .point_evaluations = b_verifier,
        .output_point = ipa_proof.result,
        .proof = ipa_proof.proof,
    };
    const ok = ipa.verifyProof(xcrs, &verifier_transcript, verifier_query);
    try std.testing.expect(ok);
}
