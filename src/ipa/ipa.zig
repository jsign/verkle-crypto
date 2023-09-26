const std = @import("std");
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const Fr = banderwagon.Fr;
const crs = @import("../crs/crs.zig");
const Transcript = @import("transcript.zig");
const assert = std.debug.assert;

pub fn IPA(comptime VectorLength: comptime_int) type {
    assert(std.math.isPowerOfTwo(VectorLength));

    return struct {
        const Self = @This();

        const NUM_STEPS = std.math.log2(VectorLength);

        pub const ProverQuery = struct {
            commitment: Element,
            A: [VectorLength]Fr,
            B: [VectorLength]Fr,
            eval_point: Fr,
        };

        pub const VerifierQuery = struct {
            commitment: Element,
            B: [VectorLength]Fr,
            eval_point: Fr,
            result: Fr,
            proof: IPAProof,
        };

        pub const IPAProof = struct {
            L: [NUM_STEPS]Element,
            R: [NUM_STEPS]Element,
            a: Fr,
        };

        const ProofResult = struct {
            result: Fr,
            proof: IPAProof,
        };

        pub fn createProof(xcrs: crs.CRS, transcript: *Transcript, query: ProverQuery) ProofResult {
            transcript.domainSep("ipa");
            transcript.appendPoint(query.commitment, "C");
            transcript.appendScalar(query.eval_point, "input point");

            var _A = query.A;
            var A: []Fr = _A[0..];
            var _B = query.B;
            var B: []Fr = _B[0..];
            const y = innerProduct(A, B);
            transcript.appendScalar(y, "output point");

            // Rescale Q.
            const w = transcript.challengeScalar("w");
            const q = xcrs.Q.scalarMul(w);

            var L: [NUM_STEPS]Element = undefined;
            var R: [NUM_STEPS]Element = undefined;
            var _basis = xcrs.Gs;
            var basis: []Element = _basis[0..];

            var step: usize = 0;
            var n: usize = VectorLength;
            while (n > 1) {
                const m = n / 2;

                const a_L = A[0..m];
                const a_R = A[m..];
                const b_L = B[0..m];
                const b_R = B[m..];
                const z_L = innerProduct(a_R, b_L);
                const z_R = innerProduct(a_L, b_R);

                var C_L: Element = undefined;
                C_L.add(banderwagon.msm(basis[0..m], a_R), q.scalarMul(z_L));
                var C_R: Element = undefined;
                C_R.add(banderwagon.msm(basis[m..], a_L), q.scalarMul(z_R));

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
                    basis[i].add(bl, Element.scalarMul(br, xinv));
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

            const z = query.eval_point;
            transcript.appendScalar(z, "input point");

            const proof = query.proof;
            const y = query.result;
            transcript.appendScalar(y, "output point");

            const w = transcript.challengeScalar("w");
            const q = xcrs.Q.scalarMul(w);

            var current_commitment: Element = undefined;
            current_commitment.add(C, q.scalarMul(y));

            var xs: [NUM_STEPS]Fr = undefined;
            var xinvs: [NUM_STEPS]Fr = undefined;

            var step: usize = 0;
            var n: usize = VectorLength;
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

                var tmp: Element = undefined;
                tmp.add(Element.scalarMul(C_L, x), Element.scalarMul(C_R, x_inv));
                current_commitment.add(current_commitment, tmp);

                n = m;
                step = step + 1;
            }

            // Do it the inefficient way (TODO: optimize)
            var _Gs = xcrs.Gs;
            var current_basis: []Element = _Gs[0..];

            var _b = query.B;
            var b: []Fr = _b[0..];
            for (0..xs.len) |j| {
                const current_basis_split = splitInHalf(Element, current_basis);
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
            var got_commitment: Element = undefined;
            got_commitment.add(Element.scalarMul(G_0, proof.a), Element.scalarMul(q, Fr.mul(proof.a, b_0)));

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
        fn foldPoints(a: []Element, b: []const Element, folding_challenge: Fr) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a[i].add(a[i], Element.scalarMul(b[i], folding_challenge));
            }
        }

        // Splits a list into two lists of equal length
        // Eg[S1, S2, S3, S4] becomes[S1, S2], [S3, S4]
        fn splitInHalf(comptime T: type, x: []T) struct { L: []T, R: []T } {
            assert(x.len % 2 == 0);
            const mid = x.len / 2;
            return .{ .L = x[0..mid], .R = x[mid..] };
        }

        fn innerProduct(a: []const Fr, b: []const Fr) Fr {
            var result = Fr.zero();
            for (a, b) |ai, bi| {
                const term = ai.mul(bi);
                result = Fr.add(result, term);
            }
            return result;
        }
    };
}

test "inner product smoke" {
    var a = [_]Fr{ Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    var b = [_]Fr{ Fr.fromInteger(10), Fr.fromInteger(12), Fr.fromInteger(13), Fr.fromInteger(14), Fr.fromInteger(15) };

    // Expected result should be 1*10 + 2*12 + 3*13 + 4*14 + 5*15
    const expected_result = Fr.fromInteger(204);
    const got_result = IPA(crs.DomainSize).innerProduct(&a, &b);
    try std.testing.expect(got_result.equal(expected_result));
}

test "basic proof" {
    const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights(crs.DomainSize, crs.Domain);
    const ipa = IPA(crs.DomainSize);

    // Test a simple IPA proof
    var weights = try PrecomputedWeights.init();

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
    const eval_point = Fr.fromInteger(2101);
    const b = try weights.barycentricFormulaConstants(eval_point);
    const output_point_check = ipa.innerProduct(&lagrange_poly, &b);
    const output_point_check_hex = std.fmt.bytesToHex(output_point_check.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("4a353e70b03c89f161de002e8713beec0d740a5e20722fd5bd68b30540a33208", &output_point_check_hex);

    var query = ipa.ProverQuery{
        .commitment = commitment,
        .A = lagrange_poly,
        .B = b,
        .eval_point = eval_point,
    };

    var ipa_proof = ipa.createProof(xcrs, &prover_transcript, query);

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");
    const p_challenge_hex = std.fmt.bytesToHex(p_challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("0a81881cbfd7d7197a54ebd67ed6a68b5867f3c783706675b34ece43e85e7306", &p_challenge_hex);

    // Verify the proof.
    var verifier_transcript = Transcript.init("test");
    const b_verifier = try weights.barycentricFormulaConstants(eval_point);
    const verifier_query = ipa.VerifierQuery{
        .commitment = commitment,
        .B = b_verifier,
        .eval_point = eval_point,
        .result = ipa_proof.result,
        .proof = ipa_proof.proof,
    };
    const ok = ipa.verifyProof(xcrs, &verifier_transcript, verifier_query);
    try std.testing.expect(ok);
}
