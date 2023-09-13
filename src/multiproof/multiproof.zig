const std = @import("std");
const Allocator = std.mem.Allocator;
const Fr = @import("../ecc/bandersnatch/bandersnatch.zig").Fr;
const Banderwagon = @import("../ecc/bandersnatch/banderwagon.zig").Banderwagon;
const LagrangeBasis = @import("../polynomial/lagrange_basis.zig").LagrangeBasis;
const Transcript = @import("../ipa/transcript.zig");
const IPA = @import("../ipa/ipa.zig");
const CRS = @import("../crs/crs.zig").CRS;
const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights;
const Quotient = @import("quotient.zig");

//TODO: FIX
fn varbase_commit(values: []const Fr, elements: []const Banderwagon) Banderwagon {
    return Banderwagon.msm(elements, values);
}

// We could store the polynomial once and just gather the queries for
// that polynomial. This way is in-efficient, however it's easier to read
const ProverQuery = struct {
    f: LagrangeBasis,
    C: Banderwagon,
    z: Fr,
    y: Fr,
};

const VerifierQuery = struct {
    C: Banderwagon,
    z: Fr,
    y: Fr,
};

const Proof = struct {
    ipa: IPA.Proof,
    D: Banderwagon,
};

const MultiProof = struct {
    precomp: PrecomputedWeights,
    crs: CRS,

    fn init(allocator: Allocator, domain: []Fr, crs: CRS) !MultiProof {
        return MultiProof{
            .precomp = try PrecomputedWeights.init(allocator, domain),
            .crs = crs, // TODO: check if assignement does copy (if needed)
        };
    }

    fn make_multiproof(
        self: MultiProof,
        allocator: Allocator,
        transcript: *Transcript,
        queries: []const ProverQuery,
    ) !Proof {
        const domain_size = self.precomp.domain.len; // TODO: comptime.

        transcript.domainSep("multiproof");

        // Add queries into transcript
        for (queries) |query| {
            transcript.appendPoint(query.C, "C");
            transcript.appendScalar(query.z, "z");
            transcript.appendScalar(query.y, "y");
        }

        // Generate challenge from queries
        const r = transcript.challengeScalar("r");

        var g = [_]Fr{Fr.zero()} ** 256; // TODO: comptime.
        var power_of_r = Fr.one();
        for (queries) |query| {
            const f = query.f;
            const index = query.z;
            const quotient = try Quotient.compute_quotient_inside_domain(allocator, self.precomp, f, index);
            for (0..domain_size) |i| {
                g[i] = Fr.add(g[i], Fr.mul(power_of_r, quotient.evaluations.items[i]));
            }

            power_of_r = Fr.mul(power_of_r, r);
        }

        const D = self.crs.commit(&g);
        transcript.appendPoint(D, "D");

        // Step 2: Compute h in evaluation form

        const t = transcript.challengeScalar("t");

        var h = [_]Fr{Fr.zero()} ** 256; // TODO
        power_of_r = Fr.one();

        for (queries) |query| {
            const f = query.f;
            const index = @as(u8, @intCast(query.z.toInteger()));
            const denominator_inv = Fr.sub(t, self.precomp.domain[index]).inv().?;
            for (0..domain_size) |i| {
                h[i] = Fr.add(
                    h[i],
                    Fr.mul(Fr.mul(power_of_r, f.evaluations.items[i]), denominator_inv),
                );
            }

            power_of_r = Fr.mul(power_of_r, r);
        }

        var h_minus_g: [256]Fr = undefined;
        for (0..256) |i| {
            h_minus_g[i] = Fr.sub(h[i], g[i]);
        }

        // Step 3: Evaluate and compute IPA proofs

        const E = self.crs.commit(&h);
        transcript.appendPoint(E, "E");

        var ipa_commitment: Banderwagon = undefined;
        ipa_commitment.sub(E, D);

        const polynomial = h_minus_g;
        const input_point = t;
        const input_point_vector = try self.precomp.barycentricFormulaConstants(allocator, input_point);

        var query = IPA.ProverQuery{
            .polynomial = &polynomial,
            .commitment = ipa_commitment,
            .point = input_point,
            .point_evaluations = input_point_vector,
        };

        // TODO: simplify proof_res
        const proof_res = try IPA.make_ipa_proof(allocator, self.crs, transcript, &query);

        return Proof{ .ipa = proof_res.proof, .D = D };
    }

    fn check_multiproof(
        self: MultiProof,
        allocator: Allocator,
        transcript: Transcript,
        queries: []VerifierQuery,
        proof: Proof,
    ) !bool {
        transcript.domain_sep("multiproof");
        const D = proof.D;
        const ipa_proof = proof.ipa;

        for (queries) |query| {
            const C_i = query.C;
            const z_i = query.z;
            const y_i = query.y;
            transcript.append_point(C_i, "C");
            transcript.append_scalar(z_i, "z");
            transcript.append_scalar(y_i, "y");
        }

        // Step 1
        const r = transcript.challenge_scalar("r");

        // Step 2
        transcript.append_point(D, "D");
        const t = transcript.challenge_scalar("t");

        var g_2_of_t = Fr.zero();
        var power_of_r = Fr.one();

        var E_coefficients = try allocator.alloc(Fr, queries.len);
        var Cs = try allocator.alloc(Banderwagon, queries.len);
        for (queries, 0..) |query, i| {
            Cs[i] = query.C;
            const z = query.z.toInteger();
            const y = query.y;
            E_coefficients[i] = Fr.mul(power_of_r, Fr.sub(t, self.precomp.domain[z]).inv().?);
            g_2_of_t = Fr.add(g_2_of_t, Fr.mul(E_coefficients[i], y));

            power_of_r = Fr.mul(power_of_r, r);
        }

        const E = varbase_commit(E_coefficients, Cs);
        transcript.append_point(E, "E");

        // Step 3 (Check IPA proofs)
        const y = g_2_of_t;
        const ipa_commitment = E - D;
        const input_point = t;
        const output_point = y;
        const input_point_vector = self.precomp.barycentricFormulaConstants(allocator, input_point);

        const query = IPA.VerifierQuery(ipa_commitment, input_point, input_point_vector, output_point, ipa_proof);
        return IPA.check_ipa_proof(allocator, self.crs, transcript, query);
    }
};

test "basic" {
    var allocator = std.testing.allocator;

    // Polynomials in lagrange basis
    const poly_eval_a = [_]Fr{
        Fr.fromInteger(1),
        Fr.fromInteger(2),
        Fr.fromInteger(3),
        Fr.fromInteger(4),
        Fr.fromInteger(5),
        Fr.fromInteger(6),
        Fr.fromInteger(7),
        Fr.fromInteger(8),
        Fr.fromInteger(9),
        Fr.fromInteger(10),
        Fr.fromInteger(11),
        Fr.fromInteger(12),
        Fr.fromInteger(13),
        Fr.fromInteger(14),
        Fr.fromInteger(15),
        Fr.fromInteger(16),
        Fr.fromInteger(17),
        Fr.fromInteger(18),
        Fr.fromInteger(19),
        Fr.fromInteger(20),
        Fr.fromInteger(21),
        Fr.fromInteger(22),
        Fr.fromInteger(23),
        Fr.fromInteger(24),
        Fr.fromInteger(25),
        Fr.fromInteger(26),
        Fr.fromInteger(27),
        Fr.fromInteger(28),
        Fr.fromInteger(29),
        Fr.fromInteger(30),
        Fr.fromInteger(31),
        Fr.fromInteger(32),
    } ** 8;
    const poly_eval_b = [_]Fr{
        Fr.fromInteger(32),
        Fr.fromInteger(31),
        Fr.fromInteger(30),
        Fr.fromInteger(29),
        Fr.fromInteger(28),
        Fr.fromInteger(27),
        Fr.fromInteger(26),
        Fr.fromInteger(25),
        Fr.fromInteger(24),
        Fr.fromInteger(23),
        Fr.fromInteger(22),
        Fr.fromInteger(21),
        Fr.fromInteger(20),
        Fr.fromInteger(19),
        Fr.fromInteger(18),
        Fr.fromInteger(17),
        Fr.fromInteger(16),
        Fr.fromInteger(15),
        Fr.fromInteger(14),
        Fr.fromInteger(13),
        Fr.fromInteger(12),
        Fr.fromInteger(11),
        Fr.fromInteger(10),
        Fr.fromInteger(9),
        Fr.fromInteger(8),
        Fr.fromInteger(7),
        Fr.fromInteger(6),
        Fr.fromInteger(5),
        Fr.fromInteger(4),
        Fr.fromInteger(3),
        Fr.fromInteger(2),
        Fr.fromInteger(1),
    } ** 8;

    const crs = try CRS.init(allocator);
    const C_a = crs.commit(&poly_eval_a);
    const C_b = crs.commit(&poly_eval_b);
    const zs = [_]Fr{ Fr.zero(), Fr.zero() };
    const ys = [_]Fr{ Fr.fromInteger(1), Fr.fromInteger(32) };
    const fs = [_][256]Fr{ poly_eval_a, poly_eval_b };
    const Cs = [_]Banderwagon{ C_a, C_b };

    var domain: [256]Fr = undefined;
    for (0..domain.len) |i| {
        domain[i] = Fr.fromInteger(i);
    }

    const query_a = ProverQuery{
        .f = try LagrangeBasis.init(allocator, &fs[0], &domain),
        .C = Cs[0],
        .z = zs[0],
        .y = ys[0],
    };
    const query_b = ProverQuery{
        .f = try LagrangeBasis.init(allocator, &fs[1], &domain),
        .C = Cs[1],
        .z = zs[1],
        .y = ys[1],
    };

    const multiproof = try MultiProof.init(allocator, &domain, crs);

    var prover_transcript = Transcript.init("test");
    const proof = try multiproof.make_multiproof(allocator, &prover_transcript, &[_]ProverQuery{ query_a, query_b });
    _ = proof;

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");

    try std.testing.expectEqualStrings(
        "eee8a80357ff74b766eba39db90797d022e8d6dee426ded71234241be504d519",
        &std.fmt.bytesToHex(p_challenge.to_bytes(), std.fmt.Case.lower),
    );

    //         verifier_transcript = Transcript(b"test")
    //         query_a = VerifierQuery(Cs[0], zs[0], ys[0])
    //         query_b = VerifierQuery(Cs[1], zs[1], ys[1])
    //         ok = multiproof.check_multiproof(
    //             verifier_transcript, [query_a, query_b], proof)

    //         self.assertTrue(ok)

    //         v_challenge = verifier_transcript.challenge_scalar(b"state")
    //         self.assertEqual(v_challenge, p_challenge)
}
