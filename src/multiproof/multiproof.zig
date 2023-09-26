const std = @import("std");
const Allocator = std.mem.Allocator;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const Fr = banderwagon.Fr;
const lagrange_basis = @import("../polynomial/lagrange_basis.zig");
const Transcript = @import("../ipa/transcript.zig");
const ipa = @import("../ipa/ipa.zig");
const crs = @import("../crs/crs.zig");
const CRS = crs.CRS;
const msm = @import("../msm/msm.zig");
const precomputed_weights = @import("../polynomial/precomputed_weights.zig");

const IPA = ipa.IPA(crs.DomainSize);
const LagrangeBasis = lagrange_basis.LagrangeBasis(crs.DomainSize, crs.Domain);
const PrecomputedWeights = precomputed_weights.PrecomputedWeights(crs.DomainSize, crs.Domain);

// TODO: maybe group (potentially) multiple openings for the same poly?
pub const ProverQuery = struct {
    f: LagrangeBasis,
    C: Element,
    z: u8,
    y: Fr,
};

pub const VerifierQuery = struct {
    C: Element,
    z: u8,
    y: Fr,
};

pub const Proof = struct {
    ipa: IPA.IPAProof,
    D: Element,
};

pub const MultiProof = struct {
    precomp: PrecomputedWeights,
    crs: CRS,

    pub fn init(vkt_crs: CRS) !MultiProof {
        return MultiProof{
            .precomp = try PrecomputedWeights.init(),
            .crs = vkt_crs,
        };
    }

    pub fn createProof(
        self: MultiProof,
        transcript: *Transcript,
        queries: []const ProverQuery,
    ) !Proof {
        transcript.domainSep("multiproof");

        // Add queries into transcript
        for (queries) |query| {
            transcript.appendPoint(query.C, "C");
            transcript.appendScalar(Fr.fromInteger(@as(u256, query.z)), "z");
            transcript.appendScalar(query.y, "y");
        }

        // Generate challenge from queries
        const r = transcript.challengeScalar("r");

        var g = [_]Fr{Fr.zero()} ** crs.DomainSize;
        var power_of_r = Fr.one();
        for (queries) |query| {
            const f = query.f;
            const index = query.z;
            const quotient = self.computeQuotientInsideDomain(f, index);
            for (0..crs.DomainSize) |i| {
                g[i] = Fr.add(g[i], Fr.mul(power_of_r, quotient.evaluations[i]));
            }

            power_of_r = Fr.mul(power_of_r, r);
        }

        const D = try CRS.commit(self.crs, &g);
        transcript.appendPoint(D, "D");

        // Step 2: Compute h in evaluation form

        const t = transcript.challengeScalar("t");

        var h = [_]Fr{Fr.zero()} ** crs.DomainSize;
        power_of_r = Fr.one();

        for (queries) |query| {
            const f = query.f;
            const index = query.z;
            const denominator_inv = Fr.sub(t, crs.Domain[index]).inv().?;
            for (0..crs.DomainSize) |i| {
                h[i] = Fr.add(
                    h[i],
                    Fr.mul(Fr.mul(power_of_r, f.evaluations[i]), denominator_inv),
                );
            }

            power_of_r = Fr.mul(power_of_r, r);
        }

        var h_minus_g: [256]Fr = undefined;
        for (0..256) |i| {
            h_minus_g[i] = Fr.sub(h[i], g[i]);
        }

        // Step 3: Evaluate and compute IPA proofs

        const E = try CRS.commit(self.crs, &h);
        transcript.appendPoint(E, "E");

        var ipa_commitment: Element = undefined;
        ipa_commitment.sub(E, D);

        const polynomial = h_minus_g;
        const eval_point = t;
        const input_point_vector = try self.precomp.barycentricFormulaConstants(eval_point);

        var query = IPA.ProverQuery{
            .commitment = ipa_commitment,
            .A = polynomial,
            .B = input_point_vector,
            .eval_point = eval_point,
        };
        const proof_res = IPA.createProof(self.crs, transcript, query);

        return Proof{ .ipa = proof_res.proof, .D = D };
    }

    pub fn verifyProof(
        self: MultiProof,
        base_allocator: Allocator,
        transcript: *Transcript,
        queries: []const VerifierQuery,
        proof: Proof,
    ) !bool {
        if (queries.len == 0) {
            return error.QueriesCountIsZero;
        }

        var arena = std.heap.ArenaAllocator.init(base_allocator);
        defer arena.deinit();
        var allocator = arena.allocator();

        transcript.domainSep("multiproof");
        const D = proof.D;
        const ipa_proof = proof.ipa;

        for (queries) |query| {
            const C_i = query.C;
            const z_i = query.z;
            const y_i = query.y;
            transcript.appendPoint(C_i, "C");
            transcript.appendScalar(Fr.fromInteger(z_i), "z");
            transcript.appendScalar(y_i, "y");
        }

        // Compute powers of r.
        const r = transcript.challengeScalar("r");
        const powers_of_r = try powersOf(allocator, r, queries.len);

        // Group by eval point.
        var grouped_evals: [crs.DomainSize]Fr = grouped_evals: {
            var ret: [crs.DomainSize]Fr = undefined;
            for (0..crs.DomainSize) |i| {
                ret[i] = Fr.zero();
            }
            for (0..queries.len) |i| {
                const z = queries[i].z;
                const numerator = Fr.mul(powers_of_r[i], queries[i].y);
                ret[z] = Fr.add(ret[z], numerator);
            }
            break :grouped_evals ret;
        };

        transcript.appendPoint(D, "D");
        const t = transcript.challengeScalar("t");

        // Compute helpers 1/(t-z_i).
        const helper_scalar_den: [crs.DomainSize]Fr = helpers: {
            var den: [crs.DomainSize]Fr = undefined;
            for (0..crs.DomainSize) |z| {
                den[z] = Fr.sub(t, Fr.fromInteger(z));
            }
            var invDen: [crs.DomainSize]Fr = undefined;
            try Fr.batchInv(&invDen, &den);

            break :helpers invDen;
        };

        // Compute g_2(t) = sum(grouped_evals / helpers)
        var g_2_of_t = Fr.zero();
        for (0..crs.DomainSize) |i| {
            if (grouped_evals[i].isZero()) {
                continue;
            }
            g_2_of_t = Fr.add(g_2_of_t, Fr.mul(grouped_evals[i], helper_scalar_den[i]));
        }

        // Compute E = sum(C_i * r^i/(t-z_i))
        var E_coefficients = try allocator.alloc(Fr, queries.len);
        var Cs = try allocator.alloc(Element, queries.len);
        for (queries, 0..) |query, i| {
            Cs[i] = query.C;
            E_coefficients[i] = Fr.mul(powers_of_r[i], helper_scalar_den[queries[i].z]);
        }
        var precomp = try msm.PrecompMSM(8, 8).init(allocator, Cs);
        const E = try precomp.msm(E_coefficients);
        transcript.appendPoint(E, "E");

        // Check IPA proof.
        var E_minus_D: Element = undefined;
        E_minus_D.sub(E, D);

        const query = IPA.VerifierQuery{
            .commitment = E_minus_D,
            .B = try self.precomp.barycentricFormulaConstants(t),
            .eval_point = t,
            .result = g_2_of_t,
            .proof = ipa_proof,
        };
        return IPA.verifyProof(self.crs, transcript, query);
    }

    fn powersOf(allocator: Allocator, x: Fr, degree_minus_one: usize) ![]const Fr {
        var res = try allocator.alloc(Fr, degree_minus_one);
        res[0] = Fr.one();
        for (1..degree_minus_one) |i| {
            res[i] = Fr.mul(res[i - 1], x);
        }
        return res;
    }

    fn computeQuotientInsideDomain(self: MultiProof, f: LagrangeBasis, index: u8) LagrangeBasis {
        const inverses = self.precomp.domain_inverses;
        const Aprime_domain = self.precomp.Aprime_DOMAIN;
        const Aprime_domain_inv = self.precomp.Aprime_DOMAIN_inv;

        const indexIsize = @as(isize, index);

        var q = [_]Fr{Fr.zero()} ** crs.DomainSize;
        const y = f.evaluations[@as(usize, @intCast(indexIsize))];
        for (0..crs.DomainSize) |i| {
            if (i != indexIsize) {
                const den = @as(isize, @intCast(inverses.len));
                var num = @as(isize, @intCast(i)) - indexIsize;
                var inv_idx = @mod(num, den);
                q[i] = Fr.mul(Fr.sub(f.evaluations[i], y), inverses[@as(usize, @intCast(inv_idx))]);

                inv_idx = @mod(-inv_idx, den);
                q[@as(u8, @intCast(indexIsize))] = Fr.add(
                    q[@as(u8, @intCast(indexIsize))],
                    Fr.mul(Fr.mul(Fr.mul(Fr.sub(f.evaluations[i], y), inverses[@as(usize, @intCast(inv_idx))]), Aprime_domain[@as(u8, @intCast(indexIsize))]), Aprime_domain_inv[i]),
                );
            }
        }

        return LagrangeBasis.init(q);
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

    const xcrs = try CRS.init(std.testing.allocator);
    defer xcrs.deinit();
    const C_a = try CRS.commit(xcrs, &poly_eval_a);
    const C_b = try CRS.commit(xcrs, &poly_eval_b);
    const zs = [_]u8{ 0, 0 };
    const ys = [_]Fr{ Fr.fromInteger(1), Fr.fromInteger(32) };
    const fs = [_][256]Fr{ poly_eval_a, poly_eval_b };
    const Cs = [_]Element{ C_a, C_b };

    const query_a = ProverQuery{
        .f = LagrangeBasis.init(fs[0]),
        .C = Cs[0],
        .z = zs[0],
        .y = ys[0],
    };
    const query_b = ProverQuery{
        .f = LagrangeBasis.init(fs[1]),
        .C = Cs[1],
        .z = zs[1],
        .y = ys[1],
    };

    const multiproof = try MultiProof.init(xcrs);

    var prover_transcript = Transcript.init("test");
    const proof = try multiproof.createProof(&prover_transcript, &[_]ProverQuery{ query_a, query_b });

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");

    try std.testing.expectEqualStrings(
        "eee8a80357ff74b766eba39db90797d022e8d6dee426ded71234241be504d519",
        &std.fmt.bytesToHex(p_challenge.toBytes(), std.fmt.Case.lower),
    );

    var verifier_transcript = Transcript.init("test");
    var vquery_a = VerifierQuery{ .C = Cs[0], .z = zs[0], .y = ys[0] };
    var vquery_b = VerifierQuery{ .C = Cs[1], .z = zs[1], .y = ys[1] };
    const ok = try multiproof.verifyProof(allocator, &verifier_transcript, &[_]VerifierQuery{ vquery_a, vquery_b }, proof);

    try std.testing.expect(ok);

    const v_challenge = verifier_transcript.challengeScalar("state");
    try std.testing.expectEqualSlices(u8, &p_challenge.toBytes(), &v_challenge.toBytes());
}
