const std = @import("std");
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Banderwagon = @import("../ecc/bandersnatch/banderwagon.zig").Banderwagon;
const Fr = Bandersnatch.Fr;
const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights;
const CRS = @import("../crs/crs.zig").CRS;
const Transcript = @import("transcript.zig");
const Common = @import("common.zig");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;

pub const ProverQuery = struct {
    polynomial: []Fr,
    commitment: Banderwagon,
    // Input point
    point: Fr,
    // If polynomial was in monomial basis
    // this would be <1, b, b^2, b^3, b^4,..., b^n>
    // TODO: Can we give this a better name?
    point_evaluations: []Fr,
};

pub const Proof = struct {
    L: ArrayList(Banderwagon),
    R: ArrayList(Banderwagon),
    a: Fr,

    pub fn deinit(self: *Proof) void {
        self.L.deinit();
        self.R.deinit();
    }
};

pub const VerifierQuery = struct {
    commitment: Banderwagon,
    point: Fr,
    // If polynomial was in monomial basis
    // this would be <1, b, b^2, b^3, b^4,..., b^n>
    point_evaluations: []Fr,
    output_point: Fr,
    proof: Proof,
};

// Proves that <a, b> = inner_product
// <a,b> is read as a inner product with b
// Furthermore, b is assumed to be public, hence we will not commit to it
//
// The caller is responsible for calling deinit() in the returned proof.
pub fn make_ipa_proof(allocator: Allocator, crs: CRS, transcript: *Transcript, query: ProverQuery) !struct { result: Fr, proof: Proof } {
    transcript.domainSep("ipa");

    var n = query.polynomial.len;
    var m = n / 2;

    var a = query.polynomial;
    var b = query.point_evaluations;
    assert(a.len == b.len);
    // TODO(jsign): some extra assertions woudn't hurt here.

    const y = Common.innerProduct(a, b);

    var proof = Proof{
        .L = try ArrayList(Banderwagon).initCapacity(allocator, n),
        .R = try ArrayList(Banderwagon).initCapacity(allocator, n),
        .a = Fr.zero(),
    };

    transcript.appendPoint(query.commitment, "C");
    transcript.appendScalar(query.point, "input point");
    transcript.appendScalar(y, "output point");
    const w = transcript.challengeScalar("w");

    const q = crs.BASIS_Q.scalarMul(w);

    var current_basis = crs.BASIS_G;

    while (n > 1) {
        // Reduction step

        const a_L = a[0..m];
        const a_R = a[m..];
        const b_L = b[0..m];
        const b_R = b[m..];
        const z_L = Common.innerProduct(a_R, b_L);
        const z_R = Common.innerProduct(a_L, b_R);

        var C_L: Banderwagon = undefined;
        var C_R: Banderwagon = undefined;
        C_L.add(&Banderwagon.msm(current_basis[0..m], a_R), &q.scalarMul(z_L));
        C_R.add(&Banderwagon.msm(current_basis[m..], a_L), &q.scalarMul(z_R));

        try proof.L.append(C_L);
        try proof.R.append(C_R);

        transcript.appendPoint(C_L, "L");
        transcript.appendPoint(C_R, "R");
        const x = transcript.challengeScalar("x");

        var xinv: Fr = x.inv().?;

        // Compute updates for next round
        for (a_L, a_R, 0..) |v, z, i| {
            a[i] = v.add(x.mul(z));
        }
        a = a[0..m];

        for (b_L, b_R, 0..) |v, z, i| {
            b[i] = v.add(xinv.mul(z));
        }
        b = b[0..m];

        for (current_basis[0..m], current_basis[m..], 0..) |*v, *z, i| {
            current_basis[i].add(v, &z.scalarMul(xinv));
        }
        current_basis = current_basis[0..m];

        n = m;
        m = n / 2;
    }

    proof.a = a[0];

    return .{ .result = y, .proof = proof };
}

// def check_ipa_proof(crs: CRS, transcript: Transcript, query: VerifierQuery):
//     transcript.domain_sep(b"ipa")
//     # TODO: We should add `n` into the transcript.
//     # TODO: this breaks compatibility, with other implementations
//     # TODO: so lets wait until reference is completed
//     n = len(query.point_evaluations)
//     m = n // 2

//     C = query.commitment
//     z = query.point
//     b = query.point_evaluations
//     proof = query.proof
//     y = query.output_point

//     transcript.append_point(C, b"C")
//     transcript.append_scalar(z, b"input point")
//     transcript.append_scalar(y, b"output point")
//     w = transcript.challenge_scalar(b"w")

//     q = crs.BASIS_Q * w

//     current_commitment = C + (q * y)

//     i = 0
//     xs = []
//     xinvs = []

//     while n > 1:
//         C_L = proof.L[i]
//         C_R = proof.R[i]
//         transcript.append_point(C_L, b"L")
//         transcript.append_point(C_R, b"R")
//         x = transcript.challenge_scalar(b"x")

//         x_inv = Fr.zero()
//         xinv = x_inv.inv(x)

//         xs.append(x)
//         xinvs.append(xinv)

//         current_commitment = current_commitment + (C_L * x) + (C_R * xinv)

//         n = m
//         m = n // 2
//         i = i + 1

//     # Do it the inefficient way
//     current_basis = crs.BASIS_G

//     for i in range(len(xs)):

//         G_L, G_R = split_points(current_basis)
//         b_L, b_R = split_scalars(b)

//         x_inv = xinvs[i]

//         b = fold_scalars(b_L, b_R, x_inv)
//         current_basis = fold_points(G_L, G_R, x_inv)

//     assert len(b) == len(current_basis)
//     assert len(b) == 1

//     b_0 = b[0]
//     G_0 = current_basis[0]

//     # G[0] * a + (a * b) * Q;
//     got_commitment = G_0 * proof.a + q * (proof.a * b_0)

//     return current_commitment == got_commitment

// # Computes c[i] = a[i] + b[i] * challenge
// def fold_list(a, b, folding_challenge: Fr):
//     assert len(a) == len(b)
//     result = []
//     for a_i, b_i in zip(a, b):
//         result.append(a_i + b_i * folding_challenge)
//     return result

// def fold_points(a: List[Banderwagon], b: List[Banderwagon], folding_challenge: Fr):
//     return fold_list(a, b, folding_challenge)

// def fold_scalars(a: List[Fr], b: List[Fr], folding_challenge: Fr):
//     return fold_list(a, b, folding_challenge)

// #  Splits a list into two lists of equal length
// #  Eg[S1, S2, S3, S4] becomes[S1, S2], [S3, S4]
// def split_list_in_half(x):
//     assert len(x) % 2 == 0

//     mid = len(x) // 2
//     return (x[:mid], x[mid:])

// def split_scalars(x: List[Fr]):
//     return split_list_in_half(x)

// def split_points(x: List[Banderwagon]):
//     return split_list_in_half(x)

const testAllocator = std.testing.allocator;
test "basic proof" {
    // Test a simple IPA proof
    var domain = try testAllocator.alloc(Fr, 256);
    defer testAllocator.free(domain);
    for (0..256) |i| {
        domain[i] = Fr.fromInteger(i);
    }
    var weights = try PrecomputedWeights.init(testAllocator, domain);
    defer weights.deinit();

    // Polynomial in lagrange basis
    var lagrange_poly = try testAllocator.alloc(Fr, 256);
    defer testAllocator.free(lagrange_poly);
    for (0..256) |i| {
        lagrange_poly[i] = Fr.fromInteger(i % 32 + 1);
    }

    // Commit to the polynomial in lagrange basis
    const crs = try CRS.init(testAllocator);
    defer crs.deinit();
    const commitment = crs.commit(lagrange_poly);

    const expected_comm = std.fmt.bytesToHex(commitment.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("1b9dff8f5ebbac250d291dfe90e36283a227c64b113c37f1bfb9e7a743cdb128", &expected_comm);

    var prover_transcript = Transcript.init("test");

    // create a opening proof for a point outside of the domain
    const input_point = Fr.fromInteger(2101);
    const b = try weights.barycentricFormulaConstants(testAllocator, input_point);
    defer testAllocator.free(b);
    const output_point_check = Common.innerProduct(lagrange_poly, b);
    const output_point_check_hex = std.fmt.bytesToHex(output_point_check.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("4a353e70b03c89f161de002e8713beec0d740a5e20722fd5bd68b30540a33208", &output_point_check_hex);

    var query = ProverQuery{
        .polynomial = lagrange_poly,
        .commitment = commitment,
        .point = input_point,
        .point_evaluations = b,
    };

    var ipa_proof = try make_ipa_proof(testAllocator, crs, &prover_transcript, query);
    defer ipa_proof.proof.deinit();

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");

    const p_challenge_hex = std.fmt.bytesToHex(p_challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("0a81881cbfd7d7197a54ebd67ed6a68b5867f3c783706675b34ece43e85e7306", &p_challenge_hex);

    //         verifier_transcript = Transcript(b"test")

    //         query = VerifierQuery(commitment, input_point, b, output_point, proof)

    //         ok = check_ipa_proof(crs, verifier_transcript, query)

    //         self.assertTrue(ok)
}
