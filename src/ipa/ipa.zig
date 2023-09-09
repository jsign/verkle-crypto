const std = @import("std");
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Banderwagon = @import("../ecc/bandersnatch/banderwagon.zig").Banderwagon;
const Fr = Bandersnatch.Fr;
const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights;
const CRS = @import("../crs/crs.zig").CRS;

// @dataclass
// class ProverQuery:
//     polynomial: List[Fr]

//     commitment: Banderwagon
//     # Input point
//     point: Fr
//     # If polynomial was in monomial basis
//     # this would be <1, b, b^2, b^3, b^4,..., b^n>
//     # TODO: Can we give this a better name?
//     point_evaluations: List[Fr]

// @dataclass
// class Proof:
//     L: List[Banderwagon]
//     R: List[Banderwagon]
//     a: Fr

// @dataclass
// class VerifierQuery:
//     commitment: Banderwagon

//     point: Fr
//     # If polynomial was in monomial basis
//     # this would be <1, b, b^2, b^3, b^4,..., b^n>
//     point_evaluations: List[Fr]

//     output_point: Fr
//     # TODO: change this from typing.Any
//     proof: Proof

// # Proves that <a, b> = inner_product
// # <a,b> is read as a inner product with b
// # Furthermore, b is assumed to be public, hence we will not commit to it
// def make_ipa_proof(crs: CRS, transcript: Transcript, query: ProverQuery):
//     transcript.domain_sep(b"ipa")

//     n = len(query.polynomial)
//     m = n // 2

//     a = query.polynomial
//     b = query.point_evaluations
//     y = inner_product(a, b)

//     proof = Proof([], [], Fr.zero())

//     transcript.append_point(query.commitment, b"C")
//     transcript.append_scalar(query.point, b"input point")
//     transcript.append_scalar(y, b"output point")
//     w = transcript.challenge_scalar(b"w")

//     q = crs.BASIS_Q * w

//     current_basis = crs.BASIS_G

//     while n > 1:
//         # Reduction step

//         a_L = a[:m]
//         a_R = a[m:]
//         b_L = b[:m]
//         b_R = b[m:]
//         z_L = inner_product(a_R, b_L)
//         z_R = inner_product(a_L, b_R)

//         C_L: Banderwagon = varbase_commit(a_R, current_basis[:m]) + (q * z_L)
//         C_R: Banderwagon = varbase_commit(a_L, current_basis[m:]) + (q * z_R)

//         proof.L.append(C_L)
//         proof.R.append(C_R)

//         transcript.append_point(C_L, b"L")
//         transcript.append_point(C_R, b"R")
//         x = transcript.challenge_scalar(b"x")

//         xinv = Fr.zero()
//         xinv.inv(x)

//         # Compute updates for next round
//         a = [(v + x * w) for v, w in zip(a_L, a_R)]
//         b = [(v + xinv * w) for v, w in zip(b_L, b_R)]

//         current_basis = [v + (w * xinv)
//                          for v, w in zip(current_basis[:m], current_basis[m:])]
//         n = m
//         m = n // 2

//     proof.a = a[0]

//     return y, proof

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

    //         prover_transcript = Transcript(b"test")

    //         # create a opening proof for a point outside of the domain
    //         input_point = Fr(2101)
    //         b = weights.barycentric_formula_constants(input_point)
    //         query = ProverQuery(lagrange_poly, commitment, input_point, b)

    //         output_point, proof = make_ipa_proof(crs, prover_transcript, query)

    //         # # Lets check the state of the transcript by squeezing out another challenge
    //         p_challenge = prover_transcript.challenge_scalar(b"state")

    //         self.assertEqual(
    //             "0a81881cbfd7d7197a54ebd67ed6a68b5867f3c783706675b34ece43e85e7306", p_challenge.to_bytes().hex())

    //         verifier_transcript = Transcript(b"test")

    //         query = VerifierQuery(commitment, input_point, b, output_point, proof)

    //         ok = check_ipa_proof(crs, verifier_transcript, query)

    //         self.assertTrue(ok)
}
