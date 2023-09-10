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
pub fn make_ipa_proof(allocator: Allocator, crs: CRS, transcript: *Transcript, query: *ProverQuery) !struct { result: Fr, proof: Proof } {
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

    var current_basis = try allocator.dupe(Banderwagon, crs.BASIS_G);

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

pub fn check_ipa_proof(base_allocator: Allocator, crs: *const CRS, transcript: *Transcript, query: *const VerifierQuery) !bool {
    var arena = std.heap.ArenaAllocator.init(base_allocator);
    defer arena.deinit();
    var allocator = arena.allocator();

    transcript.domainSep("ipa");
    // TODO: We should add `n` into the transcript.
    // TODO: this breaks compatibility, with other implementations
    // TODO: so lets wait until reference is completed
    var n = query.point_evaluations.len;
    var m = n / 2;

    const C = query.commitment;
    const z = query.point;
    var b = query.point_evaluations;
    const proof = query.proof;
    const y = query.output_point;

    transcript.appendPoint(C, "C");
    transcript.appendScalar(z, "input point");
    transcript.appendScalar(y, "output point");
    const w = transcript.challengeScalar("w");

    const q = crs.BASIS_Q.scalarMul(w);

    var current_commitment: Banderwagon = undefined;
    current_commitment.add(&C, &q.scalarMul(y));

    var i: u8 = 0;
    var xs = ArrayList(Fr).init(allocator);
    var xinvs = ArrayList(Fr).init(allocator);

    while (n > 1) {
        const C_L = proof.L.items[i];
        const C_R = proof.R.items[i];
        transcript.appendPoint(C_L, "L");
        transcript.appendPoint(C_R, "R");

        const x = transcript.challengeScalar("x");
        const x_inv = x.inv().?;

        try xs.append(x);
        try xinvs.append(x_inv);

        var tmp: Banderwagon = undefined;
        tmp.add(&C_L.scalarMul(x), &C_R.scalarMul(x_inv));
        current_commitment.add(&current_commitment, &tmp);

        n = m;
        m = n / 2;
        i = i + 1;
    }

    // Do it the inefficient way (TODO: optimize)
    var current_basis = crs.BASIS_G;

    for (0..xs.items.len) |j| {
        const g_split = split_list_in_half(Banderwagon, current_basis);
        const b_split = split_list_in_half(Fr, b);

        const x_inv = xinvs.items[j];

        b = try fold_scalars(allocator, b_split.L, b_split.R, x_inv);
        current_basis = try fold_points(allocator, g_split.L, g_split.R, x_inv);
    }

    assert(b.len == current_basis.len);
    assert(b.len == 1);

    const b_0 = b[0];
    const G_0 = current_basis[0];

    // G[0] * a + (a * b) * Q;
    var got_commitment: Banderwagon = undefined;
    got_commitment.add(&G_0.scalarMul(proof.a), &q.scalarMul(Fr.mul(proof.a, b_0)));

    return current_commitment.eq(&got_commitment);
}

// Computes c[i] = a[i] + b[i] * challenge
fn fold_scalars(allocator: Allocator, a: []const Fr, b: []const Fr, folding_challenge: Fr) ![]Fr {
    assert(a.len == b.len);
    var result = try allocator.alloc(Fr, a.len);
    for (a, b, 0..) |a_i, b_i, i| {
        result[i] = a_i.add(b_i.mul(folding_challenge));
    }
    return result;
}

fn fold_points(allocator: Allocator, a: []const Banderwagon, b: []const Banderwagon, folding_challenge: Fr) ![]Banderwagon {
    assert(a.len == b.len);
    var result = try allocator.alloc(Banderwagon, a.len);
    for (a, b, 0..) |a_i, b_i, i| {
        result[i].add(&a_i, &b_i.scalarMul(folding_challenge));
    }
    return result;
}

// Splits a list into two lists of equal length
// Eg[S1, S2, S3, S4] becomes[S1, S2], [S3, S4]
fn split_list_in_half(comptime T: type, x: []const T) struct { L: []const T, R: []const T } {
    assert(x.len % 2 == 0);
    const mid = x.len / 2;
    return .{ .L = x[0..mid], .R = x[mid..] };
}

const test_allocator = std.testing.allocator;
test "basic proof" {
    // TODO: use base allocator.
    var arena = std.heap.ArenaAllocator.init(test_allocator);
    defer arena.deinit();
    var allocator = arena.allocator();

    // Test a simple IPA proof
    var domain = try allocator.alloc(Fr, 256);
    defer allocator.free(domain);
    for (0..256) |i| {
        domain[i] = Fr.fromInteger(i);
    }
    var weights = try PrecomputedWeights.init(allocator, domain);
    defer weights.deinit();

    // Polynomial in lagrange basis
    var lagrange_poly = try allocator.alloc(Fr, 256);
    defer allocator.free(lagrange_poly);
    for (0..256) |i| {
        lagrange_poly[i] = Fr.fromInteger(i % 32 + 1);
    }

    // Commit to the polynomial in lagrange basis
    const crs = try CRS.init(allocator);
    defer crs.deinit();
    const commitment = crs.commit(lagrange_poly);

    const expected_comm = std.fmt.bytesToHex(commitment.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("1b9dff8f5ebbac250d291dfe90e36283a227c64b113c37f1bfb9e7a743cdb128", &expected_comm);

    var prover_transcript = Transcript.init("test");

    // create a opening proof for a point outside of the domain
    const input_point = Fr.fromInteger(2101);
    var b = try weights.barycentricFormulaConstants(allocator, input_point);
    defer allocator.free(b);
    const output_point_check = Common.innerProduct(lagrange_poly, b);
    const output_point_check_hex = std.fmt.bytesToHex(output_point_check.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("4a353e70b03c89f161de002e8713beec0d740a5e20722fd5bd68b30540a33208", &output_point_check_hex);

    var query = ProverQuery{
        .polynomial = lagrange_poly,
        .commitment = commitment,
        .point = input_point,
        .point_evaluations = b,
    };

    var ipa_proof = try make_ipa_proof(allocator, crs, &prover_transcript, &query);
    defer ipa_proof.proof.deinit();

    // Lets check the state of the transcript by squeezing out another challenge
    const p_challenge = prover_transcript.challengeScalar("state");
    const p_challenge_hex = std.fmt.bytesToHex(p_challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("0a81881cbfd7d7197a54ebd67ed6a68b5867f3c783706675b34ece43e85e7306", &p_challenge_hex);

    // Verify the proof.
    var verifier_transcript = Transcript.init("test");
    const b_verifier = try weights.barycentricFormulaConstants(allocator, input_point);
    defer allocator.free(b_verifier);
    const verifier_query = VerifierQuery{
        .commitment = commitment,
        .point = input_point,
        .point_evaluations = b_verifier,
        .output_point = ipa_proof.result,
        .proof = ipa_proof.proof,
    };
    const ok = try check_ipa_proof(allocator, &crs, &verifier_transcript, &verifier_query);
    try std.testing.expect(ok);
}
