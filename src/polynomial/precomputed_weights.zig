const std = @import("std");
const assert = std.debug.assert;
const Allocator = std.mem.Allocator;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;
const MonomialBasis = @import("monomial_basis.zig").MonomialBasis;
const LagrangeBasis = @import("lagrange_basis.zig").LagrangeBasis; // TODO(jsign): reconsider having the wrapper.

// TODO(jsign): comptime?
pub const PrecomputedWeights = struct {
    // TODO: We don't _need_ to store the vanishing polynomial
    // TODO: we only need to store its derivative and whenever we need to evaluate
    // TODO: the vanishing polynomial, it can be done via the domain

    // Vanishing polynomial
    A: MonomialBasis,
    // Derivative of the vanishing polynomial
    Aprime: MonomialBasis,
    // Aprime evaluated on the domain
    Aprime_DOMAIN: []Fr,
    // Aprime evaluated on the domain and then inverted
    Aprime_DOMAIN_inv: []Fr,
    // Domain
    domain: []const Fr,
    //Inverse of the domain
    domain_inverses: []const Fr,

    // TODO(jsign): avoid allocators.
    allocator: Allocator,

    pub fn init(allocator: Allocator, domain: []Fr) !PrecomputedWeights {
        assert(checkDomainIsContinousAndIncreasing(domain));

        var dom = try allocator.alloc(Fr, domain.len);
        for (domain, 0..) |x, i| {
            dom[i] = x;
        }
        const A = try MonomialBasis.vanishingPoly(allocator, dom);
        const Aprime = try MonomialBasis.formalDerivative(allocator, A);
        const Aprime_domain = try allocator.alloc(Fr, domain.len);
        const Aprime_domain_inv = try allocator.alloc(Fr, domain.len);

        for (0..domain.len) |i| {
            Aprime_domain[i] = Aprime.evaluate(Fr.fromInteger(i));
            Aprime_domain_inv[i] = Fr.inv(Aprime_domain[i]).?; // TODO(jsign): could do batching.
        }

        // This is not fully correct as the first element will be the inverse of 0
        // We keep it this way for now because it is what the research code did
        // TODO: refactor this to make it more readable
        // If domain size is 4 for example, the output would be:
        // [1/0, 1/1, 1/2, 1/3, -1/3, -1/2,-1/1]
        const inverses = try allocator.alloc(Fr, 2 * domain.len);
        for (0..domain.len) |d| {
            inverses[d] = Fr.inv(Fr.fromInteger(d)) orelse Fr.zero(); // It should only happen for 0.
            inverses[domain.len + d] = Fr.inv(Fr.fromInteger(Fr.MODULO - d)) orelse Fr.zero();
        }

        return PrecomputedWeights{
            .allocator = allocator,
            .A = A,
            .Aprime = Aprime,
            .Aprime_DOMAIN = Aprime_domain,
            .Aprime_DOMAIN_inv = Aprime_domain_inv,
            .domain = dom,
            .domain_inverses = inverses,
        };
    }

    pub fn deinit(self: *PrecomputedWeights) void {
        self.A.deinit();
        self.Aprime.deinit();
        self.allocator.free(self.Aprime_DOMAIN);
        self.allocator.free(self.Aprime_DOMAIN_inv);
        self.allocator.free(self.domain);
        self.allocator.free(self.domain_inverses);
    }

    // barycentricFormularConstants returns a slice with the constants to be used when evaluating a polynomial at z.
    // b_i = A(z) / A'(DOMAIN[i]) * 1 / (z - DOMAIN[i])
    // The caller is responsible for freeing the returned slice.
    pub fn barycentricFormulaConstants(self: *PrecomputedWeights, alloc: Allocator, z: Fr) ![]Fr {
        const Az = self.A.evaluate(z);

        var inverses = try alloc.alloc(Fr, self.domain.len);
        defer alloc.free(inverses);
        for (self.domain, 0..) |x, i| {
            // TODO(jsign): batching.
            inverses[i] = Fr.inv(Fr.sub(z, x)).?;
        }

        var r = try self.allocator.alloc(Fr, self.domain.len);
        for (0..self.domain.len) |i| {
            r[i] = Fr.mul(Az, Fr.mul(self.Aprime_DOMAIN_inv[i], inverses[i]));
        }

        return r;
    }
};

fn checkDomainIsContinousAndIncreasing(domain: []Fr) bool {
    for (1..domain.len) |i| {
        if (!Fr.sub(domain[i], domain[i - 1]).isOne()) {
            return false;
        }
    }
    return true;
}

// TODO: this test is commented in the spec.
// #     def test_domain_correctness(self):
// #         """
// #             Test that the domain is continuos and increasing. If this is not the case
// #             then the setup for precomputed weights, will need to be changed
// #         """
// #         domain = [Fr(0), Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)]
// #         domain_is_correct = check_domain_is_continuous_and_increasing(domain)
// #         self.assertTrue(domain_is_correct)

// #         # This domain has a gap between 0 and 2
// #         domain = [Fr(0), Fr(2), Fr(3), Fr(4), Fr(5)]
// #         domain_is_correct = check_domain_is_continuous_and_increasing(domain)
// #         self.assertFalse(domain_is_correct)

// #         # This domain is not increasing
// #         domain = [Fr(5), Fr(4), Fr(3), Fr(2), Fr(1)]
// #         domain_is_correct = check_domain_is_continuous_and_increasing(domain)
// #         self.assertFalse(domain_is_correct)
