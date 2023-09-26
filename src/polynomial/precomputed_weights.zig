const std = @import("std");
const assert = std.debug.assert;
const crs = @import("../crs/crs.zig");
const Fr = @import("../banderwagon/banderwagon.zig").Fr;
const monomial_basis = @import("monomial_basis.zig");
const lagrange_basis = @import("lagrange_basis.zig");

pub fn PrecomputedWeights(
    comptime DomainSize: comptime_int,
    comptime domain: [DomainSize]Fr,
) type {
    return struct {
        const Self = @This();
        // TODO: We don't _need_ to store the vanishing polynomial
        //       we only need to store its derivative and whenever we need to evaluate
        //       the vanishing polynomial, it can be done via the domain

        // Vanishing polynomial
        A: monomial_basis.MonomialBasis(DomainSize),
        // Derivative of the vanishing polynomial
        Aprime: monomial_basis.MonomialBasis(DomainSize),
        // Aprime evaluated on the domain
        Aprime_DOMAIN: [crs.DomainSize]Fr,
        // Aprime evaluated on the domain and then inverted
        Aprime_DOMAIN_inv: [crs.DomainSize]Fr,
        //Inverse of the domain
        domain_inverses: [2 * crs.DomainSize]Fr,

        pub fn init() !Self {
            const _A = monomial_basis.MonomialBasis(DomainSize).vanishingPoly(domain);
            const _Aprime = _A.formalDerivative();
            var _Aprime_domain: [DomainSize]Fr = undefined;

            for (0..DomainSize) |i| {
                _Aprime_domain[i] = _Aprime.evaluate(Fr.fromInteger(i));
            }
            var _Aprime_domain_inv: [DomainSize]Fr = undefined;
            try Fr.batchInv(&_Aprime_domain_inv, &_Aprime_domain);

            // This is not fully correct as the first element will be the inverse of 0
            // We keep it this way for now because it is what the research code did
            // If domain size is 4 for example, the output would be:
            // [1/0, 1/1, 1/2, 1/3, -1/3, -1/2,-1/1]
            var inverses: [2 * DomainSize]Fr = undefined;
            inverses[0] = Fr.zero();
            for (1..DomainSize) |d| {
                inverses[d] = Fr.inv(Fr.fromInteger(d)) orelse Fr.zero();
                inverses[inverses.len - d] = Fr.inv(Fr.fromInteger(Fr.MODULO - d)) orelse Fr.zero();
            }
            return .{
                .A = _A,
                .Aprime = _Aprime,
                .Aprime_DOMAIN = _Aprime_domain,
                .Aprime_DOMAIN_inv = _Aprime_domain_inv,
                .domain_inverses = inverses,
            };
        }

        // barycentricFormularConstants returns a slice with the constants to be used when evaluating a polynomial at z.
        // b_i = A(z) / A'(DOMAIN[i]) * 1 / (z - DOMAIN[i])
        // The caller is responsible for freeing the returned slice.
        pub fn barycentricFormulaConstants(self: Self, z: Fr) ![DomainSize]Fr {
            std.debug.assert(z.toInteger() >= DomainSize);

            const Az = self.A.evaluate(z);

            var zSubX: [DomainSize]Fr = undefined;
            for (crs.Domain, 0..) |x, i| {
                zSubX[i] = Fr.sub(z, x);
            }
            var zSubXInvs: [DomainSize]Fr = undefined;
            try Fr.batchInv(&zSubXInvs, &zSubX);

            var r: [DomainSize]Fr = undefined;
            for (0..DomainSize) |i| {
                r[i] = Fr.mul(Az, Fr.mul(self.Aprime_DOMAIN_inv[i], zSubXInvs[i]));
            }

            return r;
        }
    };
}
