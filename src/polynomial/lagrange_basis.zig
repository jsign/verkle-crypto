const std = @import("std");
const ArrayList = std.ArrayList;
const Bandersnatch = @import("ecc.zig");
const Fr = Bandersnatch.Fr;

pub const OperationError = error{LengthMismatch};
// from dataclasses import dataclass
// from typing import List
// from ecc import Fr
// from copy import deepcopy
// from .monomial_basis import MonomialBasis

pub const LagrangeBasis = struct {
    evaluations: ArrayList(Fr),
    domain: ArrayList(Fr),

    pub fn empty() LagrangeBasis {
        return .{
            .evaluations = ArrayList(Fr).init(),
            .domain = ArrayList(Fr).init(),
        };
    }

    pub fn values(self: *LagrangeBasis) ArrayList(Fr) {
        return try self.evaluations.clone();
    }

    pub fn add(self: *LagrangeBasis, lhs: *LagrangeBasis, rhs: *LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        self.evaluations = try ArrayList(Fr).initCapacity(lhs.evaluations.len);
        for (lhs.evaluations.items, 0..) |lhs_i, i| {
            self.evaluations.appendAssumeCapacity(lhs_i + rhs.evaluations.items[i]);
        }
        self.domain = try lhs.domain.clone();
    }

    pub fn sub(self: *LagrangeBasis, lhs: *LagrangeBasis, rhs: *LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        self.evaluations = try ArrayList(Fr).initCapacity(lhs.evaluations.len);
        for (lhs.evaluations.items, 0..) |lhs_i, i| {
            self.evaluations.appendAssumeCapacity(lhs_i - rhs.evaluations.items[i]);
        }
        self.domain = try lhs.domain.clone();
    }

    pub fn mul(self: *LagrangeBasis, lhs: *LagrangeBasis, rhs: *LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        self.evaluations = try ArrayList(Fr).initCapacity(lhs.evaluations.len);
        for (lhs.evaluations.items, 0..) |lhs_i, i| {
            self.evaluations.appendAssumeCapacity(lhs_i * rhs.evaluations.items[i]);
        }
        self.domain = try lhs.domain.clone();
    }

    pub fn scale(self: *LagrangeBasis, poly: *LagrangeBasis, constant: Fr) !void {
        self.evaluations = try ArrayList(Fr).initCapacity(poly.evaluations.len);
        for (poly.evaluations.items) |eval| {
            self.evaluations.appendAssumeCapacity(eval * constant);
        }
        self.domain = try poly.domain.clone();
    }

    //     # TODO: we cannot add the type PrecomputedWeights because it
    //     # TODO: will trigger a circular import
    //     # TODO: we could take as a parameter: Aprime_DOMAIN_inv
    //     def evaluate_outside_domain(self, precomputed_weights, z: Fr):

    //         r = Fr.zero()
    //         A = MonomialBasis.vanishing_poly(self.domain)
    //         Az = A.evaluate(z)

    //         if Az.is_zero() == True:
    //             raise Exception(
    //                 "vanishing polynomial evaluated to zero. z is therefore a point on the domain")

    //         inverses = Fr.multi_inv([z - x for x in self.domain])

    //         for i, x in enumerate(inverses):
    //             r += self[i] * precomputed_weights.Aprime_DOMAIN_inv[i] * x

    //         r = r * Az

    //         return r

    pub fn interpolate(self: *LagrangeBasis) !MonomialBasis {
                const xs = self.domain;
                const ys = self.evaluations;

                // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
                root = MonomialBasis.vanishing_poly(xs)
                assert len(root) == len(ys) + 1

                // Generate per-value numerator polynomials, eg. for x=x2,
                // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
                // polynomial back by each x coordinate
                nums = [root / MonomialBasis([-x, Fr.one()]) for x in xs]

                //  Generate denominators by evaluating numerator polys at each x
                denoms = [nums[i].evaluate(xs[i]) for i in range(len(xs))]
                invdenoms = Fr.multi_inv(denoms)

                // Generate output polynomial, which is the sum of the per-value numerator
                // polynomials rescaled to have the right y values
                b = [Fr.zero()] * len(ys)
                for i in range(len(xs)):
                    yslice = ys[i] * invdenoms[i]
                    for j in range(len(ys)):
                        if nums[i][j] and ys[i]:
                            b[j] += nums[i][j] * yslice

                // Remove zero terms from the highest degree until
                // we get to a non-zero term
                while b[-1].is_zero():
                    b.pop()

                return MonomialBasis(b)
    }

    //     def equal(self, other):
    //         assert(isinstance(other, LagrangeBasis))
    //         for lhs_i, rhs_i in zip(self.evaluations, other.evaluations):
    //             if lhs_i != rhs_i:
    //                 return False
    //         return True
};
