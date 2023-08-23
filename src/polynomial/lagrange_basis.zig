const std = @import("std");
const ArrayList = std.ArrayList;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;
const MonomialBasis = @import("monomial_basis.zig").MonomialBasis;

pub const OperationError = error{LengthMismatch};

// TODO(jsign): see if we can avoid using an allocator.
pub const LagrangeBasis = struct {
    evaluations: ArrayList(Fr),
    domain: ArrayList(Fr),

    pub fn init(gpa: std.mem.Allocator, evaluations: []const Fr, domain: []const Fr) !LagrangeBasis {
        var evals = try ArrayList(Fr).initCapacity(gpa, evaluations.len);
        evals.appendSliceAssumeCapacity(evaluations);
        var dom = try ArrayList(Fr).initCapacity(gpa, domain.len);
        dom.appendSliceAssumeCapacity(domain);

        return .{
            .evaluations = evals,
            .domain = dom,
        };
    }

    pub fn deinit(self: LagrangeBasis) void {
        self.evaluations.deinit();
        self.domain.deinit();
    }

    pub fn values(self: *LagrangeBasis) ArrayList(Fr) {
        return try self.evaluations;
    }

    pub fn add(self: *LagrangeBasis, lhs: *const LagrangeBasis, rhs: *const LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        try self.evaluations.resize(lhs.evaluations.items.len);
        for (0..lhs.evaluations.items.len) |i| {
            self.evaluations.items[i] = lhs.evaluations.items[i].add(rhs.evaluations.items[i]);
        }
        try self.domain.resize(lhs.evaluations.items.len);
        std.mem.copy(Fr, self.domain.items, lhs.domain.items);
    }

    pub fn sub(self: *LagrangeBasis, lhs: *const LagrangeBasis, rhs: *const LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }

        try self.evaluations.resize(lhs.evaluations.items.len);
        for (0..lhs.evaluations.items.len) |i| {
            self.evaluations.items[i] = lhs.evaluations.items[i].sub(rhs.evaluations.items[i]);
        }

        try self.domain.resize(lhs.evaluations.items.len);
        std.mem.copy(Fr, self.domain.items, lhs.domain.items);
    }

    pub fn mul(self: *LagrangeBasis, lhs: *LagrangeBasis, rhs: *LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        self.evaluations.resize(lhs.evaluations.len);
        for (0..lhs.evaluations.items.len) |i| {
            self.evaluations = lhs.evaluations.items[i].mul(rhs.evaluations.items[i]);
        }
        try self.domain.resize(lhs.evaluations.items.len);
        std.mem.copy(Fr, self.domain.items, lhs.domain.items);
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

    pub fn interpolate(self: *LagrangeBasis, gpa: std.mem.Allocator) !MonomialBasis {
        const xs = self.domain;
        const ys = self.evaluations;

        // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
        const root = try MonomialBasis.vanishingPoly(xs);
        std.debug.assert(root.coeffs.len == ys.len + 1);

        // Generate per-value numerator polynomials, eg. for x=x2,
        // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
        // polynomial back by each x coordinate
        var nums = try ArrayList(MonomialBasis).initCapacity(gpa, xs.len);
        for (xs.items) |x| {
            const num = MonomialBasis.div(gpa, root, MonomialBasis.fromCoeffsFr(&[_]Fr{ -x, Fr.one() }));
            nums.appendAssumeCapacity(num);
        }

        //  Generate denominators by evaluating numerator polys at each x
        var denoms = try ArrayList(MonomialBasis).initCapacity(gpa, xs.len);
        for (nums.items, 0..) |num, i| {
            denoms.appendAssumeCapacity(num.evaluate(xs[i]));
        }
        const invdenoms = Fr.multiInv(denoms);

        // Generate output polynomial, which is the sum of the per-value numerator
        // polynomials rescaled to have the right y values
        var b = try ArrayList(Fr).initCapacity(gpa, ys.len);
        for (0..ys.len) |_| {
            b.appendAssumeCapacity(Fr.zero());
        }
        for (0..xs.items.len) |i| {
            const yslice = ys[i].mul(invdenoms[i]);
            for (0..ys.len) |j| {
                //if nums[i][j] and ys[i]:
                b[j] = b.add(nums[i][j].mul(yslice));
            }
        }
        // Remove zero terms from the highest degree until
        // we get to a non-zero term
        while (b[b.len - 1].IsZero()) {
            b.pop();
        }

        return MonomialBasis.fromCoeffsFr(b);
    }

    pub fn eq(lhs: *const LagrangeBasis, rhs: *const LagrangeBasis) bool {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            std.debug.print("{} != {}\n", .{ lhs.evaluations.items.len, rhs.evaluations.items.len });
            return false;
        }
        for (lhs.evaluations.items, 0..) |lhs_i, i| {
            if (!lhs_i.eq(rhs.evaluations.items[i])) {
                return false;
            }
        }
        // TODO(jsign): upstream suggestion to spec?
        if (lhs.domain.items.len != rhs.domain.items.len) {
            return false;
        }
        for (lhs.domain.items, 0..) |lhs_i, i| {
            if (!lhs_i.eq(rhs.domain.items[i])) {
                return false;
            }
        }
        return true;
    }
};

var allocator_test = std.testing.allocator;

test "add & sub" {
    const domain = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };

    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    const a = try LagrangeBasis.init(allocator_test, &x_squared, &domain);
    defer a.deinit();

    const x_plus_2 = [_]Fr{ Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5), Fr.fromInteger(6), Fr.fromInteger(7) };
    const b = try LagrangeBasis.init(allocator_test, &x_plus_2, &domain);
    defer b.deinit();

    var result = try LagrangeBasis.init(allocator_test, &[_]Fr{}, &[_]Fr{});
    defer result.deinit();
    try result.add(&a, &b);

    const expected_evaluations = [_]Fr{ Fr.fromInteger(2), Fr.fromInteger(4), Fr.fromInteger(8), Fr.fromInteger(14), Fr.fromInteger(22), Fr.fromInteger(32) };
    var expected_result = try LagrangeBasis.init(allocator_test, &expected_evaluations, &domain);
    defer expected_result.deinit();

    try std.testing.expect(expected_result.eq(&result));

    try expected_result.sub(&expected_result, &b);
    try std.testing.expect(expected_result.eq(&a));
}

//     def test_mul(self):
//         domain = [Fr(0), Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)]

//         # Evaluations
//         # x^2
//         x_squared = [Fr(0), Fr(1), Fr(4), Fr(9), Fr(16), Fr(25)]
//         # x^4
//         x_pow_4 = [Fr(0), Fr(1), Fr(16), Fr(81), Fr(256), Fr(625)]

//         a = Polynomial(x_squared, domain)

//         result = a * a

//         expected_result = Polynomial(x_pow_4, domain)

//         self.assertEqual(expected_result, result)

//     def test_scale(self):
//         domain = [Fr(0), Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)]

//         # Evaluations
//         # x^2
//         x_squared = [Fr(0), Fr(1), Fr(4), Fr(9), Fr(16), Fr(25)]
//         constant = Fr(10)

//         a = Polynomial(x_squared, domain)

//         result = a * constant

//         expected_evaluations = [
//             Fr(0), Fr(10), Fr(40), Fr(90), Fr(160), Fr(250)]

//         expected_result = Polynomial(expected_evaluations, domain)

//         self.assertEqual(expected_result, result)

//     def test_interpolation(self):
//         domain = [Fr(0), Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)]

//         # Evaluations
//         # x^2
//         x_squared = [Fr(0), Fr(1), Fr(4), Fr(9), Fr(16), Fr(25)]

//         x_squared_lagrange = Polynomial(x_squared, domain)
//         x_squared_coeff = x_squared_lagrange.interpolate()
//         expected_x_squared_coeff = MonomialBasis(
//             [Fr.zero(), Fr.zero(), Fr.one()])

//         self.assertEqual(expected_x_squared_coeff, x_squared_coeff)
