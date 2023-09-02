const std = @import("std");
const ArrayList = std.ArrayList;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;
const MonomialBasis = @import("monomial_basis.zig").MonomialBasis;

pub const OperationError = error{LengthMismatch};

// TODO(jsign): see if we can avoid using an allocator.
pub const LagrangeBasis = struct {
    evaluations: ArrayList(Fr),
    domain: ArrayList(Fr), // TODO(jsign): slice

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

    pub fn empty(gpa: std.mem.Allocator) LagrangeBasis {
        var evals = ArrayList(Fr).init(gpa);
        var dom = ArrayList(Fr).init(gpa);

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

    pub fn mul(self: *LagrangeBasis, lhs: *const LagrangeBasis, rhs: *const LagrangeBasis) !void {
        if (lhs.evaluations.items.len != rhs.evaluations.items.len) {
            return OperationError.LengthMismatch;
        }
        try self.evaluations.resize(lhs.evaluations.items.len);
        for (0..lhs.evaluations.items.len) |i| {
            self.evaluations.items[i] = lhs.evaluations.items[i].mul(rhs.evaluations.items[i]);
        }
        try self.domain.resize(lhs.evaluations.items.len);
        std.mem.copy(Fr, self.domain.items, lhs.domain.items);
    }

    pub fn scale(self: *LagrangeBasis, poly: *const LagrangeBasis, constant: Fr) !void {
        try self.evaluations.resize(poly.evaluations.items.len);
        for (0..poly.evaluations.items.len) |i| {
            self.evaluations.items[i] = poly.evaluations.items[i].mul(constant);
        }
        try self.domain.resize(poly.evaluations.items.len);
        std.mem.copy(Fr, self.domain.items, poly.domain.items);
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

    pub fn interpolate(self: *const LagrangeBasis, gpa: std.mem.Allocator) !MonomialBasis {
        const xs = self.domain;
        const ys = self.evaluations;

        // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
        var root = try MonomialBasis.vanishingPoly(gpa, xs.items);
        defer root.deinit();

        std.debug.assert(root.coeffs.items.len == ys.items.len + 1);

        // Generate per-value numerator polynomials, eg. for x=x2,
        // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
        // polynomial back by each x coordinate
        var nums = try ArrayList(MonomialBasis).initCapacity(gpa, xs.items.len);
        defer {
            for (nums.items) |num| {
                num.deinit();
            }
            nums.deinit();
        }

        for (xs.items) |x| {
            var temp = try MonomialBasis.fromCoeffsFr(gpa, &[_]Fr{ x.neg(), Fr.one() });
            defer temp.deinit();
            const num = try MonomialBasis.div(gpa, root, temp);
            nums.appendAssumeCapacity(num);
        }

        //  Generate denominators by evaluating numerator polys at each x
        var denoms = try ArrayList(Fr).initCapacity(gpa, xs.items.len);
        defer denoms.deinit();
        for (nums.items, 0..) |num, i| {
            denoms.appendAssumeCapacity(num.evaluate(xs.items[i]));
        }
        const invdenoms = try Fr.multiInv(gpa, denoms.items);
        defer invdenoms.deinit();

        // Generate output polynomial, which is the sum of the per-value numerator
        // polynomials rescaled to have the right y values
        var b = try ArrayList(Fr).initCapacity(gpa, ys.items.len);
        defer b.deinit();

        for (0..ys.items.len) |_| {
            b.appendAssumeCapacity(Fr.zero());
        }
        for (0..xs.items.len) |i| {
            const yslice = ys.items[i].mul(invdenoms.items[i]);
            for (0..ys.items.len) |j| {
                //if nums[i][j] and ys[i]:
                b.items[j] = b.items[j].add(nums.items[i].coeffs.items[j].mul(yslice));
            }
        }
        // Remove zero terms from the highest degree until
        // we get to a non-zero term
        while (b.items[b.items.len - 1].isZero()) {
            _ = b.pop();
        }

        return try MonomialBasis.fromCoeffsFr(gpa, b.items);
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

test "mul" {
    // TODO(jsign): create some helper to do this.
    const domain = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };

    //# Evaluations
    //# x^2
    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    //# x^4
    const x_pow_4 = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(16), Fr.fromInteger(81), Fr.fromInteger(256), Fr.fromInteger(625) };

    const a = try LagrangeBasis.init(allocator_test, &x_squared, &domain);
    defer a.deinit();

    // TODO(jsign): create an .empty() API.
    var result = try LagrangeBasis.init(allocator_test, &[_]Fr{}, &[_]Fr{});
    defer result.deinit();

    try result.mul(&a, &a);

    const expected_result = try LagrangeBasis.init(allocator_test, &x_pow_4, &domain);
    defer expected_result.deinit();

    try std.testing.expect(expected_result.eq(&result));
}

test "scale" {
    const domain = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };

    // Evaluations
    // x^2
    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    const constant = Fr.fromInteger(10);

    const a = try LagrangeBasis.init(allocator_test, &x_squared, &domain);
    defer a.deinit();

    var result = LagrangeBasis.empty(allocator_test);
    defer result.deinit();
    try result.scale(&a, constant);

    const expected_evaluations = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(10), Fr.fromInteger(40), Fr.fromInteger(90), Fr.fromInteger(160), Fr.fromInteger(250) };

    const expected_result = try LagrangeBasis.init(allocator_test, &expected_evaluations, &domain);
    defer expected_result.deinit();

    try std.testing.expect(expected_result.eq(&result));
}

test "interpolation" {
    const domain = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };

    // Evaluations
    // x^2
    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };

    const x_squared_lagrange = try LagrangeBasis.init(allocator_test, &x_squared, &domain);
    defer x_squared_lagrange.deinit();

    const x_squared_coeff = try x_squared_lagrange.interpolate(allocator_test);
    defer x_squared_coeff.deinit();

    const expected_x_squared_coeff = try MonomialBasis.fromCoeffsFr(allocator_test, &[_]Fr{ Fr.zero(), Fr.zero(), Fr.one() });
    defer expected_x_squared_coeff.deinit();

    try std.testing.expect(expected_x_squared_coeff.eq(x_squared_coeff));
}
