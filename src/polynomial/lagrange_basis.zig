const std = @import("std");
const ArrayList = std.ArrayList;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;

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
