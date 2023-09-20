const std = @import("std");
const ArrayList = std.ArrayList;
const Fr = @import("../banderwagon/banderwagon.zig").Fr;

pub fn LagrangeBasis(comptime domain_size: comptime_int, comptime eval_domain: [domain_size]Fr) type {
    return struct {
        const Self = @This();

        const domain = eval_domain;
        evaluations: [domain_size]Fr,

        pub fn init(evaluations: [domain_size]Fr) Self {
            return Self{ .evaluations = evaluations };
        }

        pub fn zero() Self {
            return Self{
                .evaluations = comptime blk: {
                    var zero_evals: [domain_size]Fr = undefined;
                    for (0..domain_size) |i| {
                        zero_evals[i] = Fr.zero();
                    }
                    break :blk zero_evals;
                },
            };
        }

        pub fn add(self: *Self, lhs: Self, rhs: Self) void {
            for (0..lhs.evaluations.len) |i| {
                self.evaluations[i] = lhs.evaluations[i].add(rhs.evaluations[i]);
            }
        }

        pub fn sub(self: *Self, lhs: Self, rhs: Self) void {
            for (0..lhs.evaluations.len) |i| {
                self.evaluations[i] = lhs.evaluations[i].sub(rhs.evaluations[i]);
            }
        }

        pub fn mul(self: *Self, lhs: Self, rhs: Self) void {
            for (0..lhs.evaluations.len) |i| {
                self.evaluations[i] = lhs.evaluations[i].mul(rhs.evaluations[i]);
            }
        }

        pub fn scale(self: *Self, poly: Self, constant: Fr) void {
            for (0..poly.evaluations.len) |i| {
                self.evaluations[i] = poly.evaluations[i].mul(constant);
            }
        }

        pub fn eq(lhs: Self, rhs: Self) bool {
            for (lhs.evaluations, 0..) |lhs_i, i| {
                if (!lhs_i.equal(rhs.evaluations[i])) {
                    return false;
                }
            }
            return true;
        }
    };
}

test "add & sub" {
    const domain = comptime [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    const TestLagrangeBasis = LagrangeBasis(domain.len, domain);

    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    const a = TestLagrangeBasis.init(x_squared);

    const x_plus_2 = [_]Fr{ Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5), Fr.fromInteger(6), Fr.fromInteger(7) };
    const b = TestLagrangeBasis.init(x_plus_2);

    var result = TestLagrangeBasis.zero();
    result.add(a, b);

    const expected_evaluations = [_]Fr{ Fr.fromInteger(2), Fr.fromInteger(4), Fr.fromInteger(8), Fr.fromInteger(14), Fr.fromInteger(22), Fr.fromInteger(32) };
    var expected_result = TestLagrangeBasis.init(expected_evaluations);

    try std.testing.expect(expected_result.eq(result));

    expected_result.sub(expected_result, b);
    try std.testing.expect(expected_result.eq(a));
}

test "mul" {
    const domain = comptime [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    const TestLagrangeBasis = LagrangeBasis(domain.len, domain);

    //# Evaluations
    //# x^2
    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    //# x^4
    const x_pow_4 = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(16), Fr.fromInteger(81), Fr.fromInteger(256), Fr.fromInteger(625) };

    const a = TestLagrangeBasis.init(x_squared);

    var result = TestLagrangeBasis.zero();
    result.mul(a, a);

    const expected_result = TestLagrangeBasis.init(x_pow_4);

    try std.testing.expect(expected_result.eq(result));
}

test "scale" {
    const domain = comptime [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    const TestLagrangeBasis = LagrangeBasis(domain.len, domain);

    // Evaluations
    // x^2
    const x_squared = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(4), Fr.fromInteger(9), Fr.fromInteger(16), Fr.fromInteger(25) };
    const constant = Fr.fromInteger(10);

    const a = TestLagrangeBasis.init(x_squared);

    var result = TestLagrangeBasis.zero();
    result.scale(a, constant);

    const expected_evaluations = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(10), Fr.fromInteger(40), Fr.fromInteger(90), Fr.fromInteger(160), Fr.fromInteger(250) };
    const expected_result = TestLagrangeBasis.init(expected_evaluations);

    try std.testing.expect(expected_result.eq(result));
}
