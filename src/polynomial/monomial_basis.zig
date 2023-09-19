const std = @import("std");
const ArrayList = std.ArrayList;
const Fr = @import("../banderwagon/banderwagon.zig").Fr;

pub fn MonomialBasis(comptime PolyDegree: comptime_int) type {
    return struct {
        // [1,2,3] = 3x^2 + 2x + 1
        coeffs: [PolyDegree + 1]Fr,

        const Self = @This();

        pub fn fromCoeffsU256(u256_coeffs: [PolyDegree + 1]u256) Self {
            var coeffs: [PolyDegree + 1]Fr = undefined;
            for (u256_coeffs, 0..) |c, i| {
                coeffs[i] = Fr.fromInteger(c);
            }
            return .{ .coeffs = coeffs };
        }

        pub fn evaluate(self: Self, x: Fr) Fr {
            var y = Fr.zero();
            var power_of_x = Fr.one();

            for (self.coeffs) |p_coeff| {
                y = Fr.add(y, Fr.mul(power_of_x, p_coeff));
                power_of_x = Fr.mul(power_of_x, x);
            }
            return y;
        }

        pub fn formalDerivative(f: Self) Self {
            var coeffs: [PolyDegree + 1]Fr = undefined;
            coeffs[coeffs.len - 1] = Fr.zero();
            for (f.coeffs[1..], 0..) |c, i| {
                coeffs[i] = Fr.mul(c, Fr.fromInteger(i + 1));
            }
            return .{ .coeffs = coeffs };
        }

        pub fn vanishingPoly(xs: [PolyDegree]Fr) Self {
            var coeffs: [PolyDegree + 1]Fr = undefined;
            coeffs[0] = Fr.one();
            for (xs, 0..) |x, i| {
                for (0..i + 1) |k| {
                    coeffs[i + 1 - k] = coeffs[i + 1 - k - 1];
                }
                coeffs[0] = Fr.zero();
                for (0..i + 1) |j| {
                    coeffs[j] = Fr.sub(coeffs[j], Fr.mul(coeffs[j + 1], x));
                }
            }
            return .{ .coeffs = coeffs };
        }

        pub fn eq(self: Self, other: Self) bool {
            for (self.coeffs, other.coeffs) |a, b| {
                if (!a.equal(b)) return false;
            }
            return true;
        }

        pub fn print(self: Self) void {
            std.debug.print("MonomialBasis(", .{});
            for (self.coeffs) |c| {
                std.debug.print("{} ", .{c.toInteger()});
            }
            std.debug.print(")\n", .{});
        }
    };
}

test "Vanishing Polynomial on domain" {
    const xs = [_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    var z = MonomialBasis(xs.len).vanishingPoly(xs);

    for (xs) |x| {
        const eval = z.evaluate(x);
        try std.testing.expect(eval.isZero());
    }
    const eval = z.evaluate(Fr.fromInteger(6));
    try std.testing.expect(!eval.isZero());
}

test "Derivative" {
    // a = 6x^4 + 5x^3 + 10x^2 + 20x + 9
    const coeffs = [_]u256{ 9, 20, 10, 5, 6 };
    const a = MonomialBasis(coeffs.len - 1).fromCoeffsU256(coeffs);

    // the derivative of a is 24x^3 + 15x^2 + 20x + 20
    const deriv_coeffs = [_]u256{ 20, 20, 15, 24, 0 };
    const expected_a_prime = MonomialBasis(deriv_coeffs.len - 1).fromCoeffsU256(deriv_coeffs);

    try std.testing.expect(expected_a_prime.eq(a.formalDerivative()));
}
