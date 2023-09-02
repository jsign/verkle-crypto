const std = @import("std");
const ArrayList = std.ArrayList;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;

// from dataclasses import dataclass
// from typing import List
// from ecc import Fr
// from copy import deepcopy

pub const MonomialBasis = struct {
    //     Represents polynomials in coefficient form
    //     The coefficient corresponding to the lowest
    //     degree monomial is stored in the lowest index
    //     ie [1,2,3] = 3x^2 + 2x + 1
    coeffs: ArrayList(Fr),

    pub fn empty(allocator: std.mem.Allocator) MonomialBasis {
        return .{ .coeffs = ArrayList(Fr).init(allocator) };
    }

    pub fn fromCoeffsFr(gpa: std.mem.Allocator, coeffs: []Fr) !MonomialBasis {
        var cfs = try ArrayList(Fr).initCapacity(gpa, coeffs.len);
        try cfs.appendSlice(coeffs);
        return .{ .coeffs = cfs };
    }

    pub fn fromCoeffsU256(gpa: std.mem.Allocator, coeffs: []const u256) !MonomialBasis {
        var cfs = try ArrayList(Fr).initCapacity(gpa, coeffs.len);
        for (coeffs) |c| {
            try cfs.append(Fr.fromInteger(c));
        }
        return .{ .coeffs = cfs };
    }

    pub fn deinit(self: *const MonomialBasis) void {
        self.coeffs.deinit();
    }

    pub fn mul(self: *MonomialBasis, allocator: std.mem.Allocator, a: MonomialBasis, b: MonomialBasis) *MonomialBasis {
        self.coeffs = ArrayList(Fr).initCapacity(allocator, a.coeffs.items.len + b.coeffs.len - 1);
        self.coeffs.appendNTimes(Fr.zero(), self.coeffs.items.capacity);

        for (a.coeffs.items, 0..) |aval, i| {
            for (b.coeffs.items, 0..) |bval, j| {
                self.coeffs.items[i + j] += aval * bval;
            }
        }
        return self;
    }
    pub fn div(allocator: std.mem.Allocator, a: MonomialBasis, b: MonomialBasis) !MonomialBasis {
        std.debug.assert(a.coeffs.items.len >= b.coeffs.items.len);

        const az = try a.coeffs.clone();
        defer az.deinit();
        var o = ArrayList(Fr).init(allocator);
        var apos = a.coeffs.items.len - 1;
        var bpos = b.coeffs.items.len - 1;
        var diff = apos - bpos;

        diffloop: while (diff >= 0) {
            const quot = try Fr.div(az.items[apos], b.coeffs.items[bpos]);
            try o.insert(0, quot);

            var i = bpos;
            blk: while (i >= 0) {
                az.items[diff + i] = Fr.sub(az.items[diff + i], Fr.mul(b.coeffs.items[i], quot));
                if (i == 0) {
                    break :blk;
                }
                i -= 1;
            }
            apos -= 1;

            if (diff == 0) {
                break :diffloop;
            }
            diff -= 1;
        }

        return .{ .coeffs = o };
    }

    pub fn evaluate(self: *const MonomialBasis, x: Fr) Fr {
        var y = Fr.zero();
        var power_of_x = Fr.one();

        for (self.coeffs.items) |p_coeff| {
            y = Fr.add(y, Fr.mul(power_of_x, p_coeff));
            power_of_x = Fr.mul(power_of_x, x);
        }
        return y;
    }

    pub fn formalDerivative(allocator: std.mem.Allocator, f: MonomialBasis) !MonomialBasis {
        var coeffs = try ArrayList(Fr).initCapacity(allocator, f.coeffs.items.len - 1);
        for (f.coeffs.items[1..], 1..) |c, n| {
            try coeffs.append(c.mul(Fr.fromInteger(n)));
        }
        return .{ .coeffs = coeffs };
    }

    pub fn vanishingPoly(allocator: std.mem.Allocator, xs: ArrayList(Fr)) !MonomialBasis {
        var root = ArrayList(Fr).init(allocator);
        try root.append(Fr.one());

        for (xs.items) |x| {
            try root.insert(0, Fr.zero());
            for (0..root.items.len - 1) |j| {
                root.items[j] = Fr.sub(root.items[j], Fr.mul(root.items[j + 1], x));
            }
        }
        return .{ .coeffs = root };
    }

    pub fn eq(self: MonomialBasis, other: MonomialBasis) bool {
        if (self.coeffs.items.len != other.coeffs.items.len) {
            return false;
        }

        for (self.coeffs.items, other.coeffs.items) |a, b| {
            if (!a.eq(b)) return false;
        }
        return true;
    }

    pub fn print(self: MonomialBasis) void {
        std.debug.print("MonomialBasis(", .{});
        for (self.coeffs.items) |c| {
            std.debug.print("{} ", .{c.toInteger()});
        }
        std.debug.print(")\n", .{});
    }
};

var allocator_test = std.testing.allocator;

test "Vanishing Polynomial on domain" {
    var xs = ArrayList(Fr).init(allocator_test);
    defer xs.deinit();
    try xs.appendSlice(&[_]Fr{ Fr.fromInteger(0), Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) });

    var z = try MonomialBasis.vanishingPoly(allocator_test, xs);
    defer z.deinit();

    for (xs.items) |x| {
        const eval = z.evaluate(x);
        try std.testing.expect(eval.isZero());
    }
    const eval = z.evaluate(Fr.fromInteger(6));
    try std.testing.expect(!eval.isZero());
}

test "Polynomial Division" {
    // a = (x+1)(x+2) = x^2 + 3x + 2
    var a = try MonomialBasis.fromCoeffsU256(allocator_test, &[_]u256{ 2, 3, 1 });
    defer a.deinit();

    // b = (x+1)
    var b = try MonomialBasis.fromCoeffsU256(allocator_test, &[_]u256{ 1, 1 });
    defer b.deinit();

    var result = try MonomialBasis.div(allocator_test, a, b);
    defer result.deinit();

    // Expected result should be (x+2)
    var expected = try MonomialBasis.fromCoeffsU256(allocator_test, &[_]u256{ 2, 1 });
    defer expected.deinit();

    try std.testing.expect(expected.eq(result));
}

test "Derivative" {
    // a = 6x^4 + 5x^3 + 10x^2 + 20x + 9
    var a = try MonomialBasis.fromCoeffsU256(allocator_test, &[_]u256{ 9, 20, 10, 5, 6 });
    defer a.deinit();

    // the derivative of a is 24x^3 + 15x^2 + 20x + 20
    var expected_a_prime = try MonomialBasis.fromCoeffsU256(allocator_test, &[_]u256{ 20, 20, 15, 24 });
    defer expected_a_prime.deinit();

    var got_a_prime = try MonomialBasis.formalDerivative(allocator_test, a);
    defer got_a_prime.deinit();

    try std.testing.expect(expected_a_prime.eq(got_a_prime));
}
