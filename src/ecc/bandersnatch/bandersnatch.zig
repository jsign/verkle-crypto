const std = @import("std");
const FieldApi = @import("fieldapi.zig");

// Bandersnatch base and scalar finite fields.
pub const Fp = FieldApi.Fp;
pub const Fr = FieldApi.Fr;

// Curve parameters.
pub const A = Fp.fromInteger(Fp.MODULO - 5);
pub const D = Fp.fromInteger(138827208126141220649022263972958607803).div(Fp.fromInteger(171449701953573178309673572579671231137)).?;

// Points.
pub const AffinePoint = @import("affinepoint.zig").AffinePoint;
pub const ExtendedPoint = @import("extendedpoint.zig").ExtendedPoint;

// Errors
pub const CurveError = error{
    NotInCurve,
};

test "fields" {
    _ = @import("fieldapi.zig");
}

test "addition" {
    const gen = ExtendedPoint.generator();
    const result_add = gen.add(gen);

    var result_double = ExtendedPoint.identity();
    result_double = gen.double();

    try std.testing.expect(result_add.equal(result_double));
}

test "equality" {
    const gen = ExtendedPoint.generator();
    const neg_gen = gen.neg();

    try std.testing.expect(gen.equal(gen));
    try std.testing.expect(!gen.equal(neg_gen));
}

test "neg" {
    const gen = ExtendedPoint.generator();
    const expected = ExtendedPoint.identity();

    const neg_gen = gen.neg();
    const result = neg_gen.add(gen);

    try std.testing.expect(expected.equal(result));
}

test "serialize gen" {
    const gen = ExtendedPoint.generator();
    const serialised_point = gen.toBytes();

    // test vector taken from the rust code (see spec reference)
    const expected = "18ae52a26618e7e1658499ad22c0792bf342be7b77113774c5340b2ccc32c129";
    const actual = std.fmt.bytesToHex(&serialised_point, std.fmt.Case.lower);
    try std.testing.expectEqualSlices(u8, expected, &actual);
}

test "scalar mul smoke" {
    const gen = ExtendedPoint.generator();

    const scalar = Fr.fromInteger(2);
    const result = gen.scalarMul(scalar);

    const twoG = ExtendedPoint.generator().double();

    try std.testing.expect(twoG.equal(result));
}

test "scalar mul minus one" {
    const gen = ExtendedPoint.generator();

    const integer = Fr.MODULO - 1;

    const scalar = Fr.fromInteger(integer);
    const result = gen.scalarMul(scalar);

    const expected = "e951ad5d98e7181e99d76452e0e343281295e38d90c602bf824892fd86742c4a";
    const actual = std.fmt.bytesToHex(result.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualSlices(u8, expected, &actual);
}
