const FieldApi = @import("fieldapi.zig");

// Bandersnatch base and scalar finite fields.
pub const Fp = FieldApi.Fp;
pub const Fr = FieldApi.Fr;

// Curve parameters.
const A = Fp.fromInteger(Fp.MODULO - 5);
const D = Fp.fromInteger(138827208126141220649022263972958607803).Div(Fp.fromInteger(171449701953573178309673572579671231137));

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
