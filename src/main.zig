const std = @import("std");

pub fn main() !void {}

test "bandersnatch" {
    _ = @import("ecc/bandersnatch/bandersnatch.zig");
    std.testing.refAllDeclsRecursive(@This());
}

test "banderwagon" {
    _ = @import("ecc/bandersnatch/banderwagon.zig");
    std.testing.refAllDeclsRecursive(@This());
}

test "polynomial" {
    _ = @import("polynomial/monomial_basis.zig");
    _ = @import("polynomial/lagrange_basis.zig");
    _ = @import("polynomial/precomputed_weights.zig");
    std.testing.refAllDeclsRecursive(@This());
}

test "crs" {
    _ = @import("crs/crs.zig");
    std.testing.refAllDeclsRecursive(@This());
}
