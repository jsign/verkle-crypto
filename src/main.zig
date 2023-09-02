const std = @import("std");
const bander = @import("ecc/bandersnatch/bandersnatch.zig");

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
    _ = @import("polynomial/precomputed_Weights.zig");
    std.testing.refAllDeclsRecursive(@This());
}
