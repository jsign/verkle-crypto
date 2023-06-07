const std = @import("std");
const bander = @import("ecc/bandersnatch/bandersnatch.zig");

pub fn main() !void {}

test "bandersnatch" {
    _ = @import("ecc/bandersnatch/bandersnatch.zig");
    std.testing.refAllDeclsRecursive(@This());
}
