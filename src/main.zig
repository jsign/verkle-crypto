const std = @import("std");
const testing = std.testing;

// Lines 7-8 prints in two different ways
pub fn main() !void {
    var a: i512 = 10;
    var b: i512 = 2;
    std.debug.print("{}", .{@divTrunc(a, b)});
}

test "ecc/bandersnatch" {
    _ = @import("ecc/bandersnatch/bandersnatch.zig");
}
