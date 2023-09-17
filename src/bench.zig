const std = @import("std");
const fieldapi = @import("ecc/bandersnatch/fieldapi.zig");
const Fp = fieldapi.Fp;

pub fn main() void {
    const a = Fp.fromInteger(12345);

    const n = 1_000;
    var start = std.time.microTimestamp();
    for (0..n) |_| {
        _ = Fp.legendre(a);
    }
    std.debug.print("Legendre symbol: {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (n))});
}
