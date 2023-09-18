const std = @import("std");
const fieldapi = @import("ecc/bandersnatch/fieldapi.zig");
const Fp = fieldapi.Fp;

pub fn main() void {
    const N = 10_000;
    const fps: [N]Fp = get_fps(N);
    var start: i64 = undefined;

    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.legendre(fps[i]);
    }
    std.debug.print("Legendre symbol: {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.inv(fps[i]).?;
    }
    std.debug.print("Inverse: {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.inv(fps[i]).?;
    }
    std.debug.print("Sqrt: {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});
}

fn get_fps(comptime N: usize) [N]Fp {
    var fps: [N]Fp = undefined;
    var i: usize = 0;
    while (i < N) : (i += 1) {
        const fe = Fp.fromInteger(i);
        if (Fp.sqrt(fe) == null) {
            continue;
        }
        fps[i] = Fp.fromInteger(i);
    }
    return fps;
}
