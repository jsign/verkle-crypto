const std = @import("std");
const Allocator = std.mem.Allocator;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const ElementNormalized = banderwagon.ElementNormalized;
const Fr = banderwagon.Fr;

pub fn Pippenger(comptime c: comptime_int) type {
    return struct {
        const window_mask = (1 << c) - 1;
        const num_windows = std.math.divCeil(u8, Fr.BitSize, c) catch unreachable;
        const num_buckets = (1 << c) - 1;

        pub fn msm(base_allocator: Allocator, basis: []const ElementNormalized, scalars_mont: []const Fr) !Element {
            std.debug.assert(basis.len >= scalars_mont.len);

            var arena = std.heap.ArenaAllocator.init(base_allocator);
            defer arena.deinit();
            var allocator = arena.allocator();

            var scalars = try allocator.alloc(u256, scalars_mont.len);
            for (0..scalars.len) |i| {
                scalars[i] = scalars_mont[i].toInteger();
            }

            var result: ?Element = null;
            var buckets: [num_buckets]?Element = std.mem.zeroes([num_buckets]?Element);
            var scalar_windows = try allocator.alloc(u16, scalars.len);
            for (0..num_windows) |w| {
                // Partition scalars.
                const w_idx = num_windows - w - 1;
                for (0..scalars.len) |i| {
                    scalar_windows[i] = @as(u16, @intCast((scalars[i] >> @as(u8, @intCast(w_idx * c))) & window_mask));
                }

                // Accumulate in buckets.
                for (0..buckets.len) |i| {
                    buckets[i] = null;
                }
                for (0..scalar_windows.len) |i| {
                    if (scalar_windows[i] == 0) {
                        continue;
                    }
                    if (buckets[scalar_windows[i] - 1] == null) {
                        buckets[scalar_windows[i] - 1] = Element.identity();
                    }
                    buckets[scalar_windows[i] - 1] = Element.mixedMsmAdd(buckets[scalar_windows[i] - 1].?, basis[i]);
                }

                // Aggregate buckets.
                var window_aggr: ?Element = null;
                var sum: ?Element = null;
                for (0..buckets.len) |i| {
                    if (window_aggr == null and buckets[buckets.len - 1 - i] == null) {
                        continue;
                    }
                    if (window_aggr == null) {
                        window_aggr = buckets[buckets.len - 1 - i];
                        sum = buckets[buckets.len - 1 - i];
                        continue;
                    }
                    if (buckets[buckets.len - 1 - i] != null) {
                        sum.?.add(sum.?, buckets[buckets.len - 1 - i].?);
                    }
                    window_aggr.?.add(window_aggr.?, sum.?);
                }

                // Aggregate into the final result.
                if (result != null) {
                    for (0..c) |_| {
                        result.?.double(result.?);
                    }
                }
                if (window_aggr != null) {
                    if (result == null) {
                        result = window_aggr.?;
                    } else {
                        result.?.add(result.?, window_aggr.?);
                    }
                }
            }

            return result orelse Element.identity();
        }
    };
}

test "correctness" {
    const crs = @import("../crs/crs.zig");
    const xcrs = try crs.CRS.init(std.testing.allocator);
    defer xcrs.deinit();

    var scalars: [crs.DomainSize]Fr = undefined;
    for (0..scalars.len) |i| {
        scalars[i] = Fr.fromInteger((i + 0x93434) *% 0x424242);
    }

    inline for (2..8) |c| {
        const pippenger = Pippenger(c);

        for (1..crs.DomainSize) |msm_length| {
            const msm_scalars = scalars[0..msm_length];

            var full_scalars: [crs.DomainSize]Fr = undefined;
            for (0..full_scalars.len) |i| {
                if (i < msm_length) {
                    full_scalars[i] = msm_scalars[i];
                    continue;
                }
                full_scalars[i] = Fr.zero();
            }
            const exp = xcrs.commitSlow(full_scalars);
            const got = try pippenger.msm(std.testing.allocator, xcrs.Gs[0..msm_length], msm_scalars);

            try std.testing.expect(Element.equal(exp, got));
        }
    }
}
