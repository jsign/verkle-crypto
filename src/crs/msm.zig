const std = @import("std");
const Allocator = std.mem.Allocator;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const Fr = banderwagon.Fr;
const bandersnatch = @import("../bandersnatch/bandersnatch.zig");
const ElementNormalized = bandersnatch.ExtendedPointNormalized;

pub fn PrecompMSM(
    comptime _t: comptime_int,
    comptime _b: comptime_int,
    comptime basis_len: comptime_int,
) type {
    return struct {
        const Self = @This();

        const b = _b;
        const t = _t;
        const window_size = 1 << b;
        const points_per_column = (Fr.BitSize + t - 1) / t;
        const num_windows = (points_per_column * basis_len + b - 1) / b;

        allocator: Allocator,
        table: []const ElementNormalized,

        pub fn init(allocator: Allocator, basis: [basis_len]Element) !Self {
            var table_basis = try allocator.alloc(Element, num_windows * basis_len);
            defer allocator.free(table_basis);
            var idx: usize = 0;
            for (0..basis_len) |hi| {
                table_basis[idx] = basis[hi];
                idx += 1;
                for (1..points_per_column) |_| {
                    table_basis[idx] = table_basis[idx - 1];
                    for (0..t) |_| {
                        table_basis[idx].double(table_basis[idx]);
                    }
                    idx += 1;
                }
            }

            var nn_table = try allocator.alloc(Element, window_size * num_windows);
            defer allocator.free(nn_table);
            for (0..num_windows) |w| {
                const start = w * b;
                var end = (w + 1) * b;
                if (end > table_basis.len) {
                    end = table_basis.len;
                }
                const window_basis = table_basis[start..end];
                fill_window(window_basis, nn_table[w * window_size .. (w + 1) * window_size]);
            }

            var table = try allocator.alloc(ElementNormalized, window_size * num_windows);
            // TODO: batch.
            for (nn_table, 0..) |p, i| {
                table[i] = ElementNormalized.fromExtendedPoint(p.point);
            }

            return Self{
                .allocator = allocator,
                .table = table,
            };
        }

        pub fn deinit(self: Self) void {
            self.allocator.free(self.table);
        }

        pub fn msm(self: Self, mont_scalars: []const Fr) !Element {
            std.debug.assert(mont_scalars.len <= basis_len);

            var scalars = try self.allocator.alloc(u256, mont_scalars.len);
            defer self.allocator.free(scalars);
            for (0..mont_scalars.len) |i| {
                scalars[i] = mont_scalars[i].toInteger();
            }

            var accum = bandersnatch.ExtendedPoint.identity();
            for (0..t) |t_i| {
                if (t_i > 0) {
                    accum = bandersnatch.ExtendedPoint.double(accum);
                }

                var curr_window_idx: usize = 0;
                var curr_window_scalar: usize = 0;
                var curr_window_b_idx: u8 = 0;
                for (0..scalars.len) |s_i| {
                    var k: u16 = 0;
                    while (k < Fr.BitSize) : (k += t) {
                        const bit = scalars[s_i] >> (@as(u8, @intCast(k + t - t_i - 1))) & 1;
                        curr_window_scalar |= @as(usize, @intCast(bit << (b - curr_window_b_idx - 1)));
                        curr_window_b_idx += 1;

                        if (curr_window_b_idx == b) {
                            if (curr_window_scalar > 0) {
                                accum = bandersnatch.ExtendedPoint.mixedAdd(accum, self.table[curr_window_idx * window_size .. (curr_window_idx + 1) * window_size][curr_window_scalar]);
                            }
                            curr_window_idx += 1;

                            curr_window_scalar = 0;
                            curr_window_b_idx = 0;
                        }
                    }
                }
            }

            return Element{ .point = accum };
        }

        fn fill_window(basis: []const Element, table: []Element) void {
            if (basis.len == 0) {
                for (0..table.len) |i| {
                    table[i] = Element.identity();
                }
                return;
            }
            fill_window(basis[1..], table[0 .. table.len / 2]);
            for (0..table.len / 2) |i| {
                table[table.len / 2 + i].add(table[i], basis[0]);
            }
        }
    };
}

test "correctness" {
    const crs = @import("crs.zig");
    const CRS = crs.CRS.init();
    var test_allocator = std.testing.allocator;

    const precomp = try PrecompMSM(2, 5, crs.DomainSize).init(test_allocator, CRS.Gs);
    defer precomp.deinit();

    var scalars: [crs.DomainSize]Fr = undefined;
    for (0..scalars.len) |i| {
        scalars[i] = Fr.fromInteger((i + 0x93434) *% 0x424242);
    }

    for (1..crs.DomainSize + 1) |msm_length| {
        const msm_scalars = scalars[0..msm_length];

        var full_scalars: [crs.DomainSize]Fr = undefined;
        for (0..full_scalars.len) |i| {
            if (i < msm_length) {
                full_scalars[i] = msm_scalars[i];
                continue;
            }
            full_scalars[i] = Fr.zero();
        }
        const exp = CRS.commit(full_scalars);
        const got = try precomp.msm(msm_scalars);

        try std.testing.expect(Element.equal(exp, got));
    }
}
