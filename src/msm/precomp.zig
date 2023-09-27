const std = @import("std");
const Allocator = std.mem.Allocator;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const Fr = banderwagon.Fr;
const bandersnatch = @import("../bandersnatch/bandersnatch.zig");
const ExtendedPoint = bandersnatch.ExtendedPoint;
const ExtendedPointNormalized = bandersnatch.ExtendedPointNormalized;

pub fn PrecompMSM(
    comptime _t: comptime_int,
    comptime _b: comptime_int,
) type {
    return struct {
        const Self = @This();

        const b = _b;
        const t = _t;
        const window_size = 1 << b;
        const points_per_column = (Fr.BitSize + t - 1) / t;

        allocator: Allocator,
        table: []const ExtendedPointNormalized,

        num_windows: usize,
        basis_len: usize,

        pub fn init(allocator: Allocator, basis: []const Element) !Self {
            const num_windows = (points_per_column * basis.len + b - 1) / b;
            var table_basis = try allocator.alloc(ExtendedPoint, points_per_column * basis.len);
            defer allocator.free(table_basis);
            var idx: usize = 0;
            for (0..basis.len) |hi| {
                table_basis[idx] = basis[hi].point;
                idx += 1;
                for (1..points_per_column) |_| {
                    table_basis[idx] = table_basis[idx - 1];
                    for (0..t) |_| {
                        table_basis[idx] = ExtendedPoint.double(table_basis[idx]);
                    }
                    idx += 1;
                }
            }

            var nn_table = try allocator.alloc(ExtendedPoint, window_size * num_windows);
            defer allocator.free(nn_table);
            for (0..num_windows) |w| {
                const start = w * b;
                var end = (w + 1) * b;
                if (end > table_basis.len) {
                    end = table_basis.len;
                }
                const window_basis = table_basis[start..end];
                fillWindow(window_basis, nn_table[w * window_size .. (w + 1) * window_size]);
            }

            var table = try allocator.alloc(ExtendedPointNormalized, window_size * num_windows);
            try ExtendedPointNormalized.fromExtendedPoints(table, nn_table);

            return Self{
                .allocator = allocator,
                .table = table,
                .num_windows = num_windows,
                .basis_len = basis.len,
            };
        }

        pub fn deinit(self: Self) void {
            self.allocator.free(self.table);
        }

        pub fn msm(self: Self, mont_scalars: []const Fr) !Element {
            if (mont_scalars.len > self.basis_len) {
                return error.ScalarLengthBiggerThanBasis;
            }

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
                        if (k + t - t_i - 1 < Fr.BitSize) {
                            const bit = scalars[s_i] >> (@as(u8, @intCast(k + t - t_i - 1))) & 1;
                            curr_window_scalar |= @as(usize, @intCast(bit << (b - curr_window_b_idx - 1)));
                        }
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
                if (curr_window_scalar > 0) {
                    accum = bandersnatch.ExtendedPoint.mixedAdd(accum, self.table[curr_window_idx * window_size .. (curr_window_idx + 1) * window_size][curr_window_scalar]);
                }
            }

            return Element{ .point = accum };
        }

        fn fillWindow(basis: []const ExtendedPoint, table: []ExtendedPoint) void {
            if (basis.len == 0) {
                for (0..table.len) |i| {
                    table[i] = ExtendedPoint.identity();
                }
                return;
            }
            fillWindow(basis[1..], table[0 .. table.len / 2]);
            for (0..table.len / 2) |i| {
                table[table.len / 2 + i] = ExtendedPoint.add(table[i], basis[0]);
            }
        }
    };
}

test "correctness" {
    const crs = @import("../crs/crs.zig");
    const xcrs = try crs.CRS.init(std.testing.allocator);
    defer xcrs.deinit();

    const precomp = try PrecompMSM(2, 5).init(std.testing.allocator, &xcrs.Gs);
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
        const exp = xcrs.commitSlow(full_scalars);
        const got = try precomp.msm(msm_scalars);

        try std.testing.expect(Element.equal(exp, got));
    }
}
