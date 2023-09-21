const std = @import("std");
const Allocator = std.mem.Allocator;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Element = banderwagon.Element;
const Fr = banderwagon.Fr;

const PrecompMSM = struct {
    allocator: Allocator,
    b: usize,
    basis_len: usize,
    table: []const Element,

    pub fn init(allocator: Allocator, basis: []const Element, b: usize) !PrecompMSM {
        std.debug.assert(basis.len % b == 0);

        const window_size = std.math.shl(usize, 1, b);
        const num_windows = basis.len / b;
        const table_num_elements = window_size * num_windows;
        var table = try allocator.alloc(Element, table_num_elements);

        for (0..num_windows) |w| {
            const window_basis = basis[w * b .. (w + 1) * b];
            fill_window(window_basis, table[w * window_size .. (w + 1) * window_size]);
        }
        return PrecompMSM{
            .allocator = allocator,
            .b = b,
            .basis_len = basis.len,
            .table = table,
        };
    }

    pub fn deinit(self: PrecompMSM) void {
        self.allocator.free(self.table);
    }

    pub fn msm(self: PrecompMSM, mont_scalars: []const Fr) !Element {
        std.debug.assert(mont_scalars.len <= self.basis_len);

        var scalars = try self.allocator.alloc(u256, mont_scalars.len);
        defer self.allocator.free(scalars);
        for (0..mont_scalars.len) |i| {
            scalars[i] = mont_scalars[i].toInteger();
        }

        const window_size = std.math.shl(usize, 1, self.b);
        const num_windows = self.basis_len / self.b;
        var accum = Element.identity();
        for (0..253) |k| {
            accum.double(accum);
            for (0..num_windows) |w| {
                if (w * self.b < scalars.len) {
                    const window_scalars = scalars[w * self.b ..];
                    var table_idx: usize = 0;
                    for (0..self.b) |i| {
                        table_idx <<= 1;
                        if (i < window_scalars.len) {
                            table_idx |= @as(u1, @truncate((window_scalars[i] >> @as(u8, @intCast(252 - k)))));
                        }
                    }
                    const window_table = self.table[w * window_size .. (w + 1) * window_size];
                    accum.add(accum, window_table[table_idx]);
                }
            }
        }

        return accum;
    }

    fn fill_window(basis: []const Element, table: []Element) void {
        if (basis.len == 0) {
            table[0] = Element.identity();
            return;
        }
        fill_window(basis[1..], table[0 .. table.len / 2]);
        for (0..table.len / 2) |i| {
            table[table.len / 2 + i].add(table[i], basis[0]);
        }
    }
};

test "correctness" {
    const crs = @import("crs.zig");
    const CRS = crs.CRS.init();
    var test_allocator = std.testing.allocator;

    const precomp = try PrecompMSM.init(test_allocator, &CRS.Gs, 8);
    defer precomp.deinit();

    var scalars: [crs.DomainSize]Fr = undefined;
    for (0..scalars.len) |i| {
        scalars[i] = Fr.fromInteger(i + 0x424242);
    }

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
        std.debug.print("For {} ...", .{msm_length});
        const exp = CRS.commit(full_scalars);
        std.debug.print("ok", .{});

        const got = try precomp.msm(msm_scalars);
        std.debug.print(" ok\n", .{});

        try std.testing.expect(Element.equal(exp, got));
    }
}
