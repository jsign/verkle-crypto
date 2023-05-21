const std = @import("std");
const fp = @import("fp.zig");

pub fn main() !void {
    const f1 = getFpFrom(u32, 2);
    const f2 = getFpFrom(u32, 3);
    var out: fp.NonMontgomeryDomainFieldElement = undefined;

    const now = std.time.nanoTimestamp();
    fp.mul(&out, f1, f2);
    const now2 = std.time.nanoTimestamp();
    std.debug.print("Time taken: {}ns\n", .{now2 - now});

    const z = fpToU64(out);
    std.debug.print("Result is {}\n", .{z});
}

fn getFpFrom(comptime T: type, num: T) fp.MontgomeryDomainFieldElement {
    var lbe: [32]u8 = [_]u8{0} ** 32;
    std.mem.writeInt(u32, lbe[0..@divExact(@typeInfo(T).Int.bits, 8)], num, std.builtin.Endian.Little);

    var nonMont: fp.NonMontgomeryDomainFieldElement = undefined;
    fp.fromBytes(&nonMont, lbe);
    var mont: fp.MontgomeryDomainFieldElement = undefined;
    fp.toMontgomery(&mont, nonMont);

    return mont;
}

fn fpToU64(f: fp.MontgomeryDomainFieldElement) u64 {
    var mont: fp.MontgomeryDomainFieldElement = undefined;
    fp.fromMontgomery(&mont, f);
    var bytes: [32]u8 = [_]u8{0} ** 32;
    fp.toBytes(bytes[0..], mont);

    return std.mem.readInt(u64, bytes[0..8], std.builtin.Endian.Little);
}
