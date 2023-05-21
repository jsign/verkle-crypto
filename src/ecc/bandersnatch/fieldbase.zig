const std = @import("std");
const gen_fp = @import("gen_fp.zig");

const MODULO: u256 = 52435875175126190479447740508185965837690552500527637822603658699938581184513;
const Q_MIN_ONE_DIV_2 = (MODULO - 1) / 2;
const BYTE_LEN = 32;

const baseZero = val: {
    var bz: gen_fp.MontgomeryDomainFieldElement = undefined;
    gen_fp.fromBytes(&bz, [_]u8{0} ** BYTE_LEN);
    break :val Fp{ .fe = bz };
};

comptime {
    std.debug.assert(@bitSizeOf(u256) == BYTE_LEN * 8);
}

const Fp = struct {
    fe: gen_fp.MontgomeryDomainFieldElement,

    fn fromInteger(num: u256) Fp {
        var lbe: [BYTE_LEN]u8 = [_]u8{0} ** BYTE_LEN;
        std.mem.writeInt(u256, lbe[0..], num % MODULO, std.builtin.Endian.Little);

        var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
        gen_fp.fromBytes(&nonMont, lbe);
        var mont: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.toMontgomery(&mont, nonMont);

        return Fp{ .fe = mont };
    }

    pub fn zero() Fp {
        return baseZero;
    }
    pub fn one() Fp {
        return comptime {
            var baseOne: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.setOne(&baseOne);
            return Fp{ .fe = baseOne };
        };
    }

    pub fn fromBytes(bytes: [BYTE_LEN]u8) Fp {
        var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.fromBytes(&ret, bytes);
        return Fp{ .fe = ret };
    }

    // TODO
    //pub fn fromBytesReduce()
    //     def from_bytes_reduce(bytes):
    //         return Fp(None, Field.from_bytes_reduce(bytes, BASE_FIELD))

    pub fn toBytes(self: Fp) [BYTE_LEN]u8 {
        var ret: [BYTE_LEN]u8 = undefined;
        gen_fp.toBytes(&ret, self.fe);
        return ret;
    }

    pub fn lexographically_largest(self: Fp) bool {
        const qminusonediv2 = comptime fromInteger(Q_MIN_ONE_DIV_2);
        for (self.fe, 0..) |elem, i| {
            if (elem > qminusonediv2.fe[i]) return true;
        }
        return false;
    }

    // TODO
    //     def multi_inv(values):
    //         result = []
    //         inverses = Field.multi_inv(values)
    //         for inv in inverses:
    //         result.append(Fp(None, inv))
    //         return result

    pub fn add(self: Fp, other: Fp) Fp {
        var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.add(&ret, self.fe, other.fe);
        return Fp{ .fe = ret };
    }

    pub fn sub(self: Fp, other: Fp) Fp {
        var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.sub(&ret, self.fe, other.fe);
        return Fp{ .fe = ret };
    }

    pub fn mul(self: Fp, other: Fp) Fp {
        var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.mul(&ret, self.fe, other.fe);
        return Fp{ .fe = ret };
    }

    pub fn neg(self: Fp) Fp {
        var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
        gen_fp.sub(&ret, baseZero.fe, self.fe);
        return Fp{ .fe = ret };
    }

    pub fn isZero(self: Fp) bool {
        return self.eq(baseZero);
    }

    // TODO: this is naive, do something better.
    pub fn pow(self: Fp, exponent: u256) Fp {
        var res: Fp = self;
        var exp = exponent;
        while (exponent > 1) {
            res = res.mul(self);
            exp -= 1;
        }
        return res;
    }

    pub fn inv(self: Fp) ?Fp {
        var r_old: i512 = MODULO;
        var t_old: i512 = 0;

        var r_new: i512 = self.toInteger();
        if (r_new != 1) {
            std.debug.print("r_new is not 1: {}\n", .{r_new});
        }
        var t_new: i512 = 1;

        while (r_new != 0) {
            const q = @divTrunc(r_old, r_new);
            const s = t_old - q * t_new;
            const r = @mod(r_old, r_new);

            r_old = r_new;
            t_old = t_new;
            r_new = r;
            t_new = s;
        }

        // Not invertible
        if (r_old > 1) {
            return null;
        }

        if (t_old < 0) {
            return Fp.fromInteger(@intCast(u256, t_old) + MODULO);
        }

        return Fp.fromInteger(@intCast(u256, t_old));
    }

    pub fn div(self: Fp, den: Fp) ?Fp {
        const denInv = den.inv() orelse return null;
        return self.mul(denInv);
    }

    pub fn eq(self: Fp, other: Fp) bool {
        return std.mem.eql(u64, &self.fe, &other.fe);
    }

    fn toInteger(self: Fp) u256 {
        var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
        gen_fp.fromMontgomery(&nonMont, self.fe);

        var bytes: [BYTE_LEN]u8 = [_]u8{0} ** BYTE_LEN;
        gen_fp.toBytes(&bytes, nonMont);

        return std.mem.readInt(u256, &bytes, std.builtin.Endian.Little);
    }
};

test "one" {
    const oneFromInteger = Fp.fromInteger(1);
    const oneFromAPI = Fp.one();

    try std.testing.expect(oneFromInteger.eq(oneFromAPI));
}

test "zero" {
    const zeroFromInteger = Fp.fromInteger(0);
    const zeroFromAPI = Fp.zero();

    try std.testing.expect(zeroFromInteger.eq(zeroFromAPI));
}

test "lexographically largest" {
    try std.testing.expect(!Fp.fromInteger(0).lexographically_largest());
    try std.testing.expect(!Fp.fromInteger(Q_MIN_ONE_DIV_2).lexographically_largest());

    try std.testing.expect(Fp.fromInteger(Q_MIN_ONE_DIV_2 + 1).lexographically_largest());
    try std.testing.expect(Fp.fromInteger(MODULO - 1).lexographically_largest());
}

test "from and to bytes" {
    const cases = [_]Fp{ Fp.fromInteger(0), Fp.fromInteger(1), Fp.fromInteger(Q_MIN_ONE_DIV_2), Fp.fromInteger(MODULO - 1) };

    for (cases) |fe| {
        const bytes = fe.toBytes();
        const fe2 = Fp.fromBytes(bytes);
        try std.testing.expect(fe.eq(fe2));

        const bytes2 = fe2.toBytes();
        try std.testing.expectEqualSlices(u8, &bytes, &bytes2);
    }
}

test "to integer" {
    try std.testing.expect(Fp.fromInteger(0).toInteger() == 0);
    try std.testing.expect(Fp.fromInteger(1).toInteger() == 1);
    try std.testing.expect(Fp.fromInteger(100).toInteger() == 100);
}

test "add sub mul neg" {
    const got = Fp.fromInteger(10).mul(Fp.fromInteger(20)).add(Fp.fromInteger(30)).sub(Fp.fromInteger(40)).add(Fp.fromInteger(MODULO));
    const want = Fp.fromInteger(190);
    try std.testing.expect(got.eq(want));

    const gotneg = got.neg();
    const wantneg = Fp.fromInteger(MODULO - 190);
    try std.testing.expect(gotneg.eq(wantneg));
}

test "inv" {
    try std.testing.expect(Fp.fromInteger(0).inv() == null);

    const one = Fp.one();
    const cases = [_]Fp{Fp.fromInteger(1)}; //, Fp.fromInteger(Q_MIN_ONE_DIV_2), Fp.fromInteger(MODULO - 1) };
    for (cases) |fe| {
        try std.testing.expect(fe.mul(fe.inv().?).eq(one));
    }
}

test "div" {}
