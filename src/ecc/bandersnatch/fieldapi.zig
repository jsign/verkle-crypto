const std = @import("std");
const gen_fp = @import("gen_fp.zig");
const gen_fr = @import("gen_fr.zig");

pub const Fp = BandersnatchField(gen_fp.MontgomeryDomainFieldElement, 52435875175126190479447740508185965837690552500527637822603658699938581184513);
pub const Fr = BandersnatchField(gen_fr.MontgomeryDomainFieldElement, 13108968793781547619861935127046491459309155893440570251786403306729687672801);

fn BandersnatchField(comptime fieldType: type, comptime mod: u256) type {
    const BYTE_LEN = 32;
    comptime {
        std.debug.assert(@bitSizeOf(u256) == BYTE_LEN * 8);
    }

    return struct {
        pub const MODULO = mod;

        const Self = @This();
        const Q_MIN_ONE_DIV_2 = (MODULO - 1) / 2;
        const baseZero = val: {
            var bz: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.fromBytes(&bz, [_]u8{0} ** BYTE_LEN);
            break :val Self{ .fe = bz };
        };

        fe: fieldType,

        pub fn fromInteger(num: u256) Self {
            var lbe: [BYTE_LEN]u8 = [_]u8{0} ** BYTE_LEN;
            std.mem.writeInt(u256, lbe[0..], num % MODULO, std.builtin.Endian.Little);

            var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
            gen_fp.fromBytes(&nonMont, lbe);
            var mont: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.toMontgomery(&mont, nonMont);

            return Self{ .fe = mont };
        }

        pub fn zero() Self {
            return baseZero;
        }

        pub fn one() Self {
            return comptime {
                var baseOne: gen_fp.MontgomeryDomainFieldElement = undefined;
                gen_fp.setOne(&baseOne);
                return Self{ .fe = baseOne };
            };
        }

        pub fn fromBytes(bytes: [BYTE_LEN]u8) Self {
            var dret: gen_fp.NonMontgomeryDomainFieldElement = undefined;
            gen_fp.fromBytes(&dret, bytes);

            var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.toMontgomery(&ret, dret);
            return Self{ .fe = ret };
        }

        pub fn toBytes(self: Self) [BYTE_LEN]u8 {
            var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
            gen_fp.fromMontgomery(&nonMont, self.fe);

            var ret: [BYTE_LEN]u8 = undefined;
            gen_fp.toBytes(&ret, nonMont);
            return ret;
        }

        // TODO
        //pub fn fromBytesReduce()
        //     def from_bytes_reduce(bytes):
        //         return Fp(None, Field.from_bytes_reduce(bytes, BASE_FIELD))

        pub fn lexographically_largest(self: Self) bool {
            const selfNonMont = self.fromMontgomery();
            const qminusonediv2 = comptime fromInteger(Q_MIN_ONE_DIV_2).fromMontgomery();
            inline for (selfNonMont, 0..) |elem, i| {
                if (elem > qminusonediv2[i]) return true;
            }
            return false;
        }

        pub fn fromMontgomery(self: Self) gen_fp.NonMontgomeryDomainFieldElement {
            var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
            gen_fp.fromMontgomery(&nonMont, self.fe);
            return nonMont;
        }

        // TODO
        //     def multi_inv(values):
        //         result = []
        //         inverses = Field.multi_inv(values)
        //         for inv in inverses:
        //         result.append(Fp(None, inv))
        //         return result

        pub fn add(self: Self, other: Self) Self {
            var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.add(&ret, self.fe, other.fe);
            return Self{ .fe = ret };
        }

        pub fn sub(self: Self, other: Self) Self {
            var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.sub(&ret, self.fe, other.fe);
            return Self{ .fe = ret };
        }

        pub fn mul(self: Self, other: Self) Self {
            var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.mul(&ret, self.fe, other.fe);
            return Self{ .fe = ret };
        }

        pub fn neg(self: Self) Self {
            var ret: gen_fp.MontgomeryDomainFieldElement = undefined;
            gen_fp.sub(&ret, baseZero.fe, self.fe);
            return Self{ .fe = ret };
        }

        pub fn isZero(self: Self) bool {
            return self.eq(baseZero);
        }

        pub fn isOne(self: Self) bool {
            return self.eq(one());
        }

        // TODO(improv): this is naive, do something better.
        pub fn pow(self: Self, exponent: u256) Self {
            var res: Self = self;
            var exp = exponent;
            while (exponent > 1) {
                res = res.mul(self);
                exp -= 1;
            }
            return res;
        }

        pub fn inv(self: Self) ?Self {
            var r: u256 = MODULO;
            var t: i768 = 0;

            var newr: u256 = self.toInteger();
            var newt: i768 = 1;

            while (newr != 0) {
                const quotient = r / newr;
                const tempt = t - quotient * newt;
                const tempr = @mod(r, newr);

                r = newr;
                t = newt;
                newr = tempr;
                newt = tempt;
            }

            // Not invertible
            if (r > 1) {
                return null;
            }

            if (t < 0) {
                return Self.fromInteger(@intCast(u256, t + MODULO));
            }

            return Self.fromInteger(@intCast(u256, t));
        }

        pub fn div(self: Self, den: Self) ?Self {
            const denInv = den.inv() orelse return null;
            return self.mul(denInv);
        }

        pub fn eq(self: Self, other: Self) bool {
            return std.mem.eql(u64, &self.fe, &other.fe);
        }

        pub fn toInteger(self: Self) u256 {
            var nonMont: gen_fp.NonMontgomeryDomainFieldElement = undefined;
            gen_fp.fromMontgomery(&nonMont, self.fe);

            var bytes: [BYTE_LEN]u8 = [_]u8{0} ** BYTE_LEN;
            gen_fp.toBytes(&bytes, nonMont);

            return std.mem.readInt(u256, &bytes, std.builtin.Endian.Little);
        }

        // TODO
        // def modular_sqrt(a, p):
        //     """ Find a quadratic residue (mod p) of 'a'. p
        //         must be an odd prime.
        //         Solve the congruence of the form:
        //             x^2 = a (mod p)
        //         And returns x. Note that p - x is also a root.
        //         0 is returned is no square root exists for
        //         these a and p.
        //         The Tonelli-Shanks algorithm is used (except
        //         for some simple cases in which the solution
        //         is known from an identity). This algorithm
        //         runs in polynomial time (unless the
        //         generalized Riemann hypothesis is false).
        //     """
        //     # Simple cases
        //     #
        //     if legendre_symbol(a, p) != 1:
        //         return None
        //     elif a == 0:
        //         return 0
        //     elif p == 2:
        //         return 0
        //     elif p % 4 == 3:
        //         return pow(a, (p + 1) // 4, p)

        //     # Partition p-1 to s * 2^e for an odd s (i.e.
        //     # reduce all the powers of 2 from p-1)
        //     #
        //     s = p - 1
        //     e = 0
        //     while s % 2 == 0:
        //         s //= 2
        //         e += 1

        //     # Find some 'n' with a legendre symbol n|p = -1.
        //     # Shouldn't take long.
        //     #
        //     n = 2
        //     while legendre_symbol(n, p) != -1:
        //         n += 1

        //     # Here be dragons!
        //     # Read the paper "Square roots from 1; 24, 51,
        //     # 10 to Dan Shanks" by Ezra Brown for more
        //     # information
        //     #

        //     # x is a guess of the square root that gets better
        //     # with each iteration.
        //     # b is the "fudge factor" - by how much we're off
        //     # with the guess. The invariant x^2 = ab (mod p)
        //     # is maintained throughout the loop.
        //     # g is used for successive powers of n to update
        //     # both a and b
        //     # r is the exponent - decreases with each update
        //     #
        //     x = pow(a, (s + 1) // 2, p)
        //     b = pow(a, s, p)
        //     g = pow(n, s, p)
        //     r = e

        //     while True:
        //         t = b
        //         m = 0
        //         for m in range(r):
        //             if t == 1:
        //                 break
        //             t = pow(t, 2, p)

        //         if m == 0:
        //             return x

        //         gs = pow(g, 2 ** (r - m - 1), p)
        //         g = (gs * gs) % p
        //         x = (x * gs) % p
        //         b = (b * g) % p
        //         r = m

        // def legendre_symbol(a, p):
        //     """ Compute the Legendre symbol a|p using
        //         Euler's criterion. p is a prime, a is
        //         relatively prime to p (if p divides
        //         a, then a|p = 0)
        //         Returns 1 if a has a square root modulo
        //         p, -1 otherwise.
        //     """
        //     ls = pow(a, (p - 1) // 2, p)
        //     return -1 if ls == p - 1 else ls
    };
}

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
    try std.testing.expect(!Fp.fromInteger(Fp.Q_MIN_ONE_DIV_2).lexographically_largest());

    try std.testing.expect(Fp.fromInteger(Fp.Q_MIN_ONE_DIV_2 + 1).lexographically_largest());
    try std.testing.expect(Fp.fromInteger(Fp.MODULO - 1).lexographically_largest());
}

test "from and to bytes" {
    const cases = [_]Fp{ Fp.fromInteger(0), Fp.fromInteger(1), Fp.fromInteger(Fp.Q_MIN_ONE_DIV_2), Fp.fromInteger(Fp.MODULO - 1) };

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
    const got = Fp.fromInteger(10).mul(Fp.fromInteger(20)).add(Fp.fromInteger(30)).sub(Fp.fromInteger(40)).add(Fp.fromInteger(Fp.MODULO));
    const want = Fp.fromInteger(190);
    try std.testing.expect(got.eq(want));

    const gotneg = got.neg();
    const wantneg = Fp.fromInteger(Fp.MODULO - 190);
    try std.testing.expect(gotneg.eq(wantneg));
}

test "inv" {
    try std.testing.expect(Fp.fromInteger(0).inv() == null);

    const one = Fp.one();
    const cases = [_]Fp{ Fp.fromInteger(1), Fp.fromInteger(42), Fp.fromInteger(Fp.MODULO - 1) };
    for (cases) |fe| {
        try std.testing.expect(fe.mul(fe.inv().?).eq(one));
    }
}
