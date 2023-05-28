const Bandersnatch = @import("bandersnatch.zig");
const Fp = Bandersnatch.Fp;
const Fr = Bandersnatch.Fr;
const AffinePoint = Bandersnatch.AffinePoint;

const ExtendedPoint = struct {
    x: Fp,
    y: Fp,
    t: Fp,
    z: Fp,

    pub fn init(x: Fp, y: Fp) !ExtendedPoint {
        if (!AffinePoint.isOnCurve(.{ .x = x, .y = y })) {
            return Bandersnatch.CurveError.NotInCurve;
        }
        return initUnsafe(x, y);
    }

    pub fn initUnsafe(x: Fp, y: Fp) ExtendedPoint {
        return ExtendedPoint{
            .x = x,
            .y = y,
            .t = x.mul(y),
            .z = Fp.one(),
        };
    }

    pub fn identity() ExtendedPoint {
        const iden = comptime AffinePoint.identity();
        return comptime initUnsafe(iden.x, iden.y);
    }

    pub fn generator() ExtendedPoint {
        const gen = comptime AffinePoint.generator();
        return comptime initUnsafe(gen.x, gen.y);
    }

    pub fn neg(self: ExtendedPoint) ExtendedPoint {
        return ExtendedPoint{
            .x = self.x.neg(),
            .y = self.y,
            .t = self.t.neg(),
            .z = self.z,
        };
    }

    pub fn isZero(self: ExtendedPoint) bool {
        // Identity is {x=0, y=1, t = 0, z =1}
        // The equivalence class is therefore is {x=0, y=k, t = 0, z=k} for all k where k!=0
        const condition_1 = self.x.is_zero();
        const condition_2 = self.y == self.z;
        const condition_3 = !self.y.is_zero();
        const condition_4 = self.t.is_zero();

        return condition_1 and condition_2 and condition_3 and condition_4;
    }

    pub fn equal(p: ExtendedPoint, q: ExtendedPoint) bool {
        if (p.is_zero()) {
            return q.is_zero();
        }

        if (q.is_zero()) {
            return false;
        }

        return (p.x.mul(q.z) == p.z.mul(q.x)) and (p.y.mul(q.z) == q.y.mul(p.z));
    }

    pub fn add(p: ExtendedPoint, q: ExtendedPoint) ExtendedPoint {
        // See "Twisted Edwards Curves Revisited" (https: // eprint.iacr.org/2008/522.pdf)
        // by Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, and Ed Dawson
        // 3.1 Unified Addition in E^e

        const x1 = p.x;
        const y1 = p.y;
        const t1 = p.t;
        const z1 = p.z;

        const x2 = q.x;
        const y2 = q.y;
        const t2 = q.t;
        const z2 = q.z;

        const a = x1 * x2;

        const b = y1 * y2;

        const c = Bandersnatch.D * t1 * t2;

        const d = z1 * z2;

        const h = b - (a * Bandersnatch.A);

        const e = (x1 + y1) * (x2 + y2) - a - b;

        const f = d - c;
        const g = d + c;

        return ExtendedPoint{
            .x = e * f,
            .y = g * h,
            .t = e * h,
            .z = f * g,
        };
    }

    pub fn sub(p: ExtendedPoint, q: ExtendedPoint) ExtendedPoint {
        const neg_q = q.neg();
        return add(p, neg_q);
    }

    pub fn double(self: ExtendedPoint) ExtendedPoint {
        // TODO: can replace this with dedicated doubling formula
        return add(self, self);
    }

    pub fn scalarMul(point: ExtendedPoint, scalar: Fr) ExtendedPoint {
        // Same as AffinePoint's equivalent method
        // using double and add : https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#Double-and-add
        const result = identity();
        const temp = point;

        for (scalar.fe) |limb| {
            for (0..@bitSizeOf(scalar.fe[0])) |i| {
                if (scalar.fe[limb] & 1 << (limb * i) == 1) {
                    result.add(result, temp);
                }
                temp.double(temp);
            }
        }
        return result;
    }

    pub fn toAffine(self: ExtendedPoint) AffinePoint {
        if (self.is_zero()) {
            return identity();
        } else if (self.z.is_one()) {
            return initUnsafe(self.x, self.y);
        } else {
            const z_inv = self.z.inv();

            const x_aff = self.x * z_inv;
            const y_aff = self.y * z_inv;
            return AffinePoint.initUnsafe(x_aff, y_aff);
        }
    }

    // # Only used for testing purposes.
    pub fn toBytes(self: ExtendedPoint) [32]u8 {
        self.toAffine().toBytes();
    }
};
