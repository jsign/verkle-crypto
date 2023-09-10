const std = @import("std");
const Bandersnatch = @import("bandersnatch.zig");
const Fp = Bandersnatch.Fp;
const Fr = Bandersnatch.Fr;
const AffinePoint = Bandersnatch.AffinePoint;
const ExtendedPoint = Bandersnatch.ExtendedPoint;
const ArrayList = std.ArrayList;
const allocator_test = std.testing.allocator;

// TODO(jsign): move to separate module.
pub const Banderwagon = struct {
    point: ExtendedPoint,

    // TODO(jsign): 32 -> const.
    pub fn initUnsafe(serlialized: [32]u8) Banderwagon {
        return Banderwagon{ .point = ExtendedPoint.initUnsafe(serlialized) };
    }

    // TODO(jsign): 32 -> const.
    pub fn fromBytes(serialised_bytes_big_endian: [32]u8) !Banderwagon {
        var bytes_le: [32]u8 = undefined;
        std.mem.copy(u8, &bytes_le, &serialised_bytes_big_endian);
        std.mem.reverse(u8, &bytes_le);

        // TODO(jsign): Should -> "Will error if the bytes are not canonical"
        const x = Fp.fromBytes(bytes_le);

        // Will error if the point is not on the curve
        const y = try AffinePoint.getYCoordinate(x, true);

        // Will error if x coordinate not a quadratic residue
        if (subgroupCheck(x) != 1) {
            return error.NotInSubgroup;
        }

        return Banderwagon{ .point = ExtendedPoint.initUnsafe(x, y) };
    }

    pub fn eq(self: *const Banderwagon, other: *const Banderwagon) bool {
        // The equals method is different for the quotient group
        //
        // Check for the (0,0) point, which is _possible_
        // given that you do not need to use the constructor to construct points
        const x1 = self.point.x;
        const y1 = self.point.y;
        const x2 = other.point.x;
        const y2 = other.point.y;

        if (x1.isZero() and y1.isZero()) {
            return false;
        }
        if (x2.isZero() and y2.isZero()) {
            return false;
        }

        const lhs = Fp.mul(x1, y2);
        const rhs = Fp.mul(x2, y1);

        return Fp.eq(lhs, rhs);
    }

    pub fn generator() Banderwagon {
        return .{ .point = ExtendedPoint.generator() };
    }

    pub fn neg(self: Banderwagon, p: Banderwagon) Banderwagon {
        self.point = -p.point;
        return self;
    }

    pub fn add(self: *Banderwagon, p: *const Banderwagon, q: *const Banderwagon) void {
        self.point = ExtendedPoint.add(p.point, q.point);
    }

    pub fn sub(self: Banderwagon, p: Banderwagon, q: Banderwagon) Banderwagon {
        self.point = p.point - q.point;
        return self;
    }

    pub fn mapToFieldBytes(self: Banderwagon) [32]u8 {
        return self.map_to_field().to_bytes();
    }

    pub fn map_to_field(self: Banderwagon) Fp {
        // The map to field function for banderwagon is x/y
        const x = self.point.x.dup();
        const y = self.point.y.dup();

        return x / y;
    }

    pub fn subgroupCheck(x: Fp) i2 {
        // Compute 1 - aX^2 and check its legendre symbol
        var res = x.mul(x);
        res = res.mul(Bandersnatch.A);
        res = res.neg();
        res = res.add(Fp.one());

        return res.legendre();
    }

    pub fn toBytes(self: Banderwagon) [32]u8 {
        const affine = self.point.toAffine();
        var x = affine.x;
        if (!affine.y.lexographicallyLargest()) {
            x = Fp.neg(x);
        }

        // Little endian.
        var bytes = x.toBytes();
        // Big endian.
        std.mem.reverse(u8, &bytes);

        return bytes;
    }

    pub fn double(self: *Banderwagon, p: Banderwagon) void {
        self.point = p.point.double();
    }

    pub fn isOnCurve(self: Banderwagon) bool {
        return self.point.toAffine().isOnCurve();
    }

    pub fn dup(self: Banderwagon) Banderwagon {
        return Banderwagon{
            .point = self.point,
        };
    }

    // TODO(jsign): weird api.
    pub fn scalarMul(element: Banderwagon, scalar: Fr) Banderwagon {
        return Banderwagon{
            .point = ExtendedPoint.scalarMul(element.point, scalar),
        };
    }

    pub fn identity() Banderwagon {
        return Banderwagon{ .point = ExtendedPoint.identity() };
    }

    pub fn twoTorsionPoint() Banderwagon {
        const point = ExtendedPoint.init(Fp.zero(), Fp.one().neg()) catch unreachable;
        return Banderwagon{ .point = point };
    }

    // Multi scalar multiplication
    pub fn msm(points: []const Banderwagon, scalars: []const Fr) Banderwagon {
        var res = Banderwagon.identity();

        for (scalars, points) |scalar, point| {
            const partial_res = point.scalarMul(scalar);
            res.add(&res, &partial_res);
        }
        return res;
    }
};

test "serialize smoke" {
    // Each successive point is a doubling of the previous one
    // The first point is the generator
    const expected_bit_strings = [_][]const u8{
        "4a2c7486fd924882bf02c6908de395122843e3e05264d7991e18e7985dad51e9",
        "43aa74ef706605705989e8fd38df46873b7eae5921fbed115ac9d937399ce4d5",
        "5e5f550494159f38aa54d2ed7f11a7e93e4968617990445cc93ac8e59808c126",
        "0e7e3748db7c5c999a7bcd93d71d671f1f40090423792266f94cb27ca43fce5c",
        "14ddaa48820cb6523b9ae5fe9fe257cbbd1f3d598a28e670a40da5d1159d864a",
        "6989d1c82b2d05c74b62fb0fbdf8843adae62ff720d370e209a7b84e14548a7d",
        "26b8df6fa414bf348a3dc780ea53b70303ce49f3369212dec6fbe4b349b832bf",
        "37e46072db18f038f2cc7d3d5b5d1374c0eb86ca46f869d6a95fc2fb092c0d35",
        "2c1ce64f26e1c772282a6633fac7ca73067ae820637ce348bb2c8477d228dc7d",
        "297ab0f5a8336a7a4e2657ad7a33a66e360fb6e50812d4be3326fab73d6cee07",
        "5b285811efa7a965bd6ef5632151ebf399115fcc8f5b9b8083415ce533cc39ce",
        "1f939fa2fd457b3effb82b25d3fe8ab965f54015f108f8c09d67e696294ab626",
        "3088dcb4d3f4bacd706487648b239e0be3072ed2059d981fe04ce6525af6f1b8",
        "35fbc386a16d0227ff8673bc3760ad6b11009f749bb82d4facaea67f58fc60ed",
        "00f29b4f3255e318438f0a31e058e4c081085426adb0479f14c64985d0b956e0",
        "3fa4384b2fa0ecc3c0582223602921daaa893a97b64bdf94dcaa504e8b7b9e5f",
    };
    var points = try ArrayList(Banderwagon).initCapacity(allocator_test, expected_bit_strings.len);
    defer points.deinit();
    var point = Banderwagon.generator();

    // Check that encoding algorithm gives expected results
    for (expected_bit_strings) |bit_string| {
        const byts = std.fmt.bytesToHex(point.toBytes(), std.fmt.Case.lower);
        try std.testing.expectEqualSlices(u8, bit_string, &byts);

        points.appendAssumeCapacity(point);
        point.double(point);
    }

    // Check that decoding algorithm is correct
    for (expected_bit_strings, 0..) |bit_string, i| {
        const expected_point = points.items[i];

        var byts: [32]u8 = undefined;
        _ = try std.fmt.hexToBytes(&byts, bit_string);
        const decoded_point = try Banderwagon.fromBytes(byts);
        try std.testing.expect(decoded_point.eq(&expected_point));
    }
}

test "two torsion" {
    // two points which differ by the order two point (0,-1) should be
    // considered the same
    const gen = Banderwagon.generator();
    const two_torsion = Banderwagon.twoTorsionPoint();

    var result = Banderwagon.identity();
    result.add(&gen, &two_torsion);

    try std.testing.expect(result.eq(&gen));
}
