const std = @import("std");
const Fr = @import("../banderwagon/banderwagon.zig").Fr;

// TODO: Methods here may be moved into different modules in the future

pub fn innerProduct(a: []const Fr, b: []const Fr) Fr {
    var result = Fr.zero();
    for (a, b) |ai, bi| {
        const term = ai.mul(bi);
        result = Fr.add(result, term);
    }
    return result;
}

test "inner product smoke" {
    var a = [_]Fr{ Fr.fromInteger(1), Fr.fromInteger(2), Fr.fromInteger(3), Fr.fromInteger(4), Fr.fromInteger(5) };
    var b = [_]Fr{ Fr.fromInteger(10), Fr.fromInteger(12), Fr.fromInteger(13), Fr.fromInteger(14), Fr.fromInteger(15) };

    // Expected result should be 1*10 + 2*12 + 3*13 + 4*14 + 5*15
    const expected_result = Fr.fromInteger(204);
    const got_result = innerProduct(&a, &b);
    try std.testing.expect(got_result.eq(expected_result));
}
