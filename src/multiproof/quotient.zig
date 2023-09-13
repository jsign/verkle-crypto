const std = @import("std");
const Allocator = std.mem.Allocator;
const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights;
const LagrangeBasis = @import("../polynomial/lagrange_basis.zig").LagrangeBasis;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const Fr = Bandersnatch.Fr;

pub fn compute_quotient_inside_domain(
    allocator: Allocator,
    precomp: PrecomputedWeights,
    f: LagrangeBasis,
    index: Fr,
) !LagrangeBasis {
    const domain_size = precomp.domain.len;
    const inverses = precomp.domain_inverses;
    const Aprime_domain = precomp.Aprime_DOMAIN;
    const Aprime_domain_inv = precomp.Aprime_DOMAIN_inv;

    const indexU256 = index.toInteger();
    if (indexU256 >= domain_size) {
        return error.IndexOutOfDomain;
    }
    // This cast is safe since we checked above.
    const indexInt = @as(u8, @intCast(indexU256));

    var q = [_]Fr{Fr.zero()} ** 256; // TODO: 256 from comptime
    const y = f.evaluations.items[indexInt];
    for (0..domain_size) |i| {
        if (i != indexInt) {
            q[i] = Fr.mul(Fr.sub(f.evaluations.items[i], y), inverses[i - indexInt]);
            q[indexInt] = Fr.add(
                q[indexInt],
                Fr.mul(Fr.mul(Fr.mul(Fr.sub(f.evaluations.items[i], y), inverses[indexInt - i]), Aprime_domain[indexInt]), Aprime_domain_inv[i]),
            );
        }
    }

    return try LagrangeBasis.init(allocator, &q, f.domain.items);
}

pub fn compute_quotient_outside_domain(
    allocator: Allocator,
    precomp: PrecomputedWeights,
    f: LagrangeBasis,
    z: Fr,
    y: Fr,
) !LagrangeBasis {
    const domain = precomp.domain;
    const domain_size = domain.len; // TODO: comptime

    var q = [_]Fr{Fr.zero()} * 256; // TODO: comptime.
    for (0..domain_size) |i| {
        q[i] = Fr.mul(Fr.sub(f[i], y), Fr.inv(Fr.sub(domain[i], z)));
    }

    return try LagrangeBasis.init(allocator, q, f.domain);
}
