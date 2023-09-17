const std = @import("std");
const Allocator = std.mem.Allocator;
const PrecomputedWeights = @import("../polynomial/precomputed_weights.zig").PrecomputedWeights;
const LagrangeBasis = @import("../polynomial/lagrange_basis.zig").LagrangeBasis;
const Bandersnatch = @import("../ecc/bandersnatch/bandersnatch.zig");
const crs = @import("../crs/crs.zig");
const Fr = Bandersnatch.Fr;

pub fn compute_quotient_inside_domain(
    allocator: Allocator,
    precomp: PrecomputedWeights(crs.DomainSize, crs.Domain),
    f: LagrangeBasis,
    index: Fr,
) !LagrangeBasis {
    const inverses = precomp.domain_inverses;
    const Aprime_domain = precomp.Aprime_DOMAIN;
    const Aprime_domain_inv = precomp.Aprime_DOMAIN_inv;

    const indexU256 = index.toInteger();
    if (indexU256 >= crs.DomainSize) {
        return error.IndexOutOfDomain;
    }
    // This cast is safe since we checked above.
    const indexInt = @as(u8, @intCast(indexU256));

    var q = [_]Fr{Fr.zero()} ** crs.DomainSize;
    const y = f.evaluations.items[indexInt];
    for (0..crs.DomainSize) |i| {
        if (i != indexInt) {
            q[i] = Fr.mul(Fr.sub(f.evaluations.items[i], y), inverses[i - indexInt]);
            q[indexInt] = Fr.add(
                q[indexInt],
                Fr.mul(Fr.mul(Fr.mul(Fr.sub(f.evaluations.items[i], y), inverses[crs.DomainSize + i - indexInt]), Aprime_domain[indexInt]), Aprime_domain_inv[i]),
            );
        }
    }

    return try LagrangeBasis.init(allocator, &q, f.domain.items);
}

pub fn compute_quotient_outside_domain(
    allocator: Allocator,
    f: LagrangeBasis,
    z: Fr,
    y: Fr,
) !LagrangeBasis {
    var q = [_]Fr{Fr.zero()} * crs.DomainSize;
    for (0..crs.DomainSize) |i| {
        q[i] = Fr.mul(Fr.sub(f[i], y), Fr.inv(Fr.sub(crs.Domain[i], z)));
    }

    return try LagrangeBasis.init(allocator, q, f.domain);
}
