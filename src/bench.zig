const std = @import("std");
const fields = @import("fields/fields.zig");
const crs = @import("crs/crs.zig");
const Fp = fields.BandersnatchFields.BaseField;
const banderwagon = @import("banderwagon/banderwagon.zig");
const Fr = banderwagon.Fr;
const multiproof = @import("multiproof/multiproof.zig");
const polynomials = @import("polynomial/lagrange_basis.zig");
const Transcript = @import("ipa/transcript.zig");

const LagrangeBasis = polynomials.LagrangeBasis(crs.DomainSize, crs.Domain);

pub fn main() !void {
    benchFields();
    try benchMultiproofs();
}

fn benchFields() void {
    std.debug.print("Setting up fields benchmark...\n", .{});
    const N = 30_000;
    const fps: [N]Fp = genBaseFieldElements(N);
    var start: i64 = undefined;

    std.debug.print("\tLegendre symbol... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.legendre(fps[i]);
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tField inverse... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.inv(fps[i]).?;
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tField square root... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.sqrt(fps[i]).?;
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});
}

fn benchMultiproofs() !void {
    std.debug.print("Setting up multiproofs benchmark...\n", .{});
    const N = 50;
    const openings = [_]u16{ 1, 10, 100, 1_000 };

    const vkt_crs = crs.CRS.init();

    const PolyOpeningSetup = struct {
        poly_evaluations: [crs.DomainSize]Fr,
        C: banderwagon.Element,
        z: Fr,
    };

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) std.testing.expect(false) catch @panic("memory leak");
    }
    var allocator = gpa.allocator();

    var eval_magic_offset: u64 = 0x424242;
    var z_magic_offset: u64 = 0x414039;
    var vec_openings = try allocator.alloc(PolyOpeningSetup, openings[openings.len - 1]);
    defer allocator.free(vec_openings);

    for (0..vec_openings.len) |i| {
        for (0..vec_openings[i].poly_evaluations.len) |j| {
            vec_openings[i].poly_evaluations[j] = Fr.fromInteger(i + j + eval_magic_offset);
        }
        vec_openings[i].z = Fr.fromInteger((i + z_magic_offset) % 256);
        vec_openings[i].C = crs.CRS.commit(vkt_crs, vec_openings[i].poly_evaluations);
    }

    const mproof = multiproof.MultiProof.init(vkt_crs);
    for (openings) |num_openings| {
        std.debug.print("\tBenchmarking {} openings...", .{num_openings});

        var accum_proving: i64 = 0;
        var accum_verifying: i64 = 0;
        for (0..N) |_| {
            // Proving.
            var prover_queries = try allocator.alloc(multiproof.ProverQuery, num_openings);
            defer allocator.free(prover_queries);
            for (0..num_openings) |i| {
                prover_queries[i] = multiproof.ProverQuery{
                    .f = LagrangeBasis.init(vec_openings[i].poly_evaluations),
                    .C = vec_openings[i].C,
                    .z = vec_openings[i].z,
                    .y = vec_openings[i].poly_evaluations[@as(usize, @intCast(vec_openings[i].z.toInteger()))],
                };
            }

            var start = std.time.milliTimestamp();
            var prover_transcript = Transcript.init("test");
            const proof = try mproof.createProof(&prover_transcript, prover_queries);
            accum_proving += std.time.milliTimestamp() - start;

            // Verifying.
            var verifier_transcript = Transcript.init("test");
            var verifier_queries = try allocator.alloc(multiproof.VerifierQuery, num_openings);
            defer allocator.free(verifier_queries);
            for (0..num_openings) |i| {
                verifier_queries[i] = multiproof.VerifierQuery{
                    .C = vec_openings[i].C,
                    .z = vec_openings[i].z,
                    .y = vec_openings[i].poly_evaluations[@as(usize, @intCast(vec_openings[i].z.toInteger()))],
                };
            }
            start = std.time.milliTimestamp();
            const ok = try mproof.verifyProof(allocator, &verifier_transcript, verifier_queries, proof);
            std.debug.assert(ok);
            accum_verifying += std.time.milliTimestamp() - start;
        }
        std.debug.print(
            "proving takes {}ms, verifying takes {}ms\n",
            .{
                @divTrunc((accum_proving), (N)),
                @divTrunc((accum_verifying), (N)),
            },
        );
    }
}

fn genBaseFieldElements(comptime N: usize) [N]Fp {
    var fps: [N]Fp = undefined;
    var i: usize = 0;
    while (i < N) : (i += 1) {
        const fe = Fp.fromInteger(i);
        if (Fp.sqrt(fe) == null) {
            continue;
        }
        fps[i] = Fp.fromInteger(i);
    }
    return fps;
}
