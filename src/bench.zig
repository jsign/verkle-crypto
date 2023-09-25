const std = @import("std");
const fields = @import("fields/fields.zig");
const crs = @import("crs/crs.zig");
const Fp = fields.BandersnatchFields.BaseField;
const banderwagon = @import("banderwagon/banderwagon.zig");
const Fr = banderwagon.Fr;
const multiproof = @import("multiproof/multiproof.zig");
const polynomials = @import("polynomial/lagrange_basis.zig");
const ipa = @import("ipa/ipa.zig");
const Transcript = @import("ipa/transcript.zig");
const precomp = @import("crs/msm.zig");

pub fn main() !void {
    try benchFields();
    try benchPedersenHash();
    try benchIPAs();
    try benchMultiproofs();
}

fn benchFields() !void {
    std.debug.print("Setting up fields benchmark...\n", .{});
    const N = 150_000;
    const set_size = 30_000;
    var fps: [set_size]Fp = genBaseFieldElements(set_size);
    var start: i64 = undefined;
    var startNano: i128 = undefined;

    std.debug.print("\tLegendre symbol... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.legendre(fps[i % set_size]);
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tField square root... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        _ = Fp.sqrt(fps[i % set_size]).?;
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tField inverse... ", .{});
    start = std.time.microTimestamp();
    for (0..N) |i| {
        fps[i % set_size] = Fp.inv(fps[i % set_size]).?;
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tField batch inverse (100 elements)... ", .{});
    start = std.time.microTimestamp();
    var inv_fps: [100]Fp = undefined;
    for (0..N) |i| {
        const start_idx = i % (set_size - 100);
        try Fp.batchInv(&inv_fps, fps[start_idx .. start_idx + 100]);
    }
    std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});

    std.debug.print("\tMul... ", .{});
    startNano = std.time.nanoTimestamp();
    for (0..N) |i| {
        fps[i % set_size] = Fp.mul(fps[i % set_size], fps[(i + 1) % set_size]);
    }
    std.debug.print("takes {}ns\n", .{@divTrunc((std.time.nanoTimestamp() - startNano), (N))});

    std.debug.print("\tAdd... ", .{});
    startNano = std.time.nanoTimestamp();
    for (0..N) |i| {
        fps[i % set_size] = Fp.add(fps[i % set_size], fps[(i + 1) % set_size]);
    }
    std.debug.print("takes {}ns\n", .{@divTrunc((std.time.nanoTimestamp() - startNano), (N))});
}

fn benchPedersenHash() !void {
    std.debug.print("Benchmarking Pedersen hashing...\n", .{});
    const xcrs = crs.CRS.init();
    const N = 5000;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) std.testing.expect(false) catch @panic("memory leak");
    }
    var allocator = gpa.allocator();

    var precomp_msm = try precomp.PrecompMSM(4, 6, crs.DomainSize).init(allocator, xcrs.Gs);
    defer precomp_msm.deinit();

    var vec_len: usize = 1;
    while (vec_len <= 256) : (vec_len <<= 1) {
        std.debug.print("\twith {} elements... ", .{vec_len});

        var vecs = try allocator.alloc([crs.DomainSize]Fr, N);
        defer allocator.free(vecs);
        for (0..vecs.len) |i| {
            for (0..vec_len) |j| {
                vecs[i][j] = Fr.fromInteger(i + j + 0x424242);
            }
            for (vec_len..vecs[i].len) |j| {
                vecs[i][j] = Fr.zero();
            }
        }

        var start = std.time.microTimestamp();
        for (0..N) |i| {
            _ = try precomp_msm.msm(vecs[i][0..vec_len]);
        }
        std.debug.print("takes {}µs\n", .{@divTrunc((std.time.microTimestamp() - start), (N))});
    }
}

fn benchIPAs() !void {
    const PrecomputedWeights = @import("polynomial/precomputed_weights.zig").PrecomputedWeights(crs.DomainSize, crs.Domain);

    std.debug.print("Setting up IPA benchmark...\n", .{});
    const N = 100;

    var weights = PrecomputedWeights.init();
    const xcrs = crs.CRS.init();
    const IPA = ipa.IPA(crs.DomainSize);

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer {
        const deinit_status = gpa.deinit();
        if (deinit_status == .leak) std.testing.expect(false) catch @panic("memory leak");
    }
    var allocator = gpa.allocator();

    var prover_queries: []IPA.ProverQuery = try allocator.alloc(IPA.ProverQuery, 16);
    defer allocator.free(prover_queries);
    const z256 = Fr.fromInteger(256);
    for (0..prover_queries.len) |i| {
        for (0..prover_queries[i].A.len) |j| {
            prover_queries[i].A[j] = Fr.fromInteger(i + j + 0x424242);
        }
        prover_queries[i].commitment = xcrs.commit(prover_queries[i].A);
        prover_queries[i].eval_point = Fr.fromInteger(i + 0x414039).add(z256);
        prover_queries[i].B = weights.barycentricFormulaConstants(prover_queries[i].eval_point);
    }

    var accum_prover: i64 = 0;
    var accum_verifier: i64 = 0;
    for (0..N) |j| {
        const i = j % prover_queries.len;
        // Prover.
        var prover_transcript = Transcript.init("test");
        var start = std.time.milliTimestamp();
        const proof = IPA.createProof(xcrs, &prover_transcript, prover_queries[i]);
        accum_prover += std.time.milliTimestamp() - start;

        // Verifier.
        start = std.time.milliTimestamp();
        var verifier_transcript = Transcript.init("test");
        const verifier_query = IPA.VerifierQuery{
            .commitment = prover_queries[i].commitment,
            .B = prover_queries[i].B,
            .eval_point = prover_queries[i].eval_point,
            .result = proof.result,
            .proof = proof.proof,
        };
        const ok = IPA.verifyProof(xcrs, &verifier_transcript, verifier_query);
        std.debug.assert(ok);
        accum_verifier += std.time.milliTimestamp() - start;
    }
    std.debug.print(
        "\tproving takes {}ms, verifying takes {}ms\n",
        .{
            @divTrunc((accum_prover), (N)),
            @divTrunc((accum_verifier), (N)),
        },
    );
}

fn benchMultiproofs() !void {
    const LagrangeBasis = polynomials.LagrangeBasis(crs.DomainSize, crs.Domain);

    std.debug.print("Setting up multiproofs benchmark...\n", .{});
    const N = 25;
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

    var vec_openings = try allocator.alloc(PolyOpeningSetup, openings[openings.len - 1]);
    defer allocator.free(vec_openings);

    for (0..vec_openings.len) |i| {
        for (0..vec_openings[i].poly_evaluations.len) |j| {
            vec_openings[i].poly_evaluations[j] = Fr.fromInteger(i + j + 0x424242);
        }
        vec_openings[i].z = Fr.fromInteger((i + 0x414039) % 256);
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
    fps[0] = Fp.fromInteger(0x4242 * 0x4242);
    fps[1] = Fp.fromInteger(0x4140 * 0x4140);
    var i: usize = 2;
    while (i < N) : (i += 1) {
        var fe = fps[i - 1].mul(fps[i - 2]);
        while (Fp.sqrt(fe) == null) {
            fe = fe.add(Fp.one());
        }
        fps[i] = fe;
    }
    return fps;
}
