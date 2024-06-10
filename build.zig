const std = @import("std");

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const verkle_crypto_module = b.createModule(.{
        .root_source_file = .{ .cwd_relative = "src/main.zig" },
    });
    try b.modules.put(b.dupe("verkle-crypto"), verkle_crypto_module);

    const lib = b.addStaticLibrary(.{
        .name = "verkle-crypto",
        .root_source_file = .{ .cwd_relative = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(lib);

    const main_tests = b.addTest(.{
        .root_source_file = .{ .cwd_relative = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    const test_step = b.step("test", "Run library tests");
    const run_test = b.addRunArtifact(main_tests);
    // run_test.has_side_effects = true;
    test_step.dependOn(&run_test.step);

    const bench = b.addExecutable(.{
        .name = "bench",
        .root_source_file = .{ .cwd_relative = "src/bench.zig" },
        .target = target,
        .optimize = optimize,
    });

    const bench_step = b.step("bench", "Run benchmarks");
    const run_a = b.addRunArtifact(bench);
    bench_step.dependOn(&run_a.step);

    b.installArtifact(bench);
}
