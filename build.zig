const std = @import("std");

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("verkle-crypto", .{
        .root_source_file = b.path("src/main.zig"),
    });

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
