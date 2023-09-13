const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const lib = b.addStaticLibrary(.{
        .name = "verkle-crypto",
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(lib);

    // _ = b.createModule(.{
    //     .source_file = .{ .path = "ecc/bandersnatch/bandersnatch.zig" },
    // });

    const main_tests = b.addTest(.{
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    const test_step = b.step("test", "Run library tests");
    const run_test = b.addRunArtifact(main_tests);
    run_test.has_side_effects = true;
    test_step.dependOn(&run_test.step);
}
