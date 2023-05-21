const FieldApi = @import("fieldapi.zig");

// Bandersnatch base and scalar finite fields.
pub const Fp = FieldApi.Fp;
pub const Fr = FieldApi.Fr;

test "fields" {
    _ = @import("fieldapi.zig");
}
