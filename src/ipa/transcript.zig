const std = @import("std");
const sha256 = std.crypto.hash.sha2.Sha256;
const banderwagon = @import("../banderwagon/banderwagon.zig");
const Fr = banderwagon.Fr;
const Element = banderwagon.Element;

state: sha256,

const Transcript = @This();

pub fn init(label: []const u8) Transcript {
    var state = sha256.init(.{});
    state.update(label);
    return .{
        .state = state,
    };
}
// Convert bytes to scalar field
fn bytesToField(bytes: [Fr.BytesSize]u8) Fr {
    return Fr.fromBytes(bytes);
}

fn appendBytes(self: *Transcript, message: []const u8, label: []const u8) void {
    self.state.update(label);
    self.state.update(message);
}

pub fn appendScalar(self: *Transcript, scalar: Fr, label: []const u8) void {
    // Serialize the scalar in little endian
    const bytes = scalar.toBytes();
    self.appendBytes(&bytes, label);
}

pub fn appendPoint(self: *Transcript, point: Element, label: []const u8) void {
    const point_as_bytes = point.toBytes();
    self.appendBytes(&point_as_bytes, label);
}

pub fn challengeScalar(self: *Transcript, label: []const u8) Fr {
    self.domainSep(label);

    // hash the transcript to produce the challenge
    const hash = self.state.finalResult();
    const challenge = bytesToField(hash);

    // Clear the sha256 state
    // This step is not completely necessary
    // This is done so it frees memory
    self.state = sha256.init(.{});

    // Add the produced challenge into the new state
    // This is done for two reasons:
    // - It is now impossible for protocols using this
    //   class to forget to add any challenges previously seen
    // - It is now secure to repeatedly call for challenges
    self.appendScalar(challenge, label);

    // Return the new challenge
    return challenge;
}

// domainSep is used to:
// - Separate between adding elements to the transcript and squeezing elements out
// - Separate sub-protocols
pub fn domainSep(self: *Transcript, label: []const u8) void {
    self.state.update(label);
}

//     def test_prover_verifier_consistency(self):
//         """
//             Test that if the prover and verifier
//             do the exact same operations, they will end up
//             at the exact same state.
//         """
//         point = Banderwagon.generator()
//         scalar = Fr(randint(0, 2**256))

//         # Prover's View
//         prover_transcript = Transcript(b"protocol_name")

//         prover_transcript.append_point(point, b"D")
//         prover_transcript.domain_sep(b"sub_protocol_name")
//         prover_transcript.append_scalar(scalar, b"r")

//         # Prover challenge
//         prover_q = prover_transcript.challenge_scalar(b"q")

//         # Verifier's View
//         verifier_transcript = Transcript(b"protocol_name")

//         verifier_transcript.append_point(point, b"D")
//         verifier_transcript.domain_sep(b"sub_protocol_name")
//         verifier_transcript.append_scalar(scalar, b"r")

//         # Verifier challenge
//         verifier_q = verifier_transcript.challenge_scalar(b"q")

//         assert prover_q == verifier_q

test "test vector" {
    // Test that squeezing out a challenge twice
    // will produce different challenges. ie it is
    // not possible to accidentally generate the same challenge
    var transcript = Transcript.init("foo");
    const first_challenge = transcript.challengeScalar("f");
    const second_challenge = transcript.challengeScalar("f");

    try std.testing.expect(!first_challenge.equal(second_challenge));
}

test "test vector 1" {
    // Test that challenge creation is consistent across implementations
    var transcript = Transcript.init("simple_protocol");
    const challenge = transcript.challengeScalar("simple_challenge");

    const got_hex = std.fmt.bytesToHex(challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("c2aa02607cbdf5595f00ee0dd94a2bbff0bed6a2bf8452ada9011eadb538d003", &got_hex);
}

test "test vector 2" {
    // Test that append scalar is consistent across implementations
    var transcript = Transcript.init("simple_protocol");
    const scalar = Fr.fromInteger(5);

    transcript.appendScalar(scalar, "five");
    transcript.appendScalar(scalar, "five again");

    const challenge = transcript.challengeScalar("simple_challenge");
    const got_hex = std.fmt.bytesToHex(challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("498732b694a8ae1622d4a9347535be589e4aee6999ffc0181d13fe9e4d037b0b", &got_hex);
}

test "test vector 3" {
    // Test that domain separation is consistent across implementations
    var transcript = Transcript.init("simple_protocol");
    const minus_one = Fr.fromInteger(Fr.MODULO - 1);
    const one = Fr.one();
    transcript.appendScalar(minus_one, "-1");
    transcript.domainSep("separate me");
    transcript.appendScalar(minus_one, "-1 again");
    transcript.domainSep("separate me again");
    transcript.appendScalar(one, "now 1");

    const challenge = transcript.challengeScalar("simple_challenge");

    const got_hex = std.fmt.bytesToHex(challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("14f59938e9e9b1389e74311a464f45d3d88d8ac96adf1c1129ac466de088d618", &got_hex);
}

test "test vector 4" {
    // Test that appending points is consistent across implementations
    var transcript = Transcript.init("simple_protocol");

    const generator = Element.generator();

    transcript.appendPoint(generator, "generator");
    const challenge = transcript.challengeScalar("simple_challenge");
    const got_hex = std.fmt.bytesToHex(challenge.toBytes(), std.fmt.Case.lower);
    try std.testing.expectEqualStrings("8c2dafe7c0aabfa9ed542bb2cbf0568399ae794fc44fdfd7dff6cc0e6144921c", &got_hex);
}
