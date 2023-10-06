# verkle-crypto

A pure Zig implementation of Verkle Tree cryptography. 

## Features

This library is feature complete:
- [X] Relevant finite field arithmetic.
- [X] Bandersnatch and Banderwagon group implementations.
- [X] Pedersen Hashing for trie-key caculations.
- [X] Proof creation and verification.

Some notes about the implementation:
- Both trie-key calculation and proof verification are very efficient.
- This library has no external dependencies.
- No assembly is used, so it can be compiled to all [supported targets](https://ziglang.org/download/0.11.0/release-notes.html#Tier-System).
- This library is single-threaded. It's planned to add multi-threading support in the future.
- Comptetitive with (or faster, single threaded) than [go-ipa](https://github.com/crate-crypto/go-ipa) or [rust-verkle](https://github.com/crate-crypto/rust-verkle).

This library isn't audited nor battle-tested, so it isn't recommended to be used in production.

## Test
```
$ zig build test --summary all
Build Summary: 3/3 steps succeeded; 48/48 tests passed
test success
└─ run test 48 passed 4s MaxRSS:344M
   └─ zig test ReleaseSafe native success 14s MaxRSS:388M

```

## Bench
`AMD Ryzen 7 3800XT`:
```
$ zig build bench -Dtarget=native -Doptimize=ReleaseFast         
Setting up fields benchmark...
        Legendre symbol... takes 9µs
        Field square root... takes 8µs
        Field inverse... takes 6µs
        Field batch inverse (100 elements)... takes 13µs
        Mul... takes 21ns
        Add... takes 5ns

Benchmarking Pedersen hashing...
        with 1 elements... takes 6µs
        with 2 elements... takes 7µs
        with 4 elements... takes 10µs
        with 8 elements... takes 15µs
        with 16 elements... takes 26µs
        with 32 elements... takes 46µs
        with 64 elements... takes 88µs
        with 128 elements... takes 170µs
        with 256 elements... takes 343µs

Setting up IPA benchmark...
        proving takes 55ms, verifying takes 4ms

Setting up multiproofs benchmark...
        Benchmarking 100 openings... proving takes 73ms, verifying takes 5ms
        Benchmarking 1000 openings... proving takes 120ms, verifying takes 12ms
        Benchmarking 5000 openings... proving takes 336ms, verifying takes 39ms
        Benchmarking 10000 openings... proving takes 607ms, verifying takes 72ms
```

## License

MIT.