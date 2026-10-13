# MIKE: Module Isogeny Key Exchange

[![CI](https://github.com/tensor-MIKE/mike_rs/actions/workflows/ci.yml/badge.svg)](https://github.com/tensor-MIKE/mike_rs/actions/workflows/ci.yml)

An efficient and constant-time implementation of the isogeny-based NIKE, MIKE.

## Parameter Sets

This implementation includes parameter sets targeting NIST levels I, III and V for both the "FastMIKE" and "RigorousMIKE". By default, the more compact and efficient parameters are selected. To compile with the larger parameter set use
the feature flag `--features conservative_parameters`

## Optional SIMD

The `simd` feature uses `fp2`'s four-lane field arithmetic within one key
agreement. It also evaluates compatible groups of four independent theta chains
in `batch_shared_secrets`, with one agreement per lane. Both paths request a
packed context through `Fq::with_packed`; when `PackedField::is_accelerated` is
false, they retain the existing scalar agreement path. SIMD requires Rust 1.89.

## Testing

Submodules as well as the key exchange have several tests, which can be run with the command:

```
cargo test
```

## Benchmarks

Benchmarks are computed using the [`criterion`](https://crates.io/crates/criterion) crate and are compiled and run using the command

```
cargo bench
```

The finite field arithmetic is computed using the [`fp2`](https://github.com/GiacomoPope/fp2/) crate, which for some CPU is optimized by passing the flag

```
RUSTFLAGS="-C target-cpu=native" cargo bench
```



