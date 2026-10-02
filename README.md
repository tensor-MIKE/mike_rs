# MIKE: Module Isogeny Key Exchange

An efficient and constant-time implementation of the isogeny-based NIKE, MIKE.

## Parameter Sets

This implementation includes parameter sets targeting NIST levels I, III and V for both the "FastMIKE" and "RigorousMIKE". By default, the more compact and efficient parameters are selected. To compile with the larger parameter set use
the feature flag `--features rigorous_parameters`

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



