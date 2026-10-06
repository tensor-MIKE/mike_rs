// I don't love this, but including the macro within src seems messy too
include!("bench_macros.rs");
use mike_rs::mike::MIKE_I;
mike_keygen_benches!(MIKE_I, "MIKE-I");
