// I don't love this, but including the macro within src seems messy too
include!("bench_macros.rs");
use mike_rs::fields::MikeVFp2;
use mike_rs::mike::MIKE_V;
mike_secret_benches!(MIKE_V, MikeVFp2, "MIKE-V");
