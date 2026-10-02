#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

// I don't love this, but including the macro within src seems messy too
include!("bench_macros.rs");
use mike_rs::mike::MIKE_V;
mike_keygen_benches!(MIKE_V, "MIKE-V");
