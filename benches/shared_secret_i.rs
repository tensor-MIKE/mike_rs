#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

// I don't love this, but including the macro within src seems messy too
include!("bench_macros.rs");
use mike_rs::fields::MikeIFp2;
use mike_rs::mike::MIKE_I;
mike_secret_benches!(MIKE_I, MikeIFp2, "MIKE-I");
