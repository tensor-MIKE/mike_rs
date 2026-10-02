#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

// I don't love this, but including the macro within src seems messy too
include!("bench_macros.rs");
use mike_rs::fields::MikeIIIFp2;
use mike_rs::mike::MIKE_III;
mike_secret_benches!(MIKE_III, MikeIIIFp2, "MIKE-III");
