#[allow(unused)]
use criterion::{Criterion, criterion_group, criterion_main};
#[allow(unused)]
use mike_rs::gluing::gluing_isogeny;
#[allow(unused)]
use rand_core::OsRng;

#[allow(unused_macros)]
macro_rules! mike_keygen_benches {
    ($instance:expr, $prefix:expr) => {
        fn bench_keygen(c: &mut Criterion) {
            c.bench_function(concat!($prefix, " keygen"), |b| {
                b.iter(|| $instance.keygen(&mut OsRng).unwrap())
            });
        }

        fn bench_keygen_seeded(c: &mut Criterion) {
            let seed = [0u8; 48];
            c.bench_function(concat!($prefix, " keygen seeded"), |b| {
                b.iter(|| $instance.keygen_seeded(&seed))
            });
        }

        criterion_group!(benches, bench_keygen, bench_keygen_seeded);
        criterion_main!(benches);
    };
}

#[allow(unused_macros)]
macro_rules! mike_secret_benches {
    ($instance:expr, $fp2:ty, $prefix:expr) => {
        #[allow(dead_code)]
        fn bench_kernel_gen(c: &mut Criterion) {
            let seed_alice = [0u8; 48];
            let seed_bob = [1u8; 48];
            let (_, sk_a) = $instance.keygen_seeded(&seed_alice).unwrap();
            let (pk_b, _) = $instance.keygen_seeded(&seed_bob).unwrap();
            let domain = pk_b.curve;
            let a_div_three = domain.A / <$fp2>::THREE;

            c.bench_function(concat!($prefix, " kernel generation"), |b| {
                b.iter(|| sk_a.generate_kernel_elements(&domain, &a_div_three))
            });
        }

        #[allow(dead_code)]
        fn bench_gluing(c: &mut Criterion) {
            let seed_alice = [0u8; 48];
            let seed_bob = [1u8; 48];
            let (_, sk_a) = $instance.keygen_seeded(&seed_alice).unwrap();
            let (pk_b, _) = $instance.keygen_seeded(&seed_bob).unwrap();
            let domain = pk_b.curve;
            let a_div_three = domain.A / <$fp2>::THREE;
            let (kernel_data, image_points) = sk_a.generate_kernel_elements(&domain, &a_div_three);

            c.bench_function(concat!($prefix, " gluing isogeny"), |b| {
                b.iter(|| gluing_isogeny(&domain, &kernel_data, &image_points))
            });
        }

        fn bench_shared_secret(c: &mut Criterion) {
            let seed_alice = [0u8; 48];
            let seed_bob = [1u8; 48];
            let (_, sk_a) = $instance.keygen_seeded(&seed_alice).unwrap();
            let (pk_b, _) = $instance.keygen_seeded(&seed_bob).unwrap();

            c.bench_function(concat!($prefix, " shared secret"), |b| {
                b.iter(|| sk_a.shared_secret(&pk_b))
            });
        }

        // criterion_group!(benches, bench_kernel_gen, bench_gluing, bench_shared_secret);
        criterion_group!(benches, bench_shared_secret);
        criterion_main!(benches);
    };
}
