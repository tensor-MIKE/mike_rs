use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use mike_rs::mike::{MIKE_I, MIKE_III, MIKE_V};
use std::hint::black_box;

// Generate independent, reproducible keys outside the timed scan. Agreement
// benchmarks use already-prepared peers; preparation is measured separately.
macro_rules! batch_secret_benches {
    ($function:ident, $instance:expr, $name:expr) => {
        fn $function(c: &mut Criterion) {
            let receiver_keys: Vec<_> = (0u8..8)
                .map(|index| {
                    let mut seed = [0u8; 48];
                    seed[0] = index;
                    $instance.keygen_seeded(&seed).unwrap().1
                })
                .collect();
            let peers: Vec<_> = (0u8..32)
                .map(|index| {
                    let mut seed = [0u8; 48];
                    seed[0] = index;
                    $instance.keygen_seeded(&seed).unwrap().0
                })
                .collect();
            let prepared: Vec<_> = peers
                .iter()
                .map(|peer| $instance.prepare_peer(peer).unwrap())
                .collect();

            let mut preparation = c.benchmark_group(concat!($name, " prepare_peer"));
            preparation.sample_size(10);
            for &peer_count in &[1usize, 2, 4, 8, 32] {
                let peer_keys = &peers[..peer_count];
                preparation.throughput(Throughput::Elements(peer_count as u64));
                preparation.bench_function(BenchmarkId::new("serial", peer_count), |b| {
                    b.iter(|| {
                        for peer in black_box(peer_keys) {
                            black_box($instance.prepare_peer(peer).unwrap());
                        }
                    })
                });
            }
            preparation.finish();

            let mut agreement = c.benchmark_group(concat!($name, " prepared agreement"));
            // Category 5, particularly with conservative parameters, has a
            // long per-sample runtime at the larger batch sizes.
            agreement.sample_size(10);
            for &peer_count in &[1usize, 2, 4, 8, 32] {
                for &key_count in &[1usize, 2, 4, 8] {
                    let keys = &receiver_keys[..key_count];
                    let prepared_peers = &prepared[..peer_count];
                    agreement.throughput(Throughput::Elements((peer_count * key_count) as u64));
                    let size = format!("{peer_count} peers, {key_count} keys");
                    agreement.bench_function(BenchmarkId::new("serial", &size), |b| {
                        b.iter(|| {
                            for peer in black_box(prepared_peers) {
                                for key in black_box(keys) {
                                    black_box(key.shared_secret_from_prepared(peer).unwrap());
                                }
                            }
                        })
                    });
                }
            }
            agreement.finish();
        }
    };
}

batch_secret_benches!(bench_i, MIKE_I, "MIKE-I");
batch_secret_benches!(bench_iii, MIKE_III, "MIKE-III");
batch_secret_benches!(bench_v, MIKE_V, "MIKE-V");
criterion_group!(benches, bench_i, bench_iii, bench_v);
criterion_main!(benches);
