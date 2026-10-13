//! Four independent MIKE theta chains evaluated in packed field lanes.
//!
//! Each packed value contains one coordinate from each of four independent
//! chains. The chain strategy and all loop bounds are public; scalar
//! multiplication and the gluing isogeny remain outside this module.

use crate::mike::isogeny_chain::MikeChainState;
use crate::theta::dim_four::{ThetaPoint, ThetaStructure};
use fp2::packed::{PackedField, PackedOperation};
use fp2::traits::Fq;

const COMPRESS: [usize; 10] = [0, 1, 2, 3, 5, 6, 7, 9, 11, 15];
const EXPAND: [usize; 16] = [0, 1, 2, 3, 2, 4, 5, 6, 1, 7, 4, 8, 3, 8, 6, 9];

type Point<C> = [<C as PackedField>::Packed; 16];

#[inline(always)]
fn pack<C: PackedField>(ctx: &C, points: &[ThetaPoint<C::Scalar>; 4]) -> Point<C> {
    core::array::from_fn(|coordinate| {
        ctx.pack(core::array::from_fn(|lane| points[lane][coordinate]))
    })
}

#[inline(always)]
fn unpack<C: PackedField>(ctx: &C, point: &Point<C>) -> [ThetaPoint<C::Scalar>; 4] {
    let coordinates: [[C::Scalar; 4]; 16] = core::array::from_fn(|i| ctx.unpack(&point[i]));
    core::array::from_fn(|lane| {
        ThetaPoint::new_from_array(core::array::from_fn(|i| coordinates[i][lane]))
    })
}

#[inline(always)]
fn hadamard<C: PackedField>(ctx: &C, point: &mut Point<C>) {
    for width in [1, 2, 4, 8] {
        for start in (0..16).step_by(2 * width) {
            for i in start..start + width {
                let a = point[i];
                let b = point[i + width];
                point[i] = ctx.add(&a, &b);
                point[i + width] = ctx.sub(&a, &b);
            }
        }
    }
}

// Only null points have these symmetries. Kernel and image points retain all
// sixteen coordinates while they are evaluated by the chain.
#[inline(always)]
fn null_square<C: PackedField>(ctx: &C, point: &Point<C>) -> Point<C> {
    let squared = COMPRESS.map(|i| ctx.square(&point[i]));
    EXPAND.map(|i| squared[i])
}

#[inline(always)]
fn null_hadamard<C: PackedField>(ctx: &C, point: &mut Point<C>) {
    let [z0, z1, z2, z3, z4, z5, z6, z7, z8, z9] = COMPRESS.map(|i| point[i]);
    let t0 = ctx.add(&z0, &z9);
    let t1 = ctx.sub(&z0, &z9);
    let t2 = ctx.add(&z5, &z7);
    let t3 = ctx.sub(&z5, &z7);
    let t4 = ctx.add(&z2, &z8);
    let t5 = ctx.sub(&z2, &z8);
    let t6 = ctx.add(&z1, &z6);
    let t7 = ctx.sub(&z1, &z6);
    let t8 = ctx.add(&z3, &z4);
    let t9 = ctx.sub(&z3, &z4);
    let s0 = ctx.add(&t0, &t2);
    let s1 = ctx.sub(&t0, &t2);
    let s2 = ctx.add(&t1, &t3);
    let s3 = ctx.sub(&t1, &t3);
    let s4 = ctx.add(&t6, &t4);
    let s5 = ctx.sub(&t6, &t4);
    let r0 = ctx.add(&t5, &t5);
    let r1 = ctx.add(&t7, &t7);
    let r2 = ctx.add(&t9, &t9);
    let r3 = ctx.add(&t8, &t8);
    let r4 = ctx.add(&s4, &s4);
    let r5 = ctx.add(&s5, &s5);
    let u0 = ctx.add(&s0, &r3);
    let u1 = ctx.sub(&s0, &r3);
    let result = [
        ctx.add(&u0, &r4),
        ctx.add(&s2, &r0),
        ctx.add(&s3, &r1),
        ctx.add(&s1, &r2),
        ctx.sub(&s1, &r2),
        ctx.add(&u1, &r5),
        ctx.sub(&s2, &r0),
        ctx.sub(&u1, &r5),
        ctx.sub(&s3, &r1),
        ctx.sub(&u0, &r4),
    ];
    *point = EXPAND.map(|i| result[i]);
}

#[inline(always)]
fn square<C: PackedField>(ctx: &C, point: &mut Point<C>) {
    for coordinate in point {
        *coordinate = ctx.square(coordinate);
    }
}

#[inline(always)]
fn multiply<C: PackedField>(ctx: &C, point: &mut Point<C>, other: &Point<C>) {
    for (coordinate, other_coordinate) in point.iter_mut().zip(other) {
        *coordinate = ctx.mul(coordinate, other_coordinate);
    }
}

#[inline(always)]
fn pseudo_inverse<C: PackedField>(ctx: &C, point: &Point<C>) -> Point<C> {
    let values = COMPRESS.map(|i| point[i]);
    let mut prefix = [ctx.pack([C::Scalar::ONE; 4]); 10];
    let mut suffix = [ctx.pack([C::Scalar::ONE; 4]); 10];
    prefix[0] = values[0];
    for i in 1..10 {
        prefix[i] = ctx.mul(&prefix[i - 1], &values[i]);
    }
    for i in 1..10 {
        suffix[i] = ctx.mul(&suffix[i - 1], &values[10 - i]);
    }
    let mut inverse = [ctx.pack([C::Scalar::ONE; 4]); 10];
    inverse[0] = suffix[9];
    for i in 1..10 {
        inverse[i] = ctx.mul(&suffix[9 - i], &prefix[i - 1]);
    }
    EXPAND.map(|i| inverse[i])
}

struct Domain<C: PackedField> {
    null: Point<C>,
    inv_null: Point<C>,
    inv_dual_sq: Point<C>,
}

#[inline(always)]
fn domain<C: PackedField>(ctx: &C, null: Point<C>) -> Domain<C> {
    let inv_null = pseudo_inverse(ctx, &null);
    let mut dual_sq = null_square(ctx, &null);
    null_hadamard(ctx, &mut dual_sq);
    let inv_dual_sq = pseudo_inverse(ctx, &dual_sq);
    Domain {
        null,
        inv_null,
        inv_dual_sq,
    }
}

#[inline(always)]
fn double<C: PackedField>(ctx: &C, domain: &Domain<C>, point: &mut Point<C>, count: usize) {
    for _ in 0..count {
        square(ctx, point);
        hadamard(ctx, point);
        square(ctx, point);
        multiply(ctx, point, &domain.inv_dual_sq);
        hadamard(ctx, point);
        multiply(ctx, point, &domain.inv_null);
    }
}

#[inline(always)]
fn isogeny<C: PackedField>(
    ctx: &C,
    mut p: Point<C>,
    mut q: Point<C>,
    images: &mut [Point<C>],
) -> Point<C> {
    square(ctx, &mut p);
    hadamard(ctx, &mut p);
    square(ctx, &mut q);
    hadamard(ctx, &mut q);

    let r1 = p[12];
    let r2 = ctx.mul(&r1, &p[5]);
    let r3 = ctx.mul(&r2, &q[1]);
    let r4 = ctx.mul(&r3, &p[9]);
    let r5 = ctx.mul(&r4, &p[11]);
    let r6 = ctx.mul(&r5, &q[15]);
    let r7 = ctx.mul(&r6, &q[14]);
    let r8 = ctx.mul(&r7, &p[6]);
    let b8 = p[0];
    let b7 = ctx.mul(&b8, &p[2]);
    let b6 = ctx.mul(&b7, &q[6]);
    let b5 = ctx.mul(&b6, &q[7]);
    let b4 = ctx.mul(&b5, &p[15]);
    let b3 = ctx.mul(&b4, &p[13]);
    let b2 = ctx.mul(&b3, &q[9]);
    let b1 = ctx.mul(&b2, &p[1]);

    let mut inv = [ctx.pack([C::Scalar::ZERO; 4]); 16];
    inv[0] = ctx.mul(&r8, &p[4]);
    inv[12] = ctx.mul(&ctx.mul(&p[8], &p[5]), &b2);
    inv[5] = ctx.mul(&r1, &b1);
    inv[1] = ctx.mul(&r2, &b2);
    inv[9] = ctx.mul(&r3, &b3);
    inv[11] = ctx.mul(&r4, &b4);
    inv[15] = ctx.mul(&r5, &b5);
    inv[14] = ctx.mul(&r6, &b6);
    inv[6] = ctx.mul(&r7, &b7);
    inv[4] = ctx.mul(&r8, &b8);
    inv[3] = inv[12];
    inv[10] = inv[5];
    inv[8] = inv[1];
    inv[13] = inv[11];
    inv[7] = inv[14];
    inv[2] = inv[4];

    for image in images {
        square(ctx, image);
        hadamard(ctx, image);
        multiply(ctx, image, &inv);
        hadamard(ctx, image);
    }
    inv
}

struct FinishOperation<F: Fq> {
    states: [MikeChainState<F>; 4],
    len: usize,
}

impl<F: Fq> PackedOperation<F> for FinishOperation<F> {
    type Output = [ThetaStructure<F>; 4];

    fn run<C: PackedField<Scalar = F>>(self, context: C) -> Self::Output {
        if !context.is_accelerated() {
            return self.states.map(|state| state.finish(self.len));
        }
        finish(&context, self.states, self.len)
    }
}

/// Finish four independent chains while keeping corresponding coordinates in
/// packed lanes. The selected fp2 backend may be serial or SIMD.
pub(crate) fn finish_four<F: Fq>(
    states: [MikeChainState<F>; 4],
    len: usize,
) -> [ThetaStructure<F>; 4] {
    F::with_packed(FinishOperation { states, len })
}

#[inline(always)]
fn finish<C: PackedField>(
    ctx: &C,
    states: [MikeChainState<C::Scalar>; 4],
    len: usize,
) -> [ThetaStructure<C::Scalar>; 4] {
    debug_assert!(states.iter().all(|state| state.orders == states[0].orders
        && state.k == states[0].k
        && state.theta_kernel.len() == states[0].theta_kernel.len()));

    ctx.vectorize(
        #[inline(always)]
        || {
            let nulls = states.each_ref().map(|state| state.domain.null_point());
            let mut packed_domain = domain(ctx, pack(ctx, &nulls));
            let mut points: Vec<Point<C>> = (0..states[0].theta_kernel.len())
                .map(|index| {
                    let theta_points = states.each_ref().map(|state| state.theta_kernel[index]);
                    pack(ctx, &theta_points)
                })
                .collect();
            let mut orders = states[0].orders.clone();
            let mut k = states[0].k;

            for step in 0..len {
                while orders[k] != 1 {
                    k += 1;
                    let m = orders[k - 1] >> 1;
                    points[2 * k] = points[2 * k - 2];
                    points[2 * k + 1] = points[2 * k - 1];
                    double(ctx, &packed_domain, &mut points[2 * k], m);
                    double(ctx, &packed_domain, &mut points[2 * k + 1], m);
                    orders[k] = orders[k - 1] - m;
                }

                let t1 = points[2 * k];
                let t2 = points[2 * k + 1];
                let inv = isogeny(ctx, t1, t2, &mut points[..2 * k]);
                if step == len - 1 || (k > 0 && orders[k - 1] > 2) {
                    let mut null = pseudo_inverse(ctx, &inv);
                    null_hadamard(ctx, &mut null);
                    packed_domain = domain(ctx, null);
                }
                for order in orders.iter_mut().take(k) {
                    *order -= 1;
                }
                k = k.saturating_sub(1);
            }

            unpack(ctx, &packed_domain.null).map(ThetaStructure::new_from_null_point)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::{MikeIFp, MikeIIIFp, MikeVFp};
    use fp2::traits::Fq;

    struct ArithmeticOperation<F>(core::marker::PhantomData<F>);

    impl<F: Fq> PackedOperation<F> for ArithmeticOperation<F> {
        type Output = ();

        fn run<C: PackedField<Scalar = F>>(self, ctx: C) {
            ctx.vectorize(|| {
                let values: [F; 4] = core::array::from_fn(|lane| {
                    let mut bytes = F::ZERO.encode();
                    bytes.as_mut()[0] = lane as u8 + 1;
                    F::decode_reduce(bytes.as_ref())
                });
                let point = ThetaPoint::new_from_array(core::array::from_fn(|i| values[i & 3]));
                let packed = pack(&ctx, &[point; 4]);
                let mut transformed = packed;
                hadamard(&ctx, &mut transformed);
                let unpacked = unpack(&ctx, &transformed);
                let expected = point.hadamard();
                for actual in unpacked {
                    assert_eq!(actual.equals(&expected), u32::MAX);
                }
                let mut squared = packed;
                square(&ctx, &mut squared);
                let expected = point.square();
                for actual in unpack(&ctx, &squared) {
                    assert_eq!(actual.equals(&expected), u32::MAX);
                }

                let null = ThetaPoint::new_from_array([F::ONE; 16]);
                let packed_null = pack(&ctx, &[null; 4]);
                let actual_null = unpack(&ctx, &null_square(&ctx, &packed_null));
                let expected_null = null.square();
                for actual in actual_null {
                    assert_eq!(actual.equals(&expected_null), u32::MAX);
                }
            });
        }
    }

    #[test]
    fn packed_theta_arithmetic_matches_scalar() {
        MikeIFp::with_packed(ArithmeticOperation(core::marker::PhantomData));
        MikeIIIFp::with_packed(ArithmeticOperation(core::marker::PhantomData));
        MikeVFp::with_packed(ArithmeticOperation(core::marker::PhantomData));
    }

    #[test]
    fn four_chain_finish_matches_independent_finishes_and_scalar() {
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                use crate::fields::MikeTestFp;
                use crate::theta::test_data::{OA_NULL, P, Q};

                fn state(len: usize, lane: usize) -> MikeChainState<MikeTestFp> {
                    let space = (usize::BITS - len.leading_zeros() + 1) as usize;
                    let scale = [
                        MikeTestFp::ONE,
                        MikeTestFp::TWO,
                        MikeTestFp::MINUS_ONE,
                        MikeTestFp::MINUS_ONE - MikeTestFp::ONE,
                    ][lane];
                    let scale_point = ThetaPoint::new_from_array([scale; 16]);
                    let mut theta_kernel = vec![ThetaPoint::default(); 2 * space];
                    theta_kernel[0] = P.coordinate_multiply(&scale_point);
                    theta_kernel[1] = Q.coordinate_multiply(&scale_point);
                    let mut orders = vec![0; space];
                    orders[0] = len;
                    MikeChainState {
                        domain: ThetaStructure::new_from_null_point(OA_NULL),
                        theta_kernel,
                        orders,
                        k: 0,
                    }
                }

                // Keep this direct differential test below the large packed stack frame;
                // full-size MIKE chains are covered through the public batch tests.
                for len in [0usize, 1, 2, 3, 7, 8, 15, 16, 31] {
                    let states = [state(len, 0), state(len, 1), state(len, 2), state(len, 3)];
                    let independent = states.clone().map(|state| state.finish(len));
                    let packed = finish_four(states, len);
                    let scalar = state(len, 0).domain.dim_four_isogeny_chain(&P, &Q, len);
                    for (lane, (actual, expected)) in packed.iter().zip(independent).enumerate() {
                        assert_eq!(
                            actual.null_point().equals(&expected.null_point()),
                            u32::MAX,
                            "independent lane {lane}, length {len}"
                        );
                        assert_eq!(
                            actual.null_point().equals(&scalar.null_point()),
                            u32::MAX,
                            "scalar lane {lane}, length {len}"
                        );
                    }
                }
            })
            .expect("large-stack packed differential test thread")
            .join()
            .expect("packed differential test panicked");
    }
}
