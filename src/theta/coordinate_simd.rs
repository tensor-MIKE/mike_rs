//! Four-lane field operations on the coordinates of a dimension-four theta point.
//!
//! All lane positions and transform schedules are fixed, public constants.

use fp2::packed::PackedField;

use super::dim_four::{ThetaPoint, ThetaStructure};

mod chain;
mod codomain;
#[cfg(feature = "simd")]
pub(crate) use chain::try_finish;

pub(super) type Point<C> = [<C as PackedField>::Packed; 4];

#[inline(always)]
fn pack<C: PackedField>(ctx: &C, point: &ThetaPoint<C::Scalar>) -> Point<C> {
    core::array::from_fn(|block| ctx.pack(core::array::from_fn(|lane| point[4 * block + lane])))
}

#[inline(always)]
fn unpack<C: PackedField>(ctx: &C, point: &Point<C>) -> ThetaPoint<C::Scalar> {
    let coordinates = point.map(|x| ctx.unpack(&x));
    ThetaPoint::new_from_array(core::array::from_fn(|i| coordinates[i / 4][i % 4]))
}

/// Keep the full transform in a separate packed-context vectorize call. Combining
/// packs before permuting their lanes avoids extra live temporary packs.
#[inline(never)]
fn hadamard<C: PackedField>(ctx: &C, point: &mut Point<C>) {
    ctx.vectorize(
        #[inline(always)]
        || {
            // The four Hadamard dimensions commute: use distances 8, 4, 2, 1.
            for width in [2, 1] {
                for base in (0..4).step_by(2 * width) {
                    for i in base..base + width {
                        let a = point[i];
                        let b = point[i + width];
                        point[i] = ctx.add(&a, &b);
                        point[i + width] = ctx.sub(&a, &b);
                    }
                }
            }
            for x in point.iter_mut() {
                *x = ctx.butterfly::<2>(x);
                *x = ctx.butterfly::<1>(x);
            }
        },
    );
}

#[inline(always)]
fn square<C: PackedField>(ctx: &C, point: &mut Point<C>) {
    for x in point {
        *x = ctx.square(x);
    }
}

#[inline(always)]
fn multiply<C: PackedField>(ctx: &C, point: &mut Point<C>, other: &Point<C>) {
    for (x, y) in point.iter_mut().zip(other) {
        *x = ctx.mul(x, y);
    }
}

#[inline(always)]
fn double<C: PackedField>(
    ctx: &C,
    point: &mut Point<C>,
    inv_null: &Point<C>,
    inv_dual_sq: &Point<C>,
    count: usize,
) {
    for _ in 0..count {
        square(ctx, point);
        hadamard(ctx, point);
        square(ctx, point);
        multiply(ctx, point, inv_dual_sq);
        hadamard(ctx, point);
        multiply(ctx, point, inv_null);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::marker::PhantomData;
    use fp2::packed::{PackedField, PackedOperation};
    use fp2::traits::Fq;

    pub(super) fn random_point<F: Fq>(state: &mut u64) -> ThetaPoint<F> {
        ThetaPoint::new_from_array(core::array::from_fn(|_| {
            let mut bytes = F::ZERO.encode();
            for byte in bytes.as_mut() {
                *state ^= *state << 13;
                *state ^= *state >> 7;
                *state ^= *state << 17;
                *byte = *state as u8;
            }
            F::decode_reduce(bytes.as_ref())
        }))
    }

    struct PointOperation<F>(PhantomData<F>);

    impl<F: Fq> PackedOperation<F> for PointOperation<F> {
        type Output = ();

        fn run<C: PackedField<Scalar = F>>(self, context: C) {
            check_points_with_context(context);
        }
    }

    fn check_points_with_context<F: Fq, C: PackedField<Scalar = F>>(ctx: C) {
        ctx.vectorize(
            #[inline(always)]
            || {
                let mut state = 0x636f_6f72_6473_696du64;
                for round in 0..32 {
                    let point = match round {
                        0 => ThetaPoint::default(),
                        1 => ThetaPoint::new_from_array(core::array::from_fn(|i| {
                            [F::ZERO, F::ONE, F::MINUS_ONE, F::MINUS_ONE - F::ONE][i % 4]
                        })),
                        _ => random_point::<F>(&mut state),
                    };
                    let packed = pack(&ctx, &point);
                    let roundtrip = unpack(&ctx, &packed);
                    for i in 0..16 {
                        assert_eq!(roundtrip[i].equals(&point[i]), u32::MAX, "roundtrip {i}");
                    }
                    let mut transformed = packed;
                    hadamard(&ctx, &mut transformed);
                    let transformed = unpack(&ctx, &transformed);
                    let expected = point.hadamard();
                    for i in 0..16 {
                        assert_eq!(
                            transformed[i].equals(&expected[i]),
                            u32::MAX,
                            "Hadamard {i}"
                        );
                    }
                    let domain = ThetaStructure::new_from_null_point(point);
                    let (a, b) = domain.doubling_constants();
                    let mut actual = packed;
                    double(&ctx, &mut actual, &pack(&ctx, &a), &pack(&ctx, &b), 3);
                    let actual = unpack(&ctx, &actual);
                    let expected = domain.double_iter(&point, 3);
                    assert_eq!(actual.equals(&expected), u32::MAX, "doubling");
                    let mut actual = packed;
                    square(&ctx, &mut actual);
                    hadamard(&ctx, &mut actual);
                    multiply(&ctx, &mut actual, &pack(&ctx, &a));
                    hadamard(&ctx, &mut actual);
                    let mut expected = point.square().hadamard();
                    expected.set_coordinate_multiply(&a);
                    expected.set_hadamard();
                    assert_eq!(
                        unpack(&ctx, &actual).equals(&expected),
                        u32::MAX,
                        "evaluation"
                    );
                }
            },
        );
    }

    fn check_points<F: Fq>() {
        F::with_packed(PointOperation(PhantomData));
    }

    #[test]
    fn point_operations_match_scalar() {
        check_points::<crate::fields::MikeIFp>();
        check_points::<crate::fields::MikeIIIFp>();
        check_points::<crate::fields::MikeVFp>();
    }
}
