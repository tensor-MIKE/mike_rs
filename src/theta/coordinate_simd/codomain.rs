use super::*;

#[cfg(test)]
use super::tests::random_point;

type CompressedPoint<C> = [<C as PackedField>::Packed; 3];

#[inline(always)]
fn compress<C: PackedField>(ctx: &C, point: &Point<C>) -> CompressedPoint<C> {
    [
        point[0],
        ctx.gather::<5, 6, 7, 9>(point),
        ctx.gather::<11, 15, 11, 15>(point),
    ]
}

#[inline(always)]
fn decompress<C: PackedField>(ctx: &C, point: &CompressedPoint<C>) -> Point<C> {
    [
        point[0],
        ctx.gather::<2, 4, 5, 6>(point),
        ctx.gather::<1, 7, 4, 8>(point),
        ctx.gather::<3, 8, 6, 9>(point),
    ]
}

/// Products of the other nine compressed coordinates, with no inversion.
#[inline(always)]
fn pseudo_inverse<C: PackedField>(ctx: &C, point: &CompressedPoint<C>) -> CompressedPoint<C> {
    // Product trees for groups of 4, 4, and 2 coordinates. The last two
    // lanes are duplicates rather than padding by a field representation of 1.
    let pairs = ctx.mul(
        &ctx.gather::<0, 2, 4, 6>(point),
        &ctx.gather::<1, 3, 5, 7>(point),
    );
    let values = [pairs, point[0], point[1], point[2]];
    let totals = ctx.mul(
        &ctx.gather::<0, 2, 12, 12>(&values),
        &ctx.gather::<1, 3, 13, 13>(&values),
    );
    let cross = ctx.mul(
        &ctx.gather::<0, 0, 1, 1>(&[totals]),
        &ctx.gather::<1, 2, 2, 2>(&[totals]),
    );
    let a_partial = ctx.mul(
        &ctx.gather::<1, 0, 3, 2>(point),
        &ctx.gather::<1, 1, 0, 0>(&[pairs]),
    );
    let a = ctx.mul(&a_partial, &ctx.gather::<2, 2, 2, 2>(&[cross]));
    let b_partial = ctx.mul(
        &ctx.gather::<5, 4, 7, 6>(point),
        &ctx.gather::<3, 3, 2, 2>(&[pairs]),
    );
    let b = ctx.mul(&b_partial, &ctx.gather::<1, 1, 1, 1>(&[cross]));
    let c = ctx.mul(
        &ctx.gather::<0, 0, 0, 0>(&[cross]),
        &ctx.gather::<9, 8, 9, 8>(point),
    );
    [a, b, c]
}

#[inline(always)]
fn compressed_hadamard<C: PackedField>(ctx: &C, point: &CompressedPoint<C>) -> CompressedPoint<C> {
    let mut full = decompress(ctx, point);
    super::hadamard(ctx, &mut full);
    compress(ctx, &full)
}

#[inline(always)]
pub(super) fn doubling_constants<C: PackedField>(
    ctx: &C,
    inv_null_dual: &Point<C>,
) -> (Point<C>, Point<C>) {
    let null_dual = pseudo_inverse(ctx, &compress(ctx, inv_null_dual));
    let null = compressed_hadamard(ctx, &null_dual);
    let inv_null = decompress(ctx, &pseudo_inverse(ctx, &null));
    let mut dual_sq = null;
    for x in &mut dual_sq {
        *x = ctx.square(x);
    }
    let dual_sq = compressed_hadamard(ctx, &dual_sq);
    let inv_dual_sq = decompress(ctx, &pseudo_inverse(ctx, &dual_sq));
    (inv_null, inv_dual_sq)
}

#[inline(always)]
pub(super) fn inv_dual_codomain<C: PackedField>(ctx: &C, t1: &Point<C>, t2: &Point<C>) -> Point<C> {
    let mut p = *t1;
    let mut q = *t2;
    super::square(ctx, &mut p);
    super::square(ctx, &mut q);
    super::hadamard(ctx, &mut p);
    super::hadamard(ctx, &mut q);

    // Eight multiplication batches implement the scalar prefix/suffix formula.
    // Every output is homogeneous of degree 9 and has P/Q bidegree (6, 3).
    let s0 = ctx.mul(
        &ctx.gather::<12, 0, 8, 1>(&p),
        &ctx.gather::<5, 2, 5, 12>(&p),
    );
    let s1 = ctx.mul(
        &ctx.gather::<0, 1, 0, 1>(&[s0]),
        &ctx.gather::<1, 6, 1, 6>(&q),
    );
    let s2 = ctx.mul(
        &ctx.gather::<0, 1, 0, 1>(&[s1]),
        &ctx.gather::<1, 7, 1, 7>(&[p[2], q[1]]),
    );
    let s3 = ctx.mul(
        &ctx.gather::<0, 1, 0, 1>(&[s2]),
        &ctx.gather::<11, 15, 11, 15>(&p),
    );
    let s4 = ctx.mul(
        &ctx.gather::<0, 1, 4, 0>(&[s3, s2]),
        &ctx.gather::<3, 5, 9, 13>(&[q[3], p[3], s3, s2]),
    );
    let s5 = ctx.mul(
        &ctx.gather::<0, 1, 4, 0>(&[s4, s1]),
        &ctx.gather::<2, 5, 9, 13>(&[q[3], q[2], s4, s1]),
    );
    let s6 = ctx.mul(
        &ctx.gather::<0, 1, 4, 0>(&[s5, s0]),
        &ctx.gather::<2, 5, 9, 13>(&[p[1], p[0], s5, s0]),
    );
    let s7 = ctx.mul(
        &ctx.gather::<0, 6, 8, 0>(&[s6, s0, p[3]]),
        &ctx.gather::<0, 5, 9, 12>(&[p[1], s5, s6, p[0]]),
    );
    let out = [s4, s5, s6, s7];
    [
        ctx.gather::<12, 10, 15, 13>(&out),
        ctx.gather::<15, 14, 11, 7>(&out),
        ctx.gather::<10, 6, 14, 2>(&out),
        ctx.gather::<13, 2, 7, 3>(&out),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::theta::dim_four::ThetaPointCompressed;
    use core::marker::PhantomData;
    use fp2::packed::{PackedField, PackedOperation};
    use fp2::traits::Fq;

    struct CodomainOperation<F>(PhantomData<F>);

    impl<F: Fq> PackedOperation<F> for CodomainOperation<F> {
        type Output = ();

        fn run<C: PackedField<Scalar = F>>(self, context: C) {
            check_with_context(context);
        }
    }

    fn check_with_context<F: Fq, C: PackedField<Scalar = F>>(ctx: C) {
        ctx.vectorize(
            #[inline(always)]
            || {
                let mut state = 0x636f_6f72_6473_696du64;
                for seed in 0..32u32 {
                    let point = match seed {
                        0 => ThetaPoint::new_from_array(core::array::from_fn(|i| {
                            [F::ZERO, F::ONE, F::MINUS_ONE, F::MINUS_ONE - F::ONE][i % 4]
                        })),
                        1..=10 => {
                            let zero_at = [0, 1, 2, 3, 5, 6, 7, 9, 11, 15][seed as usize - 1];
                            ThetaPoint::new_from_array(core::array::from_fn(|i| {
                                if i == zero_at { F::ZERO } else { F::MINUS_ONE }
                            }))
                        }
                        11 => ThetaPoint::default(),
                        _ => random_point::<F>(&mut state),
                    };
                    let other = if seed == 0 {
                        ThetaPoint::new_from_array(core::array::from_fn(|i| {
                            [F::ZERO, F::MINUS_ONE - F::ONE, F::ONE, F::MINUS_ONE][(i + 1) % 4]
                        }))
                    } else {
                        random_point::<F>(&mut state)
                    };

                    let expected_pseudo = ThetaPointCompressed::from_theta_point(&point)
                        .proj_batch_pseudo_inversion()
                        .to_theta_point();
                    let expected_pseudo_compressed =
                        ThetaPointCompressed::from_theta_point(&expected_pseudo);
                    let expected_inv =
                        ThetaStructure::two_isogeny_evaluate(&point, &other, &mut []);
                    let expected_domain = ThetaStructure::new_from_inv_null_dual(point);
                    let (expected_inv_null, expected_inv_dual_sq) =
                        expected_domain.doubling_constants();

                    let packed_point = pack(&ctx, &point);
                    let packed_other = pack(&ctx, &other);
                    let pseudo = pseudo_inverse(&ctx, &compress(&ctx, &packed_point));
                    let actual_pseudo = unpack(&ctx, &decompress(&ctx, &pseudo));
                    let actual_inv =
                        unpack(&ctx, &inv_dual_codomain(&ctx, &packed_point, &packed_other));
                    let constants = doubling_constants(&ctx, &packed_point);
                    let actual_constants = (unpack(&ctx, &constants.0), unpack(&ctx, &constants.1));

                    let actual_pseudo_compressed =
                        ThetaPointCompressed::from_theta_point(&actual_pseudo);
                    for i in 0..10 {
                        assert_eq!(
                            actual_pseudo_compressed[i].equals(&expected_pseudo_compressed[i]),
                            u32::MAX,
                            "pseudo inverse exact coordinate {i}, seed {seed}"
                        );
                    }
                    assert_eq!(
                        actual_pseudo.equals(&expected_pseudo),
                        u32::MAX,
                        "pseudo inverse projective seed {seed}"
                    );
                    assert_eq!(
                        actual_inv.equals(&expected_inv),
                        u32::MAX,
                        "inverse dual codomain seed {seed}"
                    );
                    assert_eq!(
                        actual_constants.0.equals(&expected_inv_null),
                        u32::MAX,
                        "inverse null seed {seed}"
                    );
                    assert_eq!(
                        actual_constants.1.equals(&expected_inv_dual_sq),
                        u32::MAX,
                        "inverse dual square seed {seed}"
                    );
                }
            },
        );
    }

    fn check<F: Fq>() {
        F::with_packed(CodomainOperation(PhantomData));
    }

    #[test]
    fn compressed_codomain_operations_match_scalar() {
        check::<crate::fields::MikeIFp>();
        check::<crate::fields::MikeIIIFp>();
        check::<crate::fields::MikeVFp>();
    }
}
