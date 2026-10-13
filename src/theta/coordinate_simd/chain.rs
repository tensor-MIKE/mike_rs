//! Evaluate the public theta-chain strategy without leaving packed arithmetic.

use fp2::packed::{PackedField, PackedOperation};
use fp2::traits::Fq;

use super::{
    ThetaPoint, ThetaStructure, codomain, double, hadamard, multiply, pack, square, unpack,
};

struct FinishOperation<'a, F: Fq> {
    domain: ThetaStructure<F>,
    points: &'a [ThetaPoint<F>],
    orders: &'a [usize],
    k: usize,
    len: usize,
}

impl<F: Fq> PackedOperation<F> for FinishOperation<'_, F> {
    type Output = Option<ThetaStructure<F>>;

    fn run<C: PackedField<Scalar = F>>(self, context: C) -> Self::Output {
        if !context.is_accelerated() {
            return None;
        }
        Some(finish(
            &context,
            self.domain,
            self.points,
            self.orders,
            self.k,
            self.len,
        ))
    }
}

/// Return `None` when the selected field has no accelerated packed backend.
/// The caller retains the existing scalar chain as its fallback.
pub(crate) fn try_finish<F: Fq>(
    domain: ThetaStructure<F>,
    points: &[ThetaPoint<F>],
    orders: &[usize],
    k: usize,
    len: usize,
) -> Option<ThetaStructure<F>> {
    F::with_packed(FinishOperation {
        domain,
        points,
        orders,
        k,
        len,
    })
}

#[inline(always)]
fn finish<C: PackedField>(
    ctx: &C,
    domain: ThetaStructure<C::Scalar>,
    points: &[ThetaPoint<C::Scalar>],
    orders: &[usize],
    mut k: usize,
    len: usize,
) -> ThetaStructure<C::Scalar> {
    if len == 0 {
        return domain;
    }

    ctx.vectorize(
        #[inline(always)]
        || {
            let mut points: Vec<_> = points.iter().map(|point| pack(ctx, point)).collect();
            let mut orders = orders.to_vec();
            let (inv_null, inv_dual) = domain.doubling_constants();
            let mut inv_null = pack(ctx, &inv_null);
            let mut inv_dual = pack(ctx, &inv_dual);
            let mut inv = points[0];

            for step in 0..len {
                while orders[k] != 1 {
                    k += 1;
                    let m = orders[k - 1] >> 1;
                    points[2 * k] = points[2 * k - 2];
                    points[2 * k + 1] = points[2 * k - 1];
                    double(ctx, &mut points[2 * k], &inv_null, &inv_dual, m);
                    double(ctx, &mut points[2 * k + 1], &inv_null, &inv_dual, m);
                    orders[k] = orders[k - 1] - m;
                }
                inv = codomain::inv_dual_codomain(ctx, &points[2 * k], &points[2 * k + 1]);
                for point in &mut points[..2 * k] {
                    square(ctx, point);
                    hadamard(ctx, point);
                    multiply(ctx, point, &inv);
                    hadamard(ctx, point);
                }

                // These decisions depend only on the public balanced strategy.
                // If the next kernel is ready, its evaluation needs no constants.
                if step != len - 1 && k > 0 && orders[k - 1] > 2 {
                    (inv_null, inv_dual) = codomain::doubling_constants(ctx, &inv);
                }
                for order in orders.iter_mut().take(k) {
                    *order -= 1;
                }
                k = k.saturating_sub(1);
            }

            ThetaStructure::new_from_inv_null_dual(unpack(ctx, &inv))
        },
    )
}

#[cfg(test)]
mod tests {
    use fp2::packed::{PackedField, PackedOperation};
    use fp2::traits::Fq;

    use super::*;
    use crate::fields::MikeTestFp;
    use crate::theta::test_data::{OA_NULL, P, Q};

    struct FinishTestOperation<'a, F: Fq> {
        domain: ThetaStructure<F>,
        points: &'a [ThetaPoint<F>],
        orders: &'a [usize],
        k: usize,
        len: usize,
    }

    impl<F: Fq> PackedOperation<F> for FinishTestOperation<'_, F> {
        type Output = ThetaStructure<F>;

        fn run<C: PackedField<Scalar = F>>(self, context: C) -> Self::Output {
            finish(
                &context,
                self.domain,
                self.points,
                self.orders,
                self.k,
                self.len,
            )
        }
    }

    #[test]
    fn unsupported_field_uses_scalar_fallback() {
        // p=127 is valid for scalar fp2 but has no all-ones low radix-32 digit.
        fp2::define_fp_core!(typename = SmallFp, modulus = [127u64],);
        let domain = ThetaStructure::<SmallFp>::new_from_null_point(ThetaPoint::default());
        assert!(try_finish(domain, &[], &[], 0, 0).is_none());
    }

    #[test]
    fn complete_chain_matches_scalar() {
        let domain = ThetaStructure::new_from_null_point(OA_NULL);
        let unchanged = MikeTestFp::with_packed(FinishTestOperation {
            domain,
            points: &[],
            orders: &[],
            k: 0,
            len: 0,
        });
        assert_eq!(
            unchanged.null_point().equals(&domain.null_point()),
            u32::MAX
        );

        for len in [1usize, 2, 3, 7, 8, 15, 16, 31, 63, 127, 241, 242] {
            let p = domain.double_iter(&P, 242 - len);
            let q = domain.double_iter(&Q, 242 - len);
            let space = (usize::BITS - len.leading_zeros() + 1) as usize;
            let mut points = vec![ThetaPoint::<MikeTestFp>::default(); 2 * space];
            points[0] = p;
            points[1] = q;
            let mut orders = vec![0; space];
            orders[0] = len;
            let dispatched = MikeTestFp::with_packed(FinishTestOperation {
                domain,
                points: &points,
                orders: &orders,
                k: 0,
                len,
            });
            let expected = domain.dim_four_isogeny_chain(&p, &q, len);
            assert_eq!(
                dispatched.null_point().equals(&expected.null_point()),
                u32::MAX,
                "length {len}"
            );
        }
    }
}
