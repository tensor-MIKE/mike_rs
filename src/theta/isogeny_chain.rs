use crate::theta::dim_four::{ThetaPoint, ThetaStructure};
use fp2::traits::Fq as FqTrait;

impl<Fq: FqTrait> ThetaStructure<Fq> {
    /// Advance the balanced strategy by doubling kernel points until
    /// the top of the stack has order 2^1. Returns the new stack depth.
    pub fn advance_strategy(
        &self,
        pts: &mut [ThetaPoint<Fq>],
        orders: &mut [usize],
        k: usize,
    ) -> usize {
        let mut k = k;
        while orders[k] != 1 {
            k += 1;
            let m = orders[k - 1] >> 1;
            pts[2 * k] = self.double_iter(&pts[2 * k - 2], m);
            pts[2 * k + 1] = self.double_iter(&pts[2 * k + 1 - 2], m);
            orders[k] = orders[k - 1].saturating_sub(m);
        }
        k
    }

    /// Compute a dim four isogeny assuming the domain is the output of the Mike gluing isogeny
    pub fn dim_four_isogeny_chain(
        &self,
        ker_1: &ThetaPoint<Fq>,
        ker_2: &ThetaPoint<Fq>,
        len: usize,
    ) -> Self {
        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - len.leading_zeros() + 1) as usize;

        // Store points of order 2^i for the balanced strategy. We need two
        // vectors here, as the first step computes with elements of type
        // ProductPoint, while every other step computes points of type
        // ThetaPoint.
        let mut kernel_pts: Vec<ThetaPoint<Fq>> = vec![ThetaPoint::default(); 2 * space];

        // The values i such that each point in stategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Include the kernel points into the strategy array
        kernel_pts[0] = *ker_1;
        kernel_pts[1] = *ker_2;

        // Initalise the orders list, points in the above vectors have order
        // 2^(orders[i] + 2), as we use the 8-torsion above.
        orders[0] = len;

        let mut domain = *self;
        let mut k = 0;
        for step in 0..len {
            // Perform doublings of the kernel elements, decreasing the values of orders
            k = domain.advance_strategy(&mut kernel_pts, &mut orders, k);

            // Extract out the kernel for this step.
            let T1 = kernel_pts[2 * k];
            let T2 = kernel_pts[2 * k + 1];

            // Perform one step of the (2,2) isogeny and push through all points.
            let inv_null_dual =
                ThetaStructure::two_isogeny_evaluate(&T1, &T2, &mut kernel_pts[..(2 * k)]);
            // Only the public strategy controls this decision. Consecutive
            // isogenies need no intervening domain construction.
            if step == len - 1 || (k > 0 && orders[k - 1] > 2) {
                domain = ThetaStructure::new_from_inv_null_dual(inv_null_dual);
            }

            // Reduce the order of the points we evaluated
            for ord in orders.iter_mut().take(k) {
                *ord -= 1;
            }
            k = k.saturating_sub(1);
        }

        domain
    }
}

#[cfg(test)]
mod tests {
    use super::ThetaStructure;
    use crate::theta::test_data::{OA_NULL, OB_NULL, P, Q};

    #[test]
    fn scheduled_domains_match_eager_chain() {
        use crate::theta::dim_four::ThetaPoint;

        let domain = ThetaStructure::new_from_null_point(OA_NULL);
        for len in [
            0usize, 1, 2, 3, 4, 5, 7, 8, 15, 16, 31, 32, 63, 64, 127, 128, 241, 242,
        ] {
            // Use kernels of the requested order from the full test chain.
            let p = domain.double_iter(&P, 242 - len);
            let q = domain.double_iter(&Q, 242 - len);
            let space = (usize::BITS - len.leading_zeros() + 1) as usize;
            let mut pts = vec![ThetaPoint::default(); 2 * space];
            let mut orders = vec![0; space];
            pts[0] = p;
            pts[1] = q;
            orders[0] = len;
            let mut eager = domain;
            let mut k = 0;
            for _ in 0..len {
                k = eager.advance_strategy(&mut pts, &mut orders, k);
                let t1 = pts[2 * k];
                let t2 = pts[2 * k + 1];
                eager = ThetaStructure::two_isogeny(&t1, &t2, &mut pts[..2 * k]);
                for order in orders.iter_mut().take(k) {
                    *order -= 1;
                }
                k = k.saturating_sub(1);
            }
            let actual = domain.dim_four_isogeny_chain(&p, &q, len);
            assert_eq!(
                actual.null_point().equals(&eager.null_point()),
                u32::MAX,
                "length {len}"
            );
        }
    }

    #[test]
    fn test_two_two_chain() {
        let oa = ThetaStructure::new_from_null_point(OA_NULL);
        let codomain = oa.dim_four_isogeny_chain(&P, &Q, 242);
        assert_eq!(codomain.null_point().equals(&OB_NULL), u32::MAX);
    }
}
