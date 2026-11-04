use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;

use crate::elliptic::curve::Curve;
use crate::elliptic::jac_point::Point;
use crate::gluing::gluing_isogeny;
use crate::mike::mike_impl::KernelPoints;
use crate::theta::dim_four::{ThetaPoint, ThetaStructure};

/// Below this order the dim one strategy stops splitting and doubles all the
/// way down, as dim one doublings are much cheaper than dim four ones (and
/// every extra stack entry has to be pushed through the gluing). Tuned by
/// benchmarking 4, 8, ..., 128: 16 is fastest for MIKE-I and MIKE-V.
const DIM_ONE_SPLIT_THRESHOLD: usize = 16;

impl<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>> Curve<Fp2> {
    /// Advance the balanced strategy by doubling kernel points until
    /// the top of the stack has order 2. Returns the new stack depth.
    fn advance_strategy(
        &self,
        a_div_three: &Fp2,
        pts: &mut [[Point<Fp2>; 2]],
        orders: &mut [usize],
        k: usize,
    ) -> usize {
        let mut k = k;
        while orders[k] != 1 {
            k += 1;

            let m = if orders[k - 1] >= DIM_ONE_SPLIT_THRESHOLD {
                orders[k - 1] >> 1
            } else {
                orders[k - 1] - 1
            };

            let [p1, p2] = pts[2 * k - 2];
            let [p3, p4] = pts[2 * k - 1];
            pts[2 * k] = [
                self.double_iter_with_ws(&p1, a_div_three, m),
                self.double_iter_with_ws(&p2, a_div_three, m),
            ];
            pts[2 * k + 1] = [
                self.double_iter_with_ws(&p3, a_div_three, m),
                self.double_iter_with_ws(&p4, a_div_three, m),
            ];
            orders[k] = orders[k - 1].saturating_sub(m);
        }
        k
    }

    /// Compute the isogeny chain (E x E^sigma) -> A using a balanced strategy used for Mike
    /// secret generation.
    ///
    /// The kernels of the three gluing isogenies are derived from the top of
    /// the dim one strategy stack and passed to `check_gluing_kernel` before
    /// the gluing; an error from it is returned.
    pub fn mike_isogeny_chain<E>(
        &self,
        a_div_three: &Fp2,
        chain_kernel: &[[Point<Fp2>; 2]; 2],
        len: usize,
        check_gluing_kernel: impl FnOnce(&KernelPoints<Fp2>) -> Result<(), E>,
    ) -> Result<ThetaStructure<Fp>, E> {
        let n = len - 2;
        let space = (usize::BITS - n.leading_zeros() + 1) as usize;

        let mut kernel_points: Vec<[Point<Fp2>; 2]> = vec![[Point::INFINITY; 2]; 2 * space];
        let mut orders: Vec<usize> = vec![0; space];
        kernel_points[0] = chain_kernel[0];
        kernel_points[1] = chain_kernel[1];
        orders[0] = n;

        let mut k = self.advance_strategy(a_div_three, &mut kernel_points, &mut orders, 0);

        // The top of the stack is the chain kernel at the level of the gluing,
        // from which we derive the kernels of all three gluing isogenies
        let level = [kernel_points[2 * k], kernel_points[2 * k + 1]];
        let gluing_kernel = KernelPoints::from_gluing_level(self, &level);
        check_gluing_kernel(&gluing_kernel)?;

        // Compute the gluing isogeny and push through all the points we doubled for the strategy
        let (mut domain, mut theta_kernel) =
            gluing_isogeny(self, &gluing_kernel, &kernel_points[..(2 * k)]);
        // Gluing isogeny creates theta_kernel of length 2*k but we need it to have a minimum size
        theta_kernel.resize(2 * space, ThetaPoint::default());

        // Reduce the order of every point pushed through the isogeny
        for ord in orders.iter_mut().take(k) {
            *ord -= 1;
        }
        k = k.saturating_sub(1);

        // The three gluing steps (Scholten, diagonal, dim four gluing) are done,
        // so e - 3 (2,2,2,2)-isogenies remain. Running the loop longer only
        // recomputes the final codomain from the same exhausted kernel.
        for _ in 0..(len - 3) {
            // Perform doublings of the kernel elements, decreasing the values of orders
            k = domain.advance_strategy(&mut theta_kernel, &mut orders, k);

            // Extract out the kernel for this step.
            let T1 = theta_kernel[2 * k];
            let T2 = theta_kernel[2 * k + 1];

            // Perform one step of the (2,2) isogeny and push through all points.
            domain = ThetaStructure::two_isogeny(&T1, &T2, &mut theta_kernel[..(2 * k)]);

            // Reduce the order of the points we evaluated
            for ord in orders.iter_mut().take(k) {
                *ord -= 1;
            }
            k = k.saturating_sub(1);
        }

        Ok(domain)
    }
}
