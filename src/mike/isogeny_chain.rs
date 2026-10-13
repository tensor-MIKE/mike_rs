use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;

use crate::elliptic::curve::Curve;
use crate::elliptic::jac_point::Point;
use crate::gluing::gluing_isogeny;
use crate::mike::mike_impl::KernelPoints;
use crate::theta::dim_four::{ThetaPoint, ThetaStructure};

// The initial elliptic strategy is public and shared by preparation and agreement.
fn initial_strategy_split(order: usize) -> usize {
    if order >= 64 { order >> 1 } else { order - 1 }
}

/// Public basis multiples at the retained levels of the initial elliptic strategy.
#[derive(Clone, Debug)]
pub(crate) struct PreparedDoublings<Fp2: Fp2Trait>(Vec<[Point<Fp2>; 2]>);

impl<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>> Curve<Fp2> {
    /// Prepare the public multiples used to reconstruct paired kernel points.
    pub(crate) fn prepare_chain_doublings(
        &self,
        a_div_three: &Fp2,
        P: &Point<Fp2>,
        Q: &Point<Fp2>,
        len: usize,
    ) -> (PreparedDoublings<Fp2>, [Point<Fp2>; 2]) {
        let mut order = len - 2;
        let mut basis = [*P, *Q];
        let mut retained = Vec::new();
        while order != 1 {
            let m = initial_strategy_split(order);
            basis = basis.map(|p| self.double_iter_with_ws(&p, a_div_three, m));
            order -= m;
            if order != 1 {
                retained.push(basis);
            }
        }
        // In total the basis was doubled len - 3 times, so the final points
        // have order 32 and can also be used for the gluing kernels.
        (PreparedDoublings(retained), basis)
    }

    /// Advance the balanced strategy by doubling kernel points until
    /// the top of the stack has order 2. Returns the new stack depth.
    fn advance_strategy(
        &self,
        a_div_three: &Fp2,
        pts: &mut [[Point<Fp2>; 2]],
        orders: &mut [usize],
        k: usize,
        prepared: Option<&PreparedDoublings<Fp2>>,
    ) -> usize {
        let mut k = k;
        while orders[k] != 1 {
            k += 1;

            let m = initial_strategy_split(orders[k - 1]);

            orders[k] = orders[k - 1].saturating_sub(m);
            if orders[k] == 1 {
                // The gluing kernel is supplied separately. This last pair is
                // excluded from the gluing images, so no points are needed here.
                break;
            }

            let [p1, p2] = pts[2 * k - 2];
            let [p3, p4] = pts[2 * k - 1];
            if let Some(prepared) = prepared {
                // The paired points differ by public Q and -P. The table is
                // indexed only by the public strategy, never by scalar bits.
                let q = self.double_iter_with_ws(&p2, a_div_three, m);
                let minus_p = self.double_iter_with_ws(&p4, a_div_three, m);
                let [public_p, public_q] = prepared.0[k - 1];
                pts[2 * k] = [self.add(&q, &public_q), q];
                pts[2 * k + 1] = [self.sub(&minus_p, &public_p), minus_p];
            } else {
                pts[2 * k] = [
                    self.double_iter_with_ws(&p1, a_div_three, m),
                    self.double_iter_with_ws(&p2, a_div_three, m),
                ];
                pts[2 * k + 1] = [
                    self.double_iter_with_ws(&p3, a_div_three, m),
                    self.double_iter_with_ws(&p4, a_div_three, m),
                ];
            }
        }
        k
    }

    /// Compute the isogeny chain (E x E^sigma) -> A using a balanced strategy used for Mike
    /// secret generation
    pub fn mike_isogeny_chain(
        &self,
        a_div_three: &Fp2,
        gluing_kernel: &KernelPoints<Fp2>,
        chain_kernel: &[[Point<Fp2>; 2]; 2],
        len: usize,
    ) -> ThetaStructure<Fp> {
        self.mike_isogeny_chain_with_prepared_doublings(
            a_div_three,
            gluing_kernel,
            chain_kernel,
            len,
            None,
        )
    }

    pub(crate) fn mike_isogeny_chain_with_prepared_doublings(
        &self,
        a_div_three: &Fp2,
        gluing_kernel: &KernelPoints<Fp2>,
        chain_kernel: &[[Point<Fp2>; 2]; 2],
        len: usize,
        prepared: Option<&PreparedDoublings<Fp2>>,
    ) -> ThetaStructure<Fp> {
        self.start_mike_isogeny_chain(a_div_three, gluing_kernel, chain_kernel, len, prepared)
            .finish(len)
    }

    /// Builds the public gluing portion of a MIKE chain.
    ///
    /// The returned state can be finished independently or together with three compatible
    /// states.
    pub(crate) fn start_mike_isogeny_chain(
        &self,
        a_div_three: &Fp2,
        gluing_kernel: &KernelPoints<Fp2>,
        chain_kernel: &[[Point<Fp2>; 2]; 2],
        len: usize,
        prepared: Option<&PreparedDoublings<Fp2>>,
    ) -> MikeChainState<Fp> {
        let n = len - 2;
        let space = (usize::BITS - n.leading_zeros() + 1) as usize;

        let mut kernel_points: Vec<[Point<Fp2>; 2]> = vec![[Point::INFINITY; 2]; 2 * space];
        let mut orders: Vec<usize> = vec![0; space];
        kernel_points[0] = chain_kernel[0];
        kernel_points[1] = chain_kernel[1];
        orders[0] = n;

        let mut k =
            self.advance_strategy(a_div_three, &mut kernel_points, &mut orders, 0, prepared);

        // Compute the gluing isogeny and push through all the points we doubled for the strategy
        let (domain, mut theta_kernel) =
            gluing_isogeny(self, gluing_kernel, &kernel_points[..(2 * k)]);
        // Gluing isogeny creates theta_kernel of length 2*k but we need it to have a minimum size
        theta_kernel.resize(2 * space, ThetaPoint::default());

        // Reduce the order of every point pushed through the isogeny
        for ord in orders.iter_mut().take(k) {
            *ord -= 1;
        }
        k = k.saturating_sub(1);

        MikeChainState {
            domain,
            theta_kernel,
            orders,
            k,
        }
    }
}

/// A chain after gluing, ready for scalar or packed theta evaluation.
#[derive(Clone)]
pub(crate) struct MikeChainState<Fp: FqTrait> {
    pub(crate) domain: ThetaStructure<Fp>,
    pub(crate) theta_kernel: Vec<ThetaPoint<Fp>>,
    pub(crate) orders: Vec<usize>,
    pub(crate) k: usize,
}

impl<Fp: FqTrait> MikeChainState<Fp> {
    /// Finishes the chain.
    ///
    /// This uses the coordinate SIMD path when available.
    pub(crate) fn finish(self, len: usize) -> ThetaStructure<Fp> {
        let Self {
            mut domain,
            mut theta_kernel,
            mut orders,
            mut k,
        } = self;

        #[cfg(feature = "simd")]
        if let Some(codomain) =
            crate::theta::coordinate_simd::try_finish(domain, &theta_kernel, &orders, k, len)
        {
            return codomain;
        }

        for step in 0..len {
            // Perform doublings of the kernel elements, decreasing the values of orders
            k = domain.advance_strategy(&mut theta_kernel, &mut orders, k);

            // Extract out the kernel for this step.
            let T1 = theta_kernel[2 * k];
            let T2 = theta_kernel[2 * k + 1];

            // Perform one step of the (2,2) isogeny and push through all points.
            let inv_null_dual =
                ThetaStructure::two_isogeny_evaluate(&T1, &T2, &mut theta_kernel[..(2 * k)]);
            // Only the public strategy controls this decision. When the next
            // kernel is already ready, its evaluation does not use `domain`.
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
