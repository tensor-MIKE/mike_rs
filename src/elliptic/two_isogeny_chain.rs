use fp2::traits::Fp2 as Fp2Trait;

use super::{curve::Curve, x_point::PointX};
use crate::utils::ct_sort::ct_find_smallest_in_array;

impl<Fp2: Fp2Trait> Curve<Fp2> {
    /// Compute the codomain of the 4-isogeny E -> E/<ker> for [2]ker != (0 : 1)
    /// Returns the codomain (A24 : C24) together with three constants (c0, c1, c1)
    /// used for computing images.
    fn four_isogeny_codomain(ker: &PointX<Fp2>) -> (Fp2, Fp2, Fp2, Fp2, Fp2) {
        let c1 = ker.X - ker.Z;
        let c2 = ker.X + ker.Z;
        let mut c0 = ker.Z.square();
        c0.set_mul2();
        let C24 = c0.square();
        c0.set_mul2();
        let mut A24: Fp2 = ker.X.square();
        A24.set_mul2();
        A24.set_square();
        (A24, C24, c0, c1, c2)
    }

    /// Evaluate a point Q in place under the action of the 4-isogeny E -> E/<ker>
    /// for [2]ker != (0 : 1)
    fn four_isogeny_eval(c0: &Fp2, c1: &Fp2, c2: &Fp2, Q: &mut PointX<Fp2>) {
        let mut t0 = Q.X + Q.Z;
        let mut t1 = Q.X - Q.Z;
        Q.X = t0 * (*c1);
        Q.Z = t1 * (*c2);
        t0 *= t1;
        t0 *= *c0;
        t1 = Q.X + Q.Z;
        Q.Z = Q.X - Q.Z;
        t1.set_square();
        Q.Z.set_square();
        Q.X = t0 + t1;
        t0 = Q.Z - t0;
        Q.X *= t1;
        Q.Z *= t0;
    }

    /// Given a point of four torsion on an elliptic curve, compute all six isomorphic
    /// Montgomery coefficients as well as their Galois conjugates and pick the
    /// lexicographically smallest one.
    /// Ensures that the output of the two isogeny chain leaks no information about the kernel
    /// and that keys which project curves which are equal under Frobenius produce the same
    /// public key
    pub fn normalize_curve(four_torsion: &PointX<Fp2>) -> Curve<Fp2> {
        // Compute the three theta points
        // (a + b : a - b)
        // (a : b)
        // (a*i + b : a + i*b)
        // Assuming the Montgomery point is (a + b : a - b)
        let (a0, b0) = four_torsion.coords();
        let (a1, b1) = (a0 + b0, a0 - b0);
        // TODO: add method to fp2 library for efficient mul by i?
        let (a2, b2) = (a1 * Fp2::ZETA + b1, b1 * Fp2::ZETA + a1);

        // Compute ai^4 and bi^4
        let a0_4 = a0.n_square(2);
        let a1_4 = a1.n_square(2);
        let a2_4 = a2.n_square(2);

        let b0_4 = b0.n_square(2);
        let b1_4 = b1.n_square(2);
        let b2_4 = b2.n_square(2);

        // Compute the three distinct A coefficients
        let mut A1 = -(a0_4 + b0_4).mul2();
        let mut A2 = -(a1_4 + b1_4).mul2();
        let mut A3 = -(a2_4 + b2_4).mul2();
        let C1 = a0_4 - b0_4;
        let C2 = a1_4 - b1_4;
        let C3 = a2_4 - b2_4;
        let mut inv = [C1, C2, C3];
        Fp2::batch_invert(&mut inv);

        // Compute the 6 normalised coefficients
        A1 *= inv[0];
        A2 *= inv[1];
        A3 *= inv[2];
        let A4 = -A1;
        let A5 = -A2;
        let A6 = -A3;

        // Encode them to canonical bytes and sort them lexiographically
        let encoded_coefficients: [Fp2::Encoded; 12] = [
            A1.encode(),
            A2.encode(),
            A3.encode(),
            A4.encode(),
            A5.encode(),
            A6.encode(),
            A1.conjugate().encode(),
            A2.conjugate().encode(),
            A3.conjugate().encode(),
            A4.conjugate().encode(),
            A5.conjugate().encode(),
            A6.conjugate().encode(),
        ];
        let smallest_coefficient = ct_find_smallest_in_array(&encoded_coefficients);

        // TODO: we create a curve here as the PublicKey wants this, but we could just return the bytes
        // and represent the output as the encoded public key, it doesn't really matter but this computation
        // is essentially to please types rather than because we need this curve object.
        let (A, ok): (Fp2, u32) = Fp2::decode(smallest_coefficient.as_ref());
        debug_assert!(ok == u32::MAX);
        let A24 = (A + Fp2::TWO).half().half();

        Curve::new_no_checks(&A, &A24)
    }

    /// Compute a 2^2 isogeny using a balanced strategy with 4-isogenies for
    /// every step, returning the codomain and a point of four torsion on the codomain.
    /// Requires as input a point with four torsion above the kernel.
    /// i.e. kernel has order 2^(e + 2)
    pub fn two_isogeny_chain(&self, kernel: &PointX<Fp2>, e: usize) -> Self {
        debug_assert_eq!(e & 1, 0); // Ensure that the length is even

        // For 4-isogenies we represent (A + 2) / 4 projectively as (A24 : C24)
        let mut A24 = self.A24;
        let mut C24 = Fp2::ONE;

        // Precompute constants from the codomain at each step for computing images.
        let mut c0;
        let mut c1;
        let mut c2;

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - e.leading_zeros()) as usize;

        // These are a set of points of order 2^i
        let mut strategy_points: Vec<PointX<Fp2>> = vec![PointX::INFINITY; space];

        // The values i such that each point in strategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Initalise the first values for the strategy
        strategy_points[0] = *kernel;
        orders[0] = e + 2;

        let mut k = 0;
        for _ in 0..(e >> 1) {
            // Get the next point of order 4
            while orders[k] != 2 {
                k += 1;
                let m = 2 * (orders[k - 1] / 4) + (orders[k - 1] & 1);
                strategy_points[k] = strategy_points[k - 1];
                Self::xdbl_proj_iter_into(&mut strategy_points[k], &A24, &C24, m);
                orders[k] = orders[k - 1] - m;
            }
            // Compute the codomain from the current step
            (A24, C24, c0, c1, c2) = Self::four_isogeny_codomain(&strategy_points[k]);

            // Push through the kernel points and reduce the stored order
            for i in 0..k {
                Self::four_isogeny_eval(&c0, &c1, &c2, &mut strategy_points[i]);
                orders[i] = orders[i].saturating_sub(2);
            }

            k = k.saturating_sub(1);
        }
        Self::normalize_curve(&strategy_points[0])
    }
}
