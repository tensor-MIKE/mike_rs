use crate::theta::dim_four::{ThetaPoint, ThetaStructure};
use fp2::traits::Fq as FqTrait;

impl<Fq: FqTrait> ThetaStructure<Fq> {
    #[inline(always)]
    fn compute_inv_dual_codomain(P8: &ThetaPoint<Fq>, Q8: &ThetaPoint<Fq>) -> ThetaPoint<Fq> {
        // Precompute the prefix chain
        // 7M
        let r1 = P8[12];
        let r2 = r1 * P8[5];
        let r3 = r2 * Q8[1];
        let r4 = r3 * P8[9];
        let r5 = r4 * P8[11];
        let r6 = r5 * Q8[15];
        let r7 = r6 * Q8[14];
        let r8 = r7 * P8[6];

        // Precompute the suffix chain
        // 7M
        let b8 = P8[0];
        let b7 = b8 * P8[2];
        let b6 = b7 * Q8[6];
        let b5 = b6 * Q8[7];
        let b4 = b5 * P8[15];
        let b3 = b4 * P8[13];
        let b2 = b3 * Q8[9];
        let b1 = b2 * P8[1];

        let mut inv_dual_theta = ThetaPoint::default();

        // First element
        // 1M
        inv_dual_theta[0] = r8 * P8[4];

        // Last element required an extra mul
        // 2M
        let tmp = P8[8] * P8[5];
        inv_dual_theta[12] = tmp * b2;

        // Remaining coordinates
        // 8M
        inv_dual_theta[5] = r1 * b1;
        inv_dual_theta[1] = r2 * b2;
        inv_dual_theta[9] = r3 * b3;
        inv_dual_theta[11] = r4 * b4;
        inv_dual_theta[15] = r5 * b5;
        inv_dual_theta[14] = r6 * b6;
        inv_dual_theta[6] = r7 * b7;
        inv_dual_theta[4] = r8 * b8;

        // Direct copies from symmetry
        inv_dual_theta[3] = inv_dual_theta[12];
        inv_dual_theta[10] = inv_dual_theta[5];
        inv_dual_theta[8] = inv_dual_theta[1];
        inv_dual_theta[13] = inv_dual_theta[11];
        inv_dual_theta[7] = inv_dual_theta[14];
        inv_dual_theta[2] = inv_dual_theta[4];

        inv_dual_theta
    }

    pub fn two_isogeny(
        T1: &ThetaPoint<Fq>,
        T2: &ThetaPoint<Fq>,
        image_points: &mut [ThetaPoint<Fq>],
    ) -> Self {
        let mut P8 = T1.square();
        let mut Q8 = T2.square();
        P8.set_hadamard();
        Q8.set_hadamard();

        let inv_null_point_dual = Self::compute_inv_dual_codomain(&P8, &Q8);

        for P in image_points.iter_mut() {
            P.set_square();
            P.set_hadamard();
            P.set_coordinate_multiply(&inv_null_point_dual);
            P.set_hadamard();
        }

        ThetaStructure::new_from_inv_null_dual(inv_null_point_dual)
    }
}
