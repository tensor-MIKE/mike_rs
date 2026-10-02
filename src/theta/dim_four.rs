use fp2::traits::Fq as FqTrait;
use std::ops::{Index, IndexMut};

/// Theta Point in with a projective representation for dimension 4 with
/// named fields x0 through to x15
#[derive(Clone, Copy, Debug)]
pub struct ThetaPoint<Fq: FqTrait> {
    coords: [Fq; 16],
}

impl<Fq: FqTrait> ThetaPoint<Fq> {
    pub const fn new_from_array(coords: [Fq; 16]) -> Self {
        Self { coords }
    }
    /// Checks all coordinates for the first non-zero. Returns
    /// (coordinate, 0xFFFFFFFF) if there is a non-zero coordinate
    /// and (0, 0) if no such coordinate is found.
    fn first_non_zero_coordinate(&self) -> (Fq, u32) {
        let mut result = Fq::ZERO;
        let mut found = 0u32;
        for c in self.coords {
            let take = !c.is_zero() & !found;
            result.set_cond(&c, take);
            found |= take;
        }
        (result, found)
    }

    /// Return 0xFFFFFFFF if self and rhs represent the same projective point,
    /// 0x00000000 otherwise.
    pub fn equals(&self, rhs: &Self) -> u32 {
        let (sp, sp_non_zero) = self.first_non_zero_coordinate();
        let (qp, sq_non_zero) = rhs.first_non_zero_coordinate();

        // Two points are projectively equal if the following holds
        // for all coordinates providing sp and sq are non zero
        let mut coordinates_equal = sp_non_zero & sq_non_zero;
        for i in 0..16 {
            coordinates_equal &= (self[i] * qp).equals(&(rhs[i] * sp));
        }

        // The points are equal if the coordinates are equal, or if
        // all coordinates are zero for both points.
        (!sp_non_zero & !sq_non_zero) | coordinates_equal
    }

    /// Compute the Hadamard transformation of the point in place.
    ///
    /// Applies 4 levels of Cooley-Tukey butterflies. Each butterfly
    /// (a, b) -> (a + b, a - b) costs 2a, with 8 butterflies per level.
    ///
    /// Cost: 4 * 8 * 2 = 64a
    #[inline(always)]
    pub fn set_hadamard(&mut self) {
        #[rustfmt::skip]
        let [
            x0, x1, x2, x3, x4, x5, x6, x7,
            x8, x9, x10,x11,x12,x13,x14,x15
        ] = self.coords;

        let (t0, t1) = (x0 + x1, x0 - x1);
        let (t2, t3) = (x2 + x3, x2 - x3);
        let (t4, t5) = (x4 + x5, x4 - x5);
        let (t6, t7) = (x6 + x7, x6 - x7);
        let (t8, t9) = (x8 + x9, x8 - x9);
        let (t10, t11) = (x10 + x11, x10 - x11);
        let (t12, t13) = (x12 + x13, x12 - x13);
        let (t14, t15) = (x14 + x15, x14 - x15);

        let (u0, u2) = (t0 + t2, t0 - t2);
        let (u1, u3) = (t1 + t3, t1 - t3);
        let (u4, u6) = (t4 + t6, t4 - t6);
        let (u5, u7) = (t5 + t7, t5 - t7);
        let (u8, u10) = (t8 + t10, t8 - t10);
        let (u9, u11) = (t9 + t11, t9 - t11);
        let (u12, u14) = (t12 + t14, t12 - t14);
        let (u13, u15) = (t13 + t15, t13 - t15);

        let (v0, v4) = (u0 + u4, u0 - u4);
        let (v1, v5) = (u1 + u5, u1 - u5);
        let (v2, v6) = (u2 + u6, u2 - u6);
        let (v3, v7) = (u3 + u7, u3 - u7);
        let (v8, v12) = (u8 + u12, u8 - u12);
        let (v9, v13) = (u9 + u13, u9 - u13);
        let (v10, v14) = (u10 + u14, u10 - u14);
        let (v11, v15) = (u11 + u15, u11 - u15);

        self.coords = [
            v0 + v8,
            v1 + v9,
            v2 + v10,
            v3 + v11,
            v4 + v12,
            v5 + v13,
            v6 + v14,
            v7 + v15,
            v0 - v8,
            v1 - v9,
            v2 - v10,
            v3 - v11,
            v4 - v12,
            v5 - v13,
            v6 - v14,
            v7 - v15,
        ];
    }

    pub fn set_apply_gluing_basis_change(&mut self) {
        let t1 = self.coords[0] + self.coords[5];
        let t2 = self.coords[0] - self.coords[5];
        let t3 = self.coords[10] + self.coords[15];
        let t4 = self.coords[10] - self.coords[15];

        let t5 = self.coords[4] + self.coords[1];
        let t6 = self.coords[4] - self.coords[1];
        let t7 = self.coords[14] + self.coords[11];
        let t8 = self.coords[14] - self.coords[11];

        let t9 = self.coords[2] + self.coords[7];
        let t10 = self.coords[2] - self.coords[7];
        let t11 = self.coords[8] + self.coords[13];
        let t12 = self.coords[8] - self.coords[13];

        let t13 = self.coords[6] + self.coords[3];
        let t14 = self.coords[6] - self.coords[3];
        let t15 = self.coords[12] + self.coords[9];
        let t16 = self.coords[12] - self.coords[9];

        self.coords[0] = t1 + t3;
        self.coords[1] = t2 + t4;
        self.coords[2] = t1 - t3;
        self.coords[3] = t2 - t4;

        self.coords[4] = t5 + t7;
        self.coords[5] = t6 + t8;
        self.coords[6] = t5 - t7;
        self.coords[7] = t6 - t8;

        self.coords[8] = t9 + t11;
        self.coords[9] = t10 + t12;
        self.coords[10] = t9 - t11;
        self.coords[11] = t10 - t12;

        self.coords[12] = t13 + t15;
        self.coords[13] = t14 + t16;
        self.coords[14] = t13 - t15;
        self.coords[15] = t14 - t16;
    }

    /// Return the Hadamard transform of this point.
    ///
    /// Cost: 64a
    pub fn hadamard(&self) -> Self {
        let mut p = *self;
        p.set_hadamard();
        p
    }

    /// Compute the squaring transform of the point in place.
    ///
    /// Cost: 16s
    #[inline(always)]
    pub fn set_square(&mut self) {
        for c in &mut self.coords {
            c.set_square()
        }
    }

    /// Return the squaring transform of this point.
    ///
    /// Cost: 16s
    pub fn square(&self) -> Self {
        let mut p = *self;
        p.set_square();
        p
    }

    /// Compute the coordinate-wise multiplication of the point in place.
    ///
    /// Cost: 16m
    #[inline(always)]
    pub fn set_coordinate_multiply(&mut self, other: &Self) {
        for i in 0..16 {
            self.coords[i] *= other.coords[i]
        }
    }

    /// Return the coordinate-wise multiplication of two points
    ///
    /// Cost: 16m
    pub fn coordinate_multiply(&self, other: &Self) -> Self {
        let mut p = *self;
        p.set_coordinate_multiply(other);
        p
    }
}

/// Default element used for initialisation
impl<Fq: FqTrait> Default for ThetaPoint<Fq> {
    fn default() -> Self {
        Self {
            coords: [Fq::ZERO; 16],
        }
    }
}

impl<Fq: FqTrait> Index<usize> for ThetaPoint<Fq> {
    type Output = Fq;

    fn index(&self, i: usize) -> &Fq {
        &self.coords[i]
    }
}

impl<Fq: FqTrait> IndexMut<usize> for ThetaPoint<Fq> {
    fn index_mut(&mut self, i: usize) -> &mut Fq {
        &mut self.coords[i]
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for ThetaPoint<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        writeln!(f, "ThetaPoint: (")?;
        for i in 0..16 {
            writeln!(f, "    {},", self[i])?;
        }
        write!(f, ")")
    }
}

/// A compressed Theta Point using only the 10 independent coordinates
/// of a dimension 4 theta null point.
///
/// The mapping from full 16 coordinates to compressed 10 coordinates is:
/// [x0, x1, x2, x3, x5, x6, x7, x9, x11, x15]
///
/// The inverse mapping (comp_to_theta) is:
/// [z0, z1, z2, z3, z2, z4, z5, z6, z1, z7, z4, z8, z3, z8, z6, z9]
#[derive(Clone, Copy, Debug)]
pub struct ThetaPointCompressed<Fq: FqTrait> {
    coords: [Fq; 10],
}

impl<Fq: FqTrait> Index<usize> for ThetaPointCompressed<Fq> {
    type Output = Fq;

    fn index(&self, i: usize) -> &Fq {
        &self.coords[i]
    }
}

impl<Fq: FqTrait> IndexMut<usize> for ThetaPointCompressed<Fq> {
    fn index_mut(&mut self, i: usize) -> &mut Fq {
        &mut self.coords[i]
    }
}

impl<Fq: FqTrait> ThetaPointCompressed<Fq> {
    /// Compress a full ThetaPoint into a ThetaPointCompressed.
    /// [x0, x1, x2, x3, x5, x6, x7, x9, x11, x15]
    ///
    pub fn from_theta_point(p: &ThetaPoint<Fq>) -> Self {
        Self {
            coords: [p[0], p[1], p[2], p[3], p[5], p[6], p[7], p[9], p[11], p[15]],
        }
    }

    /// Expand a ThetaPointCompressed into a full ThetaPoint.
    /// [z0, z1, z2, z3, z2, z4, z5, z6, z1, z7, z4, z8, z3, z8, z6, z9]
    ///
    pub fn to_theta_point(self) -> ThetaPoint<Fq> {
        ThetaPoint {
            coords: [
                self[0], self[1], self[2], self[3], self[2], self[4], self[5], self[6], self[1],
                self[7], self[4], self[8], self[3], self[8], self[6], self[9],
            ],
        }
    }

    /// Checks all coordinates for the first non-zero. Returns
    /// (coordinate, 0xFFFFFFFF) if there is a non-zero coordinate
    /// and (0, 0) if no such coordinate is found.
    fn first_non_zero_coordinate(&self) -> (Fq, u32) {
        let mut result = Fq::ZERO;
        let mut found = 0u32;
        for coord in self.coords {
            let take = !coord.is_zero() & !found;
            result.set_cond(&coord, take);
            found |= take;
        }
        (result, found)
    }

    /// Return 0xFFFFFFFF if self and rhs represent the same projective point,
    /// 0x00000000 otherwise.
    pub fn equals(&self, rhs: &Self) -> u32 {
        let (sp, sp_non_zero) = self.first_non_zero_coordinate();
        let (qp, sq_non_zero) = rhs.first_non_zero_coordinate();

        // Two points are projectively equal if the following holds
        // for all coordinates providing sp and sq are non zero
        let mut coordinates_equal = sp_non_zero & sq_non_zero;
        for i in 0..10 {
            coordinates_equal &= (self[i] * qp).equals(&(rhs[i] * sp));
        }

        // The points are equal if the coordinates are equal, or if
        // all coordinates are zero for both points.
        (!sp_non_zero & !sq_non_zero) | coordinates_equal
    }

    /// Compute the squaring transform in place.
    ///
    /// Cost: 10s
    #[inline(always)]
    pub fn set_square(&mut self) {
        for c in &mut self.coords {
            c.set_square();
        }
    }

    /// Return the squaring transform of this point.
    ///
    /// Cost: 10s
    pub fn square(&self) -> Self {
        let mut p = *self;
        p.set_square();
        p
    }

    /// Compute the projective pseudo-inverse in place.
    ///
    /// Cost: 3(n-1)M = 27M
    #[inline(always)]
    pub fn set_proj_batch_pseudo_inversion(&mut self) {
        let mut multiples = [Fq::ONE; 10];
        multiples[0] = self[0];
        for i in 1..10 {
            multiples[i] = multiples[i - 1] * self[i];
        }

        let mut inverses = [Fq::ONE; 10];
        for i in 1..10 {
            inverses[i] = inverses[i - 1] * self[10 - i];
        }

        self[0] = inverses[9];
        for i in 1..10 {
            self[i] = inverses[9 - i] * multiples[i - 1];
        }
    }

    /// Return the projective pseudo-inverse of this point.
    ///
    /// Cost: 27M
    pub fn proj_batch_pseudo_inversion(&self) -> Self {
        let mut p = *self;
        p.set_proj_batch_pseudo_inversion();
        p
    }

    /// Compute the Hadamard transform of a compressed theta null point in place.
    ///
    /// This implements the optimised transform due to Sabrina, exploiting the
    /// symmetry structure of the compressed representation to reduce the cost
    /// from 64a (full point) to 34a.
    ///
    /// Cost: 10a + 6a + 6a + 2a + 10a = 34a
    #[inline(always)]
    pub fn set_hadamard(&mut self) {
        // Unpack the coordinates
        let [z0, z1, z2, z3, z4, z5, z6, z7, z8, z9] = self.coords;

        // 10a
        let (t0, t1) = (z0 + z9, z0 - z9);
        let (t2, t3) = (z5 + z7, z5 - z7);
        let (t4, t5) = (z2 + z8, z2 - z8);
        let (t6, t7) = (z1 + z6, z1 - z6);
        let (t8, t9) = (z3 + z4, z3 - z4);

        // 6a
        let (s0, s1) = (t0 + t2, t0 - t2);
        let (s2, s3) = (t1 + t3, t1 - t3);
        let (s4, s5) = (t6 + t4, t6 - t4);

        // 6 double
        let r0 = t5.mul2();
        let r1 = t7.mul2();
        let r2 = t9.mul2();
        let r3 = t8.mul2();
        let r4 = s4.mul2();
        let r5 = s5.mul2();

        // 2a
        let (u0, u1) = (s0 + r3, s0 - r3);

        // 10a
        self.coords = [
            u0 + r4,
            s2 + r0,
            s3 + r1,
            s1 + r2,
            s1 - r2,
            u1 + r5,
            s2 - r0,
            u1 - r5,
            s3 - r1,
            u0 - r4,
        ];
    }

    /// Return the Hadamard transform of this compressed point.
    ///
    /// Cost: 34a
    pub fn hadamard(&self) -> Self {
        let mut p = *self;
        p.set_hadamard();
        p
    }
}

impl<Fq: FqTrait> Default for ThetaPointCompressed<Fq> {
    fn default() -> Self {
        Self {
            coords: [Fq::ZERO; 10],
        }
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for ThetaPointCompressed<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        writeln!(f, "ThetaPointCompressed: (")?;
        for i in 0..10 {
            writeln!(f, "    {},", self[i])?;
        }
        write!(f, ")")
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThetaStructure<Fq: FqTrait> {
    null_point: ThetaPoint<Fq>,
    inv_null_point: ThetaPoint<Fq>,
    inv_null_point_dual_sq: ThetaPoint<Fq>,
}

impl<Fq: FqTrait> ThetaStructure<Fq> {
    pub fn new_from_null_point(null_point: ThetaPoint<Fq>) -> Self {
        // Compress and do all arithmetic in the cheaper 10-coordinate representation
        let null_point_comp = ThetaPointCompressed::from_theta_point(&null_point);

        let inv_null_point_comp = null_point_comp.proj_batch_pseudo_inversion();
        let mut inv_null_point_dual_sq_comp = null_point_comp.square();
        inv_null_point_dual_sq_comp.set_hadamard();
        inv_null_point_dual_sq_comp.set_proj_batch_pseudo_inversion();

        // Decompress all three points at the end
        Self {
            null_point,
            inv_null_point: inv_null_point_comp.to_theta_point(),
            inv_null_point_dual_sq: inv_null_point_dual_sq_comp.to_theta_point(),
        }
    }

    pub fn new_from_inv_null_dual(inv_null_point_dual: ThetaPoint<Fq>) -> Self {
        // Compress and do all arithmetic in the cheaper 10-coordinate representation
        let inv_null_point_dual_comp = ThetaPointCompressed::from_theta_point(&inv_null_point_dual);

        let null_point_dual_comp = inv_null_point_dual_comp.proj_batch_pseudo_inversion();
        let null_point_comp = null_point_dual_comp.hadamard();

        // Can we speed this up using the inverse we start with?
        let mut inv_null_point_dual_sq_comp = null_point_comp.square();
        inv_null_point_dual_sq_comp.set_hadamard();
        inv_null_point_dual_sq_comp.set_proj_batch_pseudo_inversion();

        // Can we speed this up using the inverse we start with?
        let inv_null_point_comp = null_point_comp.proj_batch_pseudo_inversion();

        // Decompress all three points at the end
        Self {
            null_point: null_point_comp.to_theta_point(),
            inv_null_point: inv_null_point_comp.to_theta_point(),
            inv_null_point_dual_sq: inv_null_point_dual_sq_comp.to_theta_point(),
        }
    }

    /// Return the null point of the theta structure
    pub fn null_point(&self) -> ThetaPoint<Fq> {
        self.null_point
    }

    /// Returns 0xFFFFFFFF if a theta point is zero, and 0 otherwise
    pub fn is_zero(&self, other: &ThetaPoint<Fq>) -> u32 {
        self.null_point.equals(other)
    }

    /// Compute the double of a point in place
    /// Cost 32S + 32M + 128a
    pub fn double_into(&self, p: &mut ThetaPoint<Fq>) {
        p.set_square();
        p.set_hadamard();
        p.set_square();
        p.set_coordinate_multiply(&self.inv_null_point_dual_sq);
        p.set_hadamard();
        p.set_coordinate_multiply(&self.inv_null_point);
    }

    /// Return the double of a given point
    /// Cost 32S + 32M + 128a
    pub fn double(&self, p: &ThetaPoint<Fq>) -> ThetaPoint<Fq> {
        let mut q = *p;
        self.double_into(&mut q);
        q
    }

    /// P <- [2^n]*P in place.
    /// Cost: n * (32S + 32M + 128a).
    pub fn double_iter_into(&self, p: &mut ThetaPoint<Fq>, n: usize) {
        for _ in 0..n {
            self.double_into(p);
        }
    }

    /// Return [2^n]*P
    /// Cost: n * (32S + 32M + 128a).
    pub fn double_iter(&self, p: &ThetaPoint<Fq>, n: usize) -> ThetaPoint<Fq> {
        let mut q = *p;
        self.double_iter_into(&mut q, n);
        q
    }
}

#[cfg(test)]
mod tests {
    use super::{ThetaPoint, ThetaPointCompressed, ThetaStructure};
    use crate::fields::MikeTestFp as Fp;
    use crate::theta::test_data::{OA_NULL, P, Q};
    use crate::utils::test_utils::drng::DRNG;
    use rand_core::{CryptoRng, RngCore};

    /// Build a random ThetaPoint with independently random coordinates.
    fn rand_theta_point<R: CryptoRng + RngCore>(rng: &mut R) -> ThetaPoint<Fp> {
        ThetaPoint {
            coords: [Fp::rand(rng); 16],
        }
    }

    fn rand_theta_null_point<R: CryptoRng + RngCore>(rng: &mut R) -> ThetaPoint<Fp> {
        let comp = ThetaPointCompressed {
            coords: [Fp::rand(rng); 10],
        };
        comp.to_theta_point()
    }

    /// H(H(P)) should be projectively equal to P (up to the scalar 2^4 = 16).
    #[test]
    fn test_hadamard_involution() {
        let mut rng = DRNG::from_seed("hadamard_twice".as_bytes());
        let p = rand_theta_point(&mut rng);
        let hh = p.hadamard().hadamard();
        assert_eq!(p.equals(&hh), u32::MAX);
    }

    /// The compressed Hadamard should agree with the full Hadamard after
    /// round-tripping through compress/decompress.
    #[test]
    fn test_hadamard_compressed_matches_full() {
        let mut rng = DRNG::from_seed("hadamard_compressed".as_bytes());
        let p = rand_theta_null_point(&mut rng);

        // Full Hadamard
        let h_full = p.hadamard();

        // Compressed path: compress → hadamard → decompress
        let h_comp = ThetaPointCompressed::from_theta_point(&p)
            .hadamard()
            .to_theta_point();

        assert_eq!(h_full.equals(&h_comp), u32::MAX);
    }

    /// The compressed square should agree with the full square after
    /// round-tripping through compress/decompress.
    #[test]
    fn test_square_compressed_matches_full() {
        let mut rng = DRNG::from_seed("squaring".as_bytes());
        let p = rand_theta_null_point(&mut rng);

        // Full square
        let sq_full = p.square();

        // Compressed path: compress → square → decompress
        let sq_comp = ThetaPointCompressed::from_theta_point(&p)
            .square()
            .to_theta_point();

        assert_eq!(sq_full.equals(&sq_comp), u32::MAX);
    }

    /// P.pseudo_inv().pseudo_inv() should be projectively equal to P.
    #[test]
    fn test_pseudo_inversion_involution_compressed() {
        let mut rng = DRNG::from_seed("inverse_compressed".as_bytes());
        let p = rand_theta_point(&mut rng);
        let p_comp = ThetaPointCompressed::from_theta_point(&p);
        let pp = p_comp
            .proj_batch_pseudo_inversion()
            .proj_batch_pseudo_inversion();
        assert_eq!(p_comp.equals(&pp), u32::MAX);
    }

    #[test]
    fn test_doubling() {
        let oa = ThetaStructure::new_from_null_point(OA_NULL);

        // These should be points of order 2
        let p_test = oa.double_iter(&P, 244 - 1);
        let q_test = oa.double_iter(&Q, 244 - 1);

        // So doubling them should give zero
        let p_zero = oa.double(&p_test);
        let q_zero = oa.double(&q_test);

        assert_eq!(oa.is_zero(&p_test), 0u32);
        assert_eq!(oa.is_zero(&q_test), 0u32);
        assert_eq!(oa.is_zero(&p_zero), u32::MAX);
        assert_eq!(oa.is_zero(&q_zero), u32::MAX);
    }
}
