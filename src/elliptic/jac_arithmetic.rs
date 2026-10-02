use crate::elliptic::x_point::PointX;

use super::{curve::Curve, jac_point::Point};
use fp2::traits::{Fq as FqTrait, FqRnd};
use rand_core::{CryptoRng, RngCore};

/// Projective representation of a modified Jacobian point (X : Y : Z : T)
/// used for fast repeated doubling
#[derive(Clone, Copy, Debug)]
pub struct PointW<Fq: FqTrait> {
    pub X: Fq,
    pub Y: Fq,
    pub Z: Fq,
    pub T: Fq,
}

impl<Fq: FqTrait> Curve<Fq> {
    /// Compute R = [2]*P in Jacobian coordinates, writing the result into P.
    /// Cost: 6S + 6M + 14a
    #[inline(always)]
    pub fn double_into(&self, P: &mut Point<Fq>) {
        let (XP, YP, ZP) = (P.X, P.Y, P.Z);

        let t0 = XP.square().mul3();
        let t1 = ZP.square();
        let mut t2 = (XP * self.A).mul2();
        t2 += t1;
        t2 *= t1;
        t2 += t0;

        let z2p = (YP * ZP).mul2();

        let mut x2p = t2.square();
        let t0: Fq = z2p.square() * self.A;
        let t1 = YP.square().mul2();
        let t3 = t1 * XP.mul2();
        x2p -= t0;
        x2p -= t3.mul2();

        let mut y2p = t3 - x2p;
        y2p *= t2;
        y2p -= t1.square().mul2();

        P.X = x2p;
        P.Y = y2p;
        P.Z = z2p;
    }

    /// Convert a point from Jacobian to modified Weierstrass
    /// NOTE: requires the precomputation of A_div_three
    fn convert_to_ws(&self, A_div_three: &Fq, P: &Point<Fq>) -> PointW<Fq> {
        // a = 1 - A^2 / 3
        let a = Fq::ONE - self.A * *A_div_three;

        // X = X + Z^2 * A/3
        // Y, Z = Y, Z
        // T = a * Z^4
        let z2 = P.Z.square();
        let z4 = z2.square();
        let x = P.X + z2 * *A_div_three;
        let t = a * z4;

        PointW {
            X: x,
            Y: P.Y,
            Z: P.Z,
            T: t,
        }
    }

    /// Convert a point from modified Weierstrass to Jacobian
    /// NOTE: requires the precomputation of A_div_three
    fn convert_from_ws(A_div_three: &Fq, P: &PointW<Fq>) -> Point<Fq> {
        // X = X - Z^2 * A / 3
        // Y, Z = Y, Z
        let x = P.X - P.Z.square() * *A_div_three;
        Point {
            X: x,
            Y: P.Y,
            Z: P.Z,
        }
    }

    /// Cost: 3M + 5S
    pub fn double_ws_into(&self, P: &mut PointW<Fq>) {
        let xx = P.X.square();
        let c = P.Y.square().mul2();
        let cc = c.square();
        let r = cc.mul2();
        let s = (P.X + c).square() - xx - cc;
        let m = xx.mul3() + P.T;
        let x3 = m.square() - s.mul2();
        let z3 = (P.Y * P.Z).mul2();
        let y3 = m * (s - x3) - r;
        let t3 = (r * P.T).mul2();

        P.X = x3;
        P.Y = y3;
        P.Z = z3;
        P.T = t3;
    }

    /// Return R = [2]*P in Jacobian coordinates.
    /// Cost: 6S + 6M + 14a
    pub fn double(&self, P: &Point<Fq>) -> Point<Fq> {
        let mut Q = *P;
        self.double_into(&mut Q);
        Q
    }

    /// Compute [2^n]*P in place. Cost n*(6M + 6S)
    pub fn double_iter_into(&self, P: &mut Point<Fq>, n: usize) {
        for _ in 0..n {
            self.double_into(P);
        }
    }

    /// Return [2^n]*P as a new point. Cost n*(6M + 6S)
    pub fn double_iter(&self, P: &Point<Fq>, n: usize) -> Point<Fq> {
        let mut P3 = *P;
        for _ in 0..n {
            self.double_into(&mut P3);
        }
        P3
    }

    /// Return [2^n]*P as a new point. Cost n*(3M + 5S) + 3S + 4M (conversion)
    pub fn double_iter_with_ws(&self, P: &Point<Fq>, A_div_three: &Fq, n: usize) -> Point<Fq> {
        // Convert to a WS point
        let mut P_WS = self.convert_to_ws(A_div_three, P);

        // Do doubling with faster formula
        for _ in 0..n {
            self.double_ws_into(&mut P_WS);
        }

        // Convert back at the end
        Self::convert_from_ws(A_div_three, &P_WS)
    }

    /// Compute R = P + Q with complete formulas, writing the result into P.
    /// Handles the cases P = 0, Q = 0, and P = Q correctly.
    /// Cost: 6S + 17M + 13a
    pub fn add_into(&self, R: &mut Point<Fq>, P: &Point<Fq>, Q: &Point<Fq>) {
        let (XP, YP, ZP) = (P.X, P.Y, P.Z);
        let (XQ, YQ, ZQ) = (Q.X, Q.Y, Q.Z);

        // Identity checks, used for constant time setting at the end
        let ctl_1 = ZP.is_zero();
        let ctl_2 = ZQ.is_zero();

        // Slope dy / dx for P ≠ Q
        let t0 = ZP.square();
        let t1 = ZQ.square();
        let v1 = t1 * ZQ;
        let t2 = t0 * ZP;
        let v1 = v1 * YP;
        let t2 = t2 * YQ;
        let mut dy = t2 - v1;
        let u2 = t0 * XQ;
        let u1 = t1 * XP;
        let mut dx = u2 - u1;

        // Slope f'(x) for P = Q
        let t1 = YP + YP;
        let mut t2 = (self.A * XP).mul2();
        t2 += t0;
        t2 *= t0;
        let t0_sq = XP.square();
        t2 += t0_sq;
        t2 += t0_sq;
        t2 += t0_sq;
        let t2 = t2 * ZQ;

        // If dx = dy = 0 then P = Q; use the tangent slope instead
        let ctl = dx.is_zero() & dy.is_zero();
        dx.set_cond(&t1, ctl);
        dy.set_cond(&t2, ctl);

        // New coordinates
        let t0 = ZP * ZQ;
        let t1 = t0.square();
        let t2 = dx.square();
        let t3 = dy.square();

        let mut XR = self.A * t1;
        XR += u1;
        XR += u2;
        XR *= t2;
        XR = t3 - XR;

        let mut YR = u1 * t2;
        YR -= XR;
        YR *= dy;
        let t3 = dx * t2;
        let t3 = t3 * v1;
        YR -= t3;

        let ZR = dx * t0;

        // update the coordinates
        R.X = XR;
        R.Y = YR;
        R.Z = ZR;
        R.set_cond(Q, ctl_1);
        R.set_cond(P, ctl_2);
    }

    /// Given distinct points P, Q and the Montgomery coefficient A, compute
    /// (u, v, w) such that:
    ///   x(P + Q) = (u - v) / w
    ///   x(P - Q) = (u + v) / w
    pub fn add_components(&self, P: &Point<Fq>, Q: &Point<Fq>) -> (Fq, Fq, Fq) {
        let (XP, YP, ZP) = (P.X, P.Y, P.Z);
        let (XQ, YQ, ZQ) = (Q.X, Q.Y, Q.Z);

        let mut t0 = ZP.square();
        let mut t1 = ZQ.square();
        let t2 = t1 * XP;
        let t3 = t0 * XQ;
        let mut t4 = YP * ZQ;
        t4 *= t1;
        let mut t5 = ZP * YQ;
        t5 *= t0;
        t0 *= t1;
        let mut t6 = t4 * t5;
        let v = t6.mul2();
        t4.set_square();
        t5.set_square();
        t4 += t5;
        t5 = t2 + t3;
        t6 = t3 + t3;
        t6 = t5 - t6;
        t6.set_square();
        t1 = self.A * t0;
        t1 += t5;
        t1 *= t6;
        let u = t4 - t1;
        let w = t6 * t0;

        (u, v, w)
    }

    /// P1 <- P1 + P2
    pub fn add_to(&self, P1: &mut Point<Fq>, P2: &Point<Fq>) {
        let mut P3 = Point::INFINITY;
        self.add_into(&mut P3, P1, P2);
        *P1 = P3;
    }

    /// Return P1 + P2 as a new point
    pub fn add(&self, P1: &Point<Fq>, P2: &Point<Fq>) -> Point<Fq> {
        let mut P3 = Point::INFINITY;
        self.add_into(&mut P3, P1, P2);
        P3
    }

    /// P3 <- P1 - P2
    pub fn sub_into(&self, P3: &mut Point<Fq>, P1: &Point<Fq>, P2: &Point<Fq>) {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.add_into(P3, P1, &nP2);
    }

    /// P1 <- P1 - P2
    pub fn sub_from(&self, P1: &mut Point<Fq>, P2: &Point<Fq>) {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.add_to(P1, &nP2);
    }

    /// Return P1 - P2 as a new point
    pub fn sub(&self, P1: &Point<Fq>, P2: &Point<Fq>) -> Point<Fq> {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.add(P1, &nP2)
    }

    /// P3 <- n*P
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored. Bits consumed
    /// start from n_start, when this is non-zero, it allows low bits to
    /// be truncated.
    pub fn mul_into(
        &self,
        P3: &mut Point<Fq>,
        P: &Point<Fq>,
        n: &[u8],
        n_start: usize,
        n_bitlen: usize,
    ) {
        // Montgomery ladder: see https://eprint.iacr.org/2017/212

        // We will need the complete 2*P at the end, to handle some
        // special cases of the formulas.
        let xP = P.to_point_x();
        let dP = self.double(P);
        let mut X0 = PointX::INFINITY;
        let mut X1 = xP;
        let mut cc = 0u32;
        if n_bitlen > 21 {
            // If n is large enough then it is worthwhile to
            // normalize the source point to affine.
            // We do not care if P = inf, since that is handled at
            // the end in the corrective steps.
            let Xp = xP.X / xP.Z;
            for i in (n_start..n_bitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_aff_into(&mut X0, &mut X1, &Xp);
                cc = ctl;
            }
        } else {
            for i in (n_start..n_bitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_into(&mut X0, &mut X1, &xP);
                cc = ctl;
            }
        }
        PointX::cond_swap(&mut X0, &mut X1, cc);

        // The Montgomery ladder is in projective, coordinates not Jacobian.
        // This means we need to change out P values in the Okeya and Sakurai
        // recovery
        let PX = P.X * P.Z;
        let PY = P.Y;
        let PZ: Fq = P.Z * xP.Z;

        // Special cases:
        //  - ladder fails if P = (0,0) (a point of order 2)
        //  - y is not reconstructed correctly if P has order 2,
        //    or if (n+1)*P = P, -P or infinity.
        let z0z = X0.Z.is_zero();
        let z1z = X1.Z.is_zero();
        let x1ex = (X1.X * PZ).equals(&(PX * X1.Z));

        // (X0/Z0) is the X coordinate of P0 = n*P
        // (X1/Z1) is the X coordinate of P1 = (n + 1)*P
        // We recompute the Y coordinate of n*P (formulas from
        // Okeya and Sakurai).
        let xxzz = (PX * X0.X) + (PZ * X0.Z);
        let xpz0 = PX * X0.Z;
        let x0zp = X0.X * PZ;
        let zz = PZ * X0.Z;
        let zzdA = self.A.mul2() * zz;
        let u = (xxzz * (xpz0 + x0zp + zzdA)) - (zzdA * zz);
        let v = PY.mul2() * zz * X1.Z;
        P3.X = X0.X * v;
        P3.Y = (u * X1.Z) - ((xpz0 - x0zp).square() * X1.X);
        P3.Z = X0.Z * v;

        // The above finds the projective coordinates (X : Y : Z),
        // with weighting (lam * X : lam * Y : lam * Z) but now we
        // want to convert to the Jacobian  coordinates with weighting
        // (lam^2 * X : lam^3 * Y : lam * Z)
        P3.X *= P3.Z;
        P3.Y *= P3.Z.square();

        // Fix result for the special cases.
        //  P = inf                          -> inf
        //  P != inf, 2*P = inf              -> inf or P (depending on n_0)
        //  2*P != inf, P0 = inf             -> inf
        //  2*P != inf, P0 != inf, P1 = inf  -> -P
        //  2*P != inf, P0 != inf, P1 = -P   -> -2*P
        let order1 = PZ.is_zero();
        let order2 = !order1 & PY.is_zero();
        let z0inf = !order1 & !order2 & z0z;
        let z1inf = !order1 & !order2 & !z0z & z1z;
        let p1mp = !order1 & !order2 & !z0z & !z1z & x1ex;

        let n_odd = ((n[0] as u32) & 1).wrapping_neg();
        P3.Z.set_cond(&Fq::ZERO, order1 | (order2 & !n_odd) | z0inf);
        P3.set_cond(P, z1inf | (order2 & n_odd));
        P3.set_cond(&dP, p1mp);
        P3.set_cond_neg(z1inf | p1mp);
    }

    /// Return n*P as a new point.
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored. Bits consumed
    /// start from n_start, when this is non-zero, it allows low bits to
    /// be truncated.
    pub fn mul(&self, P: &Point<Fq>, n: &[u8], n_start: usize, n_bitlen: usize) -> Point<Fq> {
        let mut P3 = Point::INFINITY;
        self.mul_into(&mut P3, P, n, n_start, n_bitlen);
        P3
    }
}

impl<Fq: FqTrait + FqRnd> Curve<Fq> {
    /// Set P to a random curve point.
    pub fn set_rand_point<R: CryptoRng + RngCore>(&self, rng: &mut R, P: &mut Point<Fq>) {
        // This function cannot actually return the point-at-infinity;
        // this is not a problem as long as the curve order is larger
        // than 2^128.
        P.Z = Fq::ONE;
        loop {
            P.X.set_rand(rng);
            P.Y = (((P.X + self.A) * P.X) * P.X) + P.X;
            if P.Y.legendre() >= 0 {
                P.Y.set_sqrt();

                // Randomly chooses the square root to use.
                let mut tmp = [0u8; 1];
                rng.fill_bytes(&mut tmp);
                let ctl = 0u32.wrapping_sub((tmp[0] as u32) & 1);
                P.Y.set_cond_neg(ctl);
                return;
            }
        }
    }

    /// Return a new random curve point.
    pub fn rand_point<R: CryptoRng + RngCore>(&self, rng: &mut R) -> Point<Fq> {
        let mut P = Point::INFINITY;
        self.set_rand_point(rng, &mut P);
        P
    }
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use super::*;

    use crate::elliptic::test_data::A1;
    use crate::fields::MikeTestFp2 as Fp2;
    use crate::utils::test_utils::drng::DRNG;

    const TEST_ITER: usize = 100;

    fn test_curve() -> Curve<Fp2> {
        let (curve, check) = Curve::new(&A1);
        assert_eq!(check, u32::MAX);
        curve
    }

    #[test]
    fn test_equals() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_equals".as_bytes());

        for _ in 0..TEST_ITER {
            let P = curve.rand_point(&mut rng);
            let Q = P.neg();

            // P != Q unless [2]P = O
            assert_eq!(P.equals(&P), u32::MAX);
            if Q.Y.is_zero() == u32::MAX {
                assert_eq!(P.equals(&Q), u32::MAX);
            } else {
                assert_eq!(P.equals(&Q), 0);
            }
        }
    }

    // Runs N random trials of the associativity test:
    //   (P + Q) + R  ==  P + (R + Q)
    // using three independently-sampled, mutually-distinct points.
    #[test]
    fn test_add_with_inf() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_add_with_inf".as_bytes());

        for _ in 0..TEST_ITER {
            let P = curve.rand_point(&mut rng);
            let O = Point::INFINITY;
            let mP = P.neg();

            // P - P = O
            let PmP = curve.add(&P, &mP);

            // P + O = P
            let PpO = curve.add(&P, &O);

            // O + P = P
            let OpP = curve.add(&O, &P);

            // O + O = P
            let OpO = curve.add(&O, &O);

            assert_eq!(PmP.equals(&O), u32::MAX);
            assert_eq!(PpO.equals(&P), u32::MAX);
            assert_eq!(OpP.equals(&P), u32::MAX);
            assert_eq!(OpO.equals(&O), u32::MAX);
        }
    }

    // Runs N random trials of the associativity test:
    //   (P + Q) + R  ==  P + (R + Q)
    // using three independently-sampled, mutually-distinct points.
    #[test]
    fn test_add_associativity_random() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_add_associativity_random".as_bytes());

        for _ in 0..TEST_ITER {
            let P = curve.rand_point(&mut rng);
            let Q = curve.rand_point(&mut rng);
            let R = curve.rand_point(&mut rng);

            // (P + Q) + R
            let pq = curve.add(&P, &Q);
            let lhs = curve.add(&pq, &R);

            // P + (R + Q)
            let rq = curve.add(&R, &Q);
            let rhs = curve.add(&P, &rq);

            assert_eq!(lhs.equals(&rhs), u32::MAX);
        }
    }

    #[test]
    fn test_add_commutativity_random() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_add_commutativity_random".as_bytes());

        for _ in 0..TEST_ITER {
            let P = curve.rand_point(&mut rng);
            let Q = curve.rand_point(&mut rng);

            // P + Q
            let pq = curve.add(&P, &Q);

            // Q + P
            let qp = curve.add(&Q, &P);

            assert_eq!(pq.equals(&qp), u32::MAX);
        }
    }

    // Check that add and dbl agree on [2]P:  P + P == dbl(P).
    #[test]
    fn test_add_equals_dbl_random() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_add_equals_dbl_random".as_bytes());

        const TRIALS: usize = 50;

        for _ in 0..TRIALS {
            let P = curve.rand_point(&mut rng);

            let via_add = curve.add(&P, &P);
            let via_dbl = curve.double(&P);

            assert_eq!(via_add.equals(&via_dbl), u32::MAX);
        }
    }

    // Ensure that the output follows:
    //   x(P + Q) = (u - v) / w
    //   x(P - Q) = (u + v) / w
    #[test]
    fn test_add_components() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_add_equals_dbl_random".as_bytes());

        const TRIALS: usize = 50;

        for _ in 0..TRIALS {
            let P = curve.rand_point(&mut rng);
            let Q = curve.rand_point(&mut rng);
            let mQ = Q.neg();

            let PpQ = curve.add(&P, &Q);
            let PmQ = curve.add(&P, &mQ);

            let (x1, _) = PpQ.to_xy();
            let (x2, _) = PmQ.to_xy();

            let (u, v, w) = curve.add_components(&P, &Q);

            let check_1 = (u - v) / w;
            let check_2 = (u + v) / w;

            assert_eq!(check_1.equals(&x1), u32::MAX);
            assert_eq!(check_2.equals(&x2), u32::MAX);
        }
    }

    #[test]
    fn test_double_iter() {
        let curve = test_curve();
        let basis = curve.torsion_basis_2e_vartime(0, &[5], 3);
        let (P, Q) = curve.lift_to_jacobian(&basis);

        let A_div_3 = curve.A / Fp2::THREE;

        let P2 = curve.double_iter(&P, 247);
        let Q2 = curve.double_iter(&Q, 247);

        assert_eq!(P2.is_zero(), 0);
        assert_eq!(Q2.is_zero(), 0);
        assert_eq!(P2.Y.is_zero(), u32::MAX);
        assert_eq!(Q2.X.is_zero(), u32::MAX);
        assert_eq!(Q2.Y.is_zero(), u32::MAX);

        let P_zero = curve.double(&P2);
        let Q_zero = curve.double(&Q2);
        assert_eq!(P_zero.is_zero(), u32::MAX);
        assert_eq!(Q_zero.is_zero(), u32::MAX);

        let P2_WS = curve.double_iter_with_ws(&P, &A_div_3, 247);
        let Q2_WS = curve.double_iter_with_ws(&Q, &A_div_3, 247);

        assert_eq!(P2_WS.equals(&P2), u32::MAX);
        assert_eq!(Q2_WS.equals(&Q2), u32::MAX);
    }

    // Ensure [n]P = O and [n+1]P = P
    #[test]
    fn test_mul_random_full_order() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_mul_random".as_bytes());

        const TRIALS: usize = 50;
        const CURVE_ORDER: [u8; 32] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 5,
        ];
        const CURVE_ORDER_PLUS_ONE: [u8; 32] = [
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 5,
        ];
        const CURVE_ORDER_BIT_LENGTH: usize = 251;

        for _ in 0..TRIALS {
            let P = curve.rand_point(&mut rng);
            let Q = curve.mul(&P, &CURVE_ORDER_PLUS_ONE, 0, CURVE_ORDER_BIT_LENGTH);
            let zero = curve.mul(&P, &CURVE_ORDER, 0, CURVE_ORDER_BIT_LENGTH);
            assert_eq!(P.equals(&Q), u32::MAX);
            assert_eq!(Q.is_zero(), 0);
            assert_eq!(zero.is_zero(), u32::MAX);
        }
    }

    // Ensure [s1] * ([s2] P) == [s2] * ([s1] P)
    #[test]
    fn test_mul_random() {
        let curve = test_curve();
        let mut rng = DRNG::from_seed("test_mul_random".as_bytes());

        const TRIALS: usize = 10;

        for _ in 0..TRIALS {
            for e in [0, 10, 20, 30] {
                let P = curve.rand_point(&mut rng);

                let mut s1 = [0; 32];
                let mut s2 = [0; 32];
                rng.fill_bytes(&mut s1);
                rng.fill_bytes(&mut s2);

                let s1_Q = curve.mul(&P, &s1, 0, 251 - e);
                let s2_s1_Q = curve.mul(&s1_Q, &s2, 0, 251 - e);

                let s2_Q = curve.mul(&P, &s2, 0, 251 - e);
                let s1_s2_Q = curve.mul(&s2_Q, &s1, 0, 251 - e);

                assert_eq!(s2_s1_Q.equals(&s1_s2_Q), u32::MAX);
            }
        }
    }
}
