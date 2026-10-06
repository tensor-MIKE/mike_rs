use crate::elliptic::jac_point::Point;

use super::basis::BasisX;
use super::curve::Curve;
use super::x_point::PointX;

use fp2::traits::Fp2 as Fp2Trait;

impl<Fp2: Fp2Trait> Curve<Fp2> {
    /// Given a potential x-coordinate, determine if it's a valid point
    /// on the curve
    fn is_on_curve(&self, x: &Fp2) -> u32 {
        let x = *x;
        let mut y = x + self.A; // y = x + A
        y *= x; // y = x^2 + A*x
        y += <Fp2>::ONE; // x^2 + A*x + 1
        y *= x; // y = x^3 + A*x^2 + x
        y.is_square()
    }

    /// Given xP and xQ as (X : Z) points, compute the point x(P ± Q) (sign is unknown).
    fn projective_difference(&self, P: &PointX<Fp2>, Q: &PointX<Fp2>) -> PointX<Fp2> {
        let mut t0 = P.X * Q.X;
        let mut t1 = P.Z * Q.Z;
        let mut bxx = t0 - t1;
        bxx.set_square(); // bxx = (PX * QX - PZ * QZ)^2
        let mut bxz = t0 + t1;
        t0 = P.X * Q.Z;
        t1 = P.Z * Q.X;
        let mut bzz = t0 + t1;
        bxz *= bzz; // (PX * QX + PZ * QZ) * (PX * QZ + PZ * QX)
        bzz = t0 - t1;
        bzz.set_square(); // bzz = (PX * QZ + PZ * QX)^2
        t0 *= t1;
        t0 *= self.A;
        t0.set_mul2();
        bxz += t0; // bzz = (PX * QX + PZ * QZ) * (PX * QZ + PZ * QX) + 2 * A * PX * QZ * PZ * QX

        // We now normalise the result by \bar{(PZ * QZ)}^2
        t0 = P.Z.conjugate();
        t1 = Q.Z.conjugate();
        t0 *= t1;
        t0.set_square();
        bxx *= t0;
        bxz *= t0;
        bzz *= t0;

        // Solve the quadratic equation
        t0 = bxz.square();
        t1 = bxx * bzz;
        t0 -= t1;
        let r = t0.set_sqrt();
        debug_assert!(r == u32::MAX);

        // Set the point (bxz + sqrt(bxz^2 - bxx*bzz) : bzz)
        PointX::new(&(bxz + t0), &bzz)
    }

    /// Given the Montgomery coefficient A which is not a square, we can
    /// compute xP as ([n]*A : 1) for some [n]. Assumes `A` is public
    fn full_even_torsion_point_from_A_vartime(&self) -> PointX<Fp2> {
        let mut x = self.A;
        while self.is_on_curve(&x) == 0 {
            x += self.A;
        }
        PointX::new(&x, &Fp2::ONE)
    }

    /// When the Montgomery coefficient A is a square, we find an element
    /// of Fp2 of the form 1 + i*b for which -A / (1 + i*b) is a valid
    /// x-coordinate on the curve and (1 + b^2) is a non-quadratic residue
    /// in the base field GF(p).
    fn full_even_torsion_point_from_nqr_vartime(&self) -> PointX<Fp2> {
        // We need to find a b such that (1 + b^2) is a not a square and
        // x = - A / (1 + i* b) is a valid x-coordinate on the curve.
        let mut is_nqr: u32 = 0;
        let mut on_curve: u32 = 0;

        let mut b = Fp2::ZERO;
        let mut z = Fp2::ONE;
        let mut h = 0;

        while on_curve == 0 {
            // Find 1 + h^2 which is not a square in GF(p)
            while is_nqr == 0 {
                b.set_x0_small(h * h + 1);
                is_nqr = !b.is_square_base_field();
                h += 1;
            }

            // We now need to determine whether -A / (1 + i*h)
            // is a point on the curve. This is the same as checking
            // whether A^2 * (z - 1) - z^2 is a non-square in GF(p^2)
            z = Fp2::ONE;
            let mut t0 = Fp2::ZERO;
            t0.set_x1_small(h - 1);
            z.set_x1_small(h - 1);

            // z is a valid coordinate providing that (A^2 * (z - 1) - z^2) is a NQR
            let t1 = self.A.square() * t0 - z.square();
            on_curve = !t1.is_square();

            // Reset is_nqr to find a new value on failure of the above.
            is_nqr = 0;
        }

        // Create a point
        PointX::new(&-self.A, &z)
    }

    /// Compute the torsion basis E\[2^e\]. For a curve with order c * 2^f, `e_diff`
    /// should be `f - e` and `cofactor` and `cofactor_bitsize` should encode `c` as
    /// an array of bytes.
    ///
    /// Warning: requires that self.A != 0
    pub fn torsion_basis_2e_vartime(
        &self,
        e_diff: usize,
        cofactor: &[u8],
        cofactor_bitsize: usize,
    ) -> BasisX<Fp2> {
        // This function cannot run with A = 0, which is not a problem for
        // MIKE as we assume all public keys have A in Fp2 (not in Fp).
        if self.A.is_zero() == u32::MAX {
            panic!("Torsion basis computation requires A != 0");
        }

        // Whether or not A is square determines how we sample P.
        // We assume that the curve itself is public, and so we can
        // branch on this outcome.
        let A_is_square = self.A.is_square();

        // Compute the point P which has full even order 2^f.
        let mut xP = if A_is_square == 0 {
            self.full_even_torsion_point_from_A_vartime()
        } else {
            self.full_even_torsion_point_from_nqr_vartime()
        };

        // We can compute another linearly independent point from xP
        // which also has full even torsion.
        let mut xPQ = if A_is_square == 0 {
            PointX::new(&-(xP.X + self.A), &Fp2::ONE)
        } else {
            PointX::new(&-(xP.X + xP.Z * self.A), &xP.Z)
        };

        // First we clear the odd cofactor using multiplication
        xP = self.xmul(&xP, cofactor, cofactor_bitsize);
        xPQ = self.xmul(&xPQ, cofactor, cofactor_bitsize);

        // We clear the 2^(f - e) order with repeated doubling.
        self.xdbl_iter_into(&mut xP, e_diff);
        self.xdbl_iter_into(&mut xPQ, e_diff);

        // Compute the difference point to get the basis x(P), x(Q) and x(P-Q)
        let xQ = self.projective_difference(&xP, &xPQ);

        // Set the basis
        BasisX::from_points(&xP, &xQ, &xPQ)
    }

    pub fn lift_to_jacobian(&self, basis: &BasisX<Fp2>) -> (Point<Fp2>, Point<Fp2>) {
        // Requires an inversion, could make a function which assumes this has been inverted?
        let PX = basis.P.x();
        let PY_sqr = PX * (PX.square() + PX * self.A + Fp2::ONE);
        let (PY, ok) = PY_sqr.sqrt();
        debug_assert!(ok == u32::MAX);
        let P = Point::new_xy(&PX, &PY);

        // TODO: clean up stack? This has just been adapted from uglier C code
        // compute QY with Okeya-Sakurai algorithm
        let (QX, QZ) = (basis.Q.X, basis.Q.Z);
        let (PQX, PQZ) = (basis.PQ.X, basis.PQ.Z);

        let v1 = PX * QZ;
        let v2 = QX + v1;
        let v3 = QX - v1;
        let v3 = v3.square();
        let v3 = v3 * PQX;
        let v1 = self.A.mul2();
        let v1 = v1 * QZ;
        let v2 = v2 + v1;
        let v4 = PX * QX;
        let v4 = v4 + QZ;
        let v2 = v2 * v4;
        let v1 = v1 * QZ;
        let v2 = v2 - v1;
        let v2 = v2 * PQZ;
        let QY = v3 - v2;
        let v1 = PY + PY;
        let v1 = v1 * QZ;
        let v1 = v1 * PQZ;
        let QX = QX * v1;
        let QZ = QZ * v1;
        let v1 = QZ * QZ;
        let QY = QY * v1;
        let QX = QX * QZ;

        let Q = Point::new(&QX, &QY, &QZ);

        (P, Q)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic::test_data::{A1, A2};
    use crate::{
        elliptic::{curve::Curve, jac_point::Point, x_point::PointX},
        fields::MikeTestFp2 as Fp2,
    };
    use fp2::traits::Fq as _;

    const COFACTOR: [u8; 1] = [5];
    const COFACTOR_BITSIZE: usize = 3;

    fn test_order_2e_xdbl(E: &Curve<Fp2>, P: &PointX<Fp2>, e: usize, mont_point: bool) {
        let tmp = E.xdbl_iter(P, e - 1);
        let zero = E.xdbl(&tmp);
        // Assert the point has order exactly 2^e
        assert_eq!(tmp.is_zero(), 0);
        assert_eq!(zero.is_zero(), u32::MAX);

        // Assert the point of 2-torsion is the point (0 : 1)
        if mont_point {
            assert_eq!(tmp.x().is_zero(), u32::MAX);
        }
    }

    fn test_order_2e_jac(E: &Curve<Fp2>, P: &Point<Fp2>, e: usize, mont_point: bool) {
        let tmp = E.double_iter(P, e - 1);
        let zero = E.double(&tmp);
        // Assert the point has order exactly 2^e
        assert_eq!(tmp.is_zero(), 0);
        assert_eq!(zero.is_zero(), u32::MAX);

        // Assert the point of 2-torsion is the point (0 : 0 : 1)
        if mont_point {
            let (x, y) = tmp.to_xy();
            assert_eq!(x.is_zero(), u32::MAX);
            assert_eq!(y.is_zero(), u32::MAX);
        }
    }

    fn test_torsion_basis(A: &Fp2) {
        let (curve, check) = Curve::<Fp2>::new(A);
        assert_eq!(check, u32::MAX);

        for e in [246, 200, 100] {
            let basis = curve.torsion_basis_2e_vartime(248 - e, &COFACTOR, COFACTOR_BITSIZE);

            // Test the x-only basis
            test_order_2e_xdbl(&curve, &basis.P, e, false);
            test_order_2e_xdbl(&curve, &basis.Q, e, true);
            test_order_2e_xdbl(&curve, &basis.PQ, e, false);

            // Test the Jacobian points
            let (P, Q) = curve.lift_to_jacobian(&basis);
            test_order_2e_jac(&curve, &P, e, false);
            test_order_2e_jac(&curve, &Q, e, true);
        }
    }

    #[test]
    fn test_torsion_basis_qr() {
        assert_eq!(A1.is_square(), u32::MAX);
        test_torsion_basis(&A1)
    }

    #[test]
    fn test_torsion_basis_nqr() {
        assert_eq!(A2.is_square(), 0);
        test_torsion_basis(&A2)
    }
}
