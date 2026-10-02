use fp2::traits::Fq as FqTrait;

/// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
/// (special case of a Montgomery curve).
#[derive(Clone, Copy, Debug)]
pub struct Curve<Fq: FqTrait> {
    pub A: Fq,   // A
    pub A24: Fq, // (A+2)/4
}

impl<Fq: FqTrait> Curve<Fq> {
    /// Create a new curve instance, with the provided constant.
    pub fn new(A: &Fq) -> (Self, u32) {
        let a = *A;
        let a_plus_two = a + <Fq>::TWO;
        let a24 = a_plus_two.half().half();

        // We check that the curve is not singular, i.e. A^2 != 4.
        let eq_two = a.equals(&<Fq>::TWO);
        let eq_neg_two = a_plus_two.is_zero();
        let valid = !(eq_two | eq_neg_two);

        let curve = Self { A: a, A24: a24 };
        (curve, valid)
    }

    /// Create a new curve instance with no checks, works at compile time
    pub const fn new_no_checks(A: &Fq, A24: &Fq) -> Self {
        Self { A: *A, A24: *A24 }
    }

    /// Compute a curve from the projective coordinates of (A + 2) / 4 = (A24 : C24)
    #[inline]
    pub fn curve_from_A24_proj(A24: &Fq, C24: &Fq) -> Self {
        // Compute A from (A24 : C24)
        let a24 = (*A24) / (*C24);
        let a = a24.mul4() - Fq::TWO;
        Self::new_no_checks(&a, &a24)
    }

    /// Compute the j-invariant of the curve.
    pub fn j_invariant(&self) -> Fq {
        let mut j = self.A.square();
        let mut t1 = Fq::ONE; // This should be C^2, but C = 1
        let mut t0 = Fq::TWO; // This should be 2C^2
        t0 = j - t0;
        t0 -= t1;
        j = t0 - t1;
        t1.set_square();
        j *= t1;
        t0.set_mul4();
        t1 = t0.square();
        t0 *= t1;
        t0.set_mul4();
        j.set_invert();
        j *= t0;

        j
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for Curve<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Elliptic Curve: y^2 = x^3 + ({})*x^2 + x", self.A)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::MikeTestFp as Fp;

    #[test]
    fn test_new_curve_valid() {
        let a = Fp::ZERO;
        let (_, valid) = Curve::new(&a);
        assert_eq!(valid, u32::MAX);
    }

    #[test]
    fn test_new_curve_singular_pos_two() {
        let a = Fp::TWO;
        let (_, valid) = Curve::new(&a);
        assert_eq!(valid, 0);
    }

    #[test]
    fn test_new_curve_singular_neg_two() {
        let a = -Fp::TWO;
        let (_, valid) = Curve::new(&a);
        assert_eq!(valid, 0);
    }
}
