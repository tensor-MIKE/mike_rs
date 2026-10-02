use crate::theta::dim_four::ThetaStructure;
use fp2::traits::Fq as FqTrait;

impl<Fq: FqTrait> ThetaStructure<Fq> {
    /// Compute the special modular invariants
    // 8a + 4S + 4a + 7a + 9a + 1I + 9M
    // = 1I + 9M + 4S + 28a
    pub fn special_modular_invariants(&self) -> [Fq; 9] {
        let null = self.null_point();

        let x0 = null[0];
        let x3 = null[3];
        let x5 = null[5];
        let x6 = null[6];
        let x9 = null[9];
        let x15 = null[15];

        // Hadamard of x0, x6, x9, x15: 8a
        let (t0, t1) = (x0 + x6, x0 - x6);
        let (t2, t3) = (x9 + x15, x9 - x15);
        let (u0, u2) = (t0 + t2, t0 - t2);
        let (u1, u3) = (t1 + t3, t1 - t3);

        // 4S
        let s0 = u0.square();
        let s1 = u1.square();
        let s2 = u2.square();
        let s3 = u3.square();

        // 4a
        let t = (x3 * x5).mul4().mul4();

        // Hadamard of s0, s1, s2, s3: 7a
        let (t0, t1) = (s0 + s1, s0 - s1);
        let (t2, t3) = (s2 + s3, s2 - s3);
        let u2 = t0 - t2;
        let (u1, u3) = (t1 + t3, t1 - t3);

        // 6a
        let (s4, s5) = (u1 + t, u1 - t);
        let (s6, s7) = (u2 + t, u2 - t);
        let (s8, s9) = (u3 + t, u3 - t);

        // 1I + 9M
        let theta_0_inv = s0.invert();
        let mut theta_n = [s1, s2, s3, s4, s5, s6, s7, s8, s9];
        for t in theta_n.iter_mut() {
            *t *= theta_0_inv;
        }

        // 3a
        theta_n[0] = theta_n[0].mul2();
        theta_n[1] = theta_n[1].mul2();
        theta_n[2] = theta_n[2].mul2();

        theta_n
    }

    /// Compute the absolute modular invariants
    // Cost: 1I + 27M + 22S + 60a
    pub fn absolute_modular_invariants(&self) -> [Fq; 4] {
        // Compute the (normalised) invariants s1...s9
        // Cost: 1I + 9M + 4S + 28a
        let theta_n = self.special_modular_invariants();

        // Compute si^j for i in range 1..9 and j in range 1..5
        // Cost 9 * (2M + 2S)
        let powers: [(Fq, Fq, Fq, Fq); 9] = theta_n.map(|s| {
            let t = s.square(); // s^2
            let r = s * t; // s^3
            let u = t.square(); // s^4
            let v = s * u; // s^5
            (t, r, u, v)
        });

        // Sum up the components for each invariant
        // Cost: 4 * 8a
        let (mut j2, mut j3, mut j4, mut j5) = powers[0];
        for &(t, r, u, v) in &powers[1..] {
            j2 += t;
            j3 += r;
            j4 += u;
            j5 += v;
        }

        // Return the invariants
        [j2, j3, j4, j5]
    }
}

#[cfg(test)]
mod tests {
    use super::ThetaStructure;
    use crate::theta::test_data::{OC_NULL, OD_NULL};

    #[test]
    fn test_special_invariants() {
        let oc = ThetaStructure::new_from_null_point(OC_NULL);
        let od = ThetaStructure::new_from_null_point(OD_NULL);
        let secret_c = oc.absolute_modular_invariants();
        let secret_d = od.absolute_modular_invariants();

        for i in 0..4 {
            let a = secret_c[i];
            let b = secret_d[i];
            assert_eq!(a.equals(&b), u32::MAX);
        }
    }
}
