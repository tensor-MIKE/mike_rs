use super::{basis::BasisX, curve::Curve, x_point::PointX};
use fp2::traits::Fq as FqTrait;

impl<Fq: FqTrait> Curve<Fq> {
    // ============================================================================
    // x-only doubling methods
    // ============================================================================

    /// P <- [2]*P (x-only variant) in place.
    /// Cost: 2S + 3M.
    #[inline(always)]
    pub fn xdbl_into(&self, P: &mut PointX<Fq>) {
        let mut V1 = (P.X + P.Z).square();
        let V2 = (P.X - P.Z).square();
        P.X = V1 * V2;
        V1 -= V2;
        P.Z = V1;
        P.Z *= self.A24;
        P.Z += V2;
        P.Z *= V1;
    }

    /// Return [2]*P (x-only variant).
    /// Cost: 2S + 3M.
    pub fn xdbl(&self, P: &PointX<Fq>) -> PointX<Fq> {
        let mut Q = *P;
        self.xdbl_into(&mut Q);
        Q
    }

    /// P <- [2^n]*P (x-only variant) in place.
    /// Cost: n * (2S + 3M).
    pub fn xdbl_iter_into(&self, P: &mut PointX<Fq>, n: usize) {
        for _ in 0..n {
            self.xdbl_into(P);
        }
    }

    /// Return [2^n]*P (x-only variant).
    /// Cost: n * (2S + 3M).
    pub fn xdbl_iter(&self, P: &PointX<Fq>, n: usize) -> PointX<Fq> {
        let mut Q = *P;
        self.xdbl_iter_into(&mut Q, n);
        Q
    }

    /// P <- [2]*P (x-only variant) in place using projective curve constant (A+2)/4 = (A24 : C24).
    /// Cost: 2S + 4M.
    #[inline(always)]
    pub fn xdbl_proj_into(P: &mut PointX<Fq>, A24: &Fq, C24: &Fq) {
        let mut t0 = P.X + P.Z;
        t0.set_square();
        let mut t1 = P.X - P.Z;
        t1.set_square();
        let t2 = t0 - t1;
        t1 *= *C24;
        P.X = t0 * t1;
        t0 = t2 * (*A24);
        t0 += t1;
        P.Z = t0 * t2;
    }

    /// P <- [2^n]*P (x-only variant) in place using projective curve constant (A+2)/4 = (A24 : C24).
    /// Cost: n * (2S + 4M).
    #[inline]
    pub fn xdbl_proj_iter_into(P: &mut PointX<Fq>, A24: &Fq, C24: &Fq, n: usize) {
        for _ in 0..n {
            Self::xdbl_proj_into(P, A24, C24);
        }
    }

    // ============================================================================
    // x-only differential addition methods
    // ============================================================================

    /// Q <- x(P + Q) in place, given x(P), x(Q) and x(P - Q) as projective points.
    /// Cost: 2S + 4M.
    #[inline(always)]
    pub fn xdiff_add_into(P: &PointX<Fq>, Q: &mut PointX<Fq>, PmQ: &PointX<Fq>) {
        let V1 = (P.X - P.Z) * (Q.X + Q.Z);
        let V2 = (P.X + P.Z) * (Q.X - Q.Z);
        Q.X = PmQ.Z * (V1 + V2).square();
        Q.Z = PmQ.X * (V1 - V2).square();
    }

    /// Return x(P + Q) given x(P), x(Q) and x(P - Q) as projective points.
    /// Cost: 2S + 4M.
    pub fn xdiff_add(P: &PointX<Fq>, Q: &PointX<Fq>, PmQ: &PointX<Fq>) -> PointX<Fq> {
        let mut R = *Q;
        Self::xdiff_add_into(P, &mut R, PmQ);
        R
    }

    /// Q <- x(P + Q) in place, given x(P), x(Q) and the affine x-coordinate x(P - Q).
    /// Cost: 2S + 3M.
    #[inline(always)]
    pub fn xdiff_add_aff_into(P: &PointX<Fq>, Q: &mut PointX<Fq>, xPmQ: &Fq) {
        let V1 = (P.X - P.Z) * (Q.X + Q.Z);
        let V2 = (P.X + P.Z) * (Q.X - Q.Z);
        Q.X = (V1 + V2).square();
        Q.Z = *xPmQ * (V1 - V2).square();
    }

    /// Return x(P + Q) given x(P), x(Q) and the affine x-coordinate x(P - Q).
    /// Cost: 2S + 3M.
    #[inline]
    pub fn xdiff_add_aff(P: &PointX<Fq>, Q: &PointX<Fq>, xPmQ: &Fq) -> PointX<Fq> {
        let mut R = *Q;
        Self::xdiff_add_aff_into(P, &mut R, xPmQ);
        R
    }

    // ============================================================================
    // x-only double and add methods
    // ============================================================================

    /// P <- [2]*P, Q <- x(P + Q) in place, given x(P), x(Q) and x(P - Q) as projective points.
    /// Cost: 4S + 7M.
    #[inline(always)]
    pub fn xdbladd_into(&self, P: &mut PointX<Fq>, Q: &mut PointX<Fq>, PmQ: &PointX<Fq>) {
        let t0 = P.X + P.Z;
        let t1 = P.X - P.Z;
        let X2P_sq = t0.square();
        let Z2P_sq = t1.square();

        let t2 = X2P_sq - Z2P_sq;
        P.X = X2P_sq * Z2P_sq;
        P.Z = t2 * (Z2P_sq + self.A24 * t2);

        let t0_Q = t0 * (Q.X - Q.Z);
        let t1_Q = t1 * (Q.X + Q.Z);

        Q.X = PmQ.Z * (t0_Q + t1_Q).square();
        Q.Z = PmQ.X * (t0_Q - t1_Q).square();
    }

    /// P <- [2]*P, Q <- x(P + Q) in place, given x(P), x(Q) and the affine x-coordinate x(P - Q).
    /// Cost: 4S + 6M.
    #[inline(always)]
    pub fn xdbladd_aff_into(&self, P: &mut PointX<Fq>, Q: &mut PointX<Fq>, xPmQ: &Fq) {
        let t0 = P.X + P.Z;
        let t1 = P.X - P.Z;
        let X2P_sq = t0.square();
        let Z2P_sq = t1.square();

        let t2 = X2P_sq - Z2P_sq;
        P.X = X2P_sq * Z2P_sq;
        P.Z = t2 * (Z2P_sq + self.A24 * t2);

        let t0_Q = t0 * (Q.X - Q.Z);
        let t1_Q = t1 * (Q.X + Q.Z);

        Q.X = (t0_Q + t1_Q).square();
        Q.Z = *xPmQ * (t0_Q - t1_Q).square();
    }

    /// Return ([2]*P, x(P + Q)) given x(P), x(Q) and x(P - Q) as projective points.
    /// Cost: 4S + 7M.
    #[inline]
    pub fn xdbladd(
        &self,
        P: &PointX<Fq>,
        Q: &PointX<Fq>,
        PmQ: &PointX<Fq>,
    ) -> (PointX<Fq>, PointX<Fq>) {
        let mut R = *P;
        let mut S = *Q;
        self.xdbladd_into(&mut R, &mut S, PmQ);
        (R, S)
    }

    // ============================================================================
    // x-only Montgomery ladder for scalar multiplication
    // ============================================================================

    /// P3 <- [n]*P (x-only variant) in place.
    /// Integer n is encoded as unsigned little-endian with length nbitlen bits.
    /// Bits beyond nbitlen are ignored.
    pub fn xmul_into(&self, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: &[u8], nbitlen: usize) {
        // Montgomery ladder: see https://eprint.iacr.org/2017/212
        if nbitlen == 0 {
            P3.X = Fq::ONE;
            P3.Z = Fq::ZERO;
            return;
        }

        let mut X0 = PointX::INFINITY;
        let mut X1 = *P;
        let mut cc = 0u32;
        if nbitlen > 21 {
            // If n is large enough then it is worthwhile to
            // normalize the source point to affine.
            // If P = inf, then this sets Xp to 0; thus, the
            // output of both xdbl() and xadd_aff() has Z = 0,
            // so we correctly get the point-at-infinity at the end.
            let Xp = P.X / P.Z;
            for i in (0..nbitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_aff_into(&mut X0, &mut X1, &Xp);
                cc = ctl;
            }
        } else {
            for i in (0..nbitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_into(&mut X0, &mut X1, P);
                cc = ctl;
            }
        }
        PointX::cond_swap(&mut X0, &mut X1, cc);

        // The ladder may fail if P = (0,0) (which is a point of
        // order 2) because in that case xadd() (and xadd_aff())
        // return Z = 0 systematically, so the result is considered
        // to be the point-at-infinity, which is wrong is n is odd.
        // We adjust the result in that case.
        let spec = P.X.is_zero() & !P.Z.is_zero() & ((n[0] as u32) & 1).wrapping_neg();
        *P3 = X0;
        P3.set_cond(&PointX::INFINITY, spec);
    }

    /// Return [n]*P (x-only variant).
    /// Integer n is encoded as unsigned little-endian with length nbitlen bits.
    /// Bits beyond nbitlen are ignored.
    pub fn xmul(&self, P: &PointX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        self.xmul_into(&mut P3, P, n, nbitlen);
        P3
    }

    /// P3 <- [n]*P (x-only variant) in place.
    /// The scalar n is a u64 and is assumed to be a public value (variable-time).
    pub fn xmul_u64_vartime_into(&self, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: u64) {
        // Handle small cases.
        match n {
            0 => {
                *P3 = PointX::INFINITY;
            }
            1 => {
                *P3 = *P;
            }
            2 => *P3 = self.xdbl(P),
            _ => {
                let nbitlen = (64 - n.leading_zeros()) as usize;

                let mut X0 = PointX::INFINITY;
                let mut X1 = *P;
                let mut cc = 0u32;
                if nbitlen > 21 {
                    // If n is large enough then it is worthwhile to
                    // normalize the source point to affine.
                    // If P = inf, then this sets Xp to 0; thus, the
                    // output of both xdbl() and xadd_aff() has Z = 0,
                    // so we correctly get the point-at-infinity at the end.
                    let Xp = P.X / P.Z;
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                        self.xdbladd_aff_into(&mut X0, &mut X1, &Xp);
                        cc = ctl;
                    }
                } else {
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        PointX::cond_swap(&mut X0, &mut X1, ctl ^ cc);
                        self.xdbladd_into(&mut X0, &mut X1, P);
                        cc = ctl;
                    }
                }
                PointX::cond_swap(&mut X0, &mut X1, cc);

                // The ladder may fail if P = (0,0) (which is a point of
                // order 2) because in that case xadd() (and xadd_aff())
                // return Z = 0 systematically, so the result is considered
                // to be the point-at-infinity, which is wrong is n is odd.
                // We adjust the result in that case.
                let spec = P.X.is_zero() & !P.Z.is_zero() & (((n & 1) as u32) & 1).wrapping_neg();
                *P3 = X0;
                P3.set_cond(&PointX::INFINITY, spec);
            }
        }
    }

    /// Return [n]*P (x-only variant).
    /// The scalar n is a u64 and is assumed to be a public value (variable-time).
    pub fn xmul_u64_vartime(&self, P: &PointX<Fq>, n: u64) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        self.xmul_u64_vartime_into(&mut P3, P, n);
        P3
    }

    /// Return ([2^e]*P, [2^e]*Q, [2^e]*(P - Q)) for the x-only basis B = (P, Q, P - Q).
    pub fn basis_double_iter(&self, B: &BasisX<Fq>, e: usize) -> BasisX<Fq> {
        // TODO: if we unrolled these xdbl_iter we might save some cycles but the
        // self.xdbl_triple_iter method would be pretty ugly!
        let P = self.xdbl_iter(&B.P, e);
        let Q = self.xdbl_iter(&B.Q, e);
        let PQ = self.xdbl_iter(&B.PQ, e);
        BasisX::from_points(&P, &Q, &PQ)
    }

    // ============================================================================
    // x-only specialised ladders
    // ============================================================================

    /// Return x(P + [n]*Q) given the x-only basis B = (x(P), x(Q), x(P - Q)).
    /// Integer n is encoded as unsigned little-endian with length nbitlen bits.
    /// Bits beyond nbitlen are ignored.
    pub fn three_point_ladder(&self, B: &BasisX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        if nbitlen == 0 {
            return B.P;
        }

        // Extract out the coordinates from the basis
        let mut X0 = B.Q;
        let mut X1 = B.P;
        let mut X2 = B.PQ;

        let mut cc = 0u32;
        for i in 0..nbitlen {
            let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
            PointX::cond_swap(&mut X1, &mut X2, ctl ^ cc);
            self.xdbladd_into(&mut X0, &mut X2, &X1);
            cc = ctl;
        }
        PointX::cond_swap(&mut X1, &mut X2, cc);
        X1
    }
}
