use fp2::traits::Fq as FqTrait;

/// Theta Point in with a projective representation for dimension 2
#[derive(Clone, Copy, Debug)]
pub struct ThetaPointDimTwo<Fp: FqTrait> {
    x: Fp,
    y: Fp,
    z: Fp,
    t: Fp,
}

impl<Fq: FqTrait> ThetaPointDimTwo<Fq> {
    /// Compile time, create a new theta point from Fq elements
    pub const fn new(x: &Fq, y: &Fq, z: &Fq, t: &Fq) -> Self {
        Self {
            x: *x,
            y: *y,
            z: *z,
            t: *t,
        }
    }

    /// Recover the coordinates of the element
    pub fn coords(&self) -> (Fq, Fq, Fq, Fq) {
        (self.x, self.y, self.z, self.t)
    }

    /// Compute the Hadamard transformation of the point in place
    /// Cost: 8a
    #[inline(always)]
    pub fn set_hadamard(&mut self) {
        let t1 = self.x + self.y;
        let t2 = self.x - self.y;
        let t3 = self.z + self.t;
        let t4 = self.z - self.t;

        self.x = t1 + t3;
        self.y = t2 + t4;
        self.z = t1 - t3;
        self.t = t2 - t4;
    }

    /// Return the Hadamard transform of this point
    /// Cost: 8a
    pub fn hadamard(&self) -> Self {
        let mut P = *self;
        P.set_hadamard();
        P
    }

    /// Compute the Scholten change of basis transformation of the point in place
    /// Cost: 8a
    #[inline(always)]
    pub fn set_scholten_change_of_basis(&mut self) {
        // xx = x + y + z - t
        // yy = x - y + z + t
        // zz = x + y - z + t
        // ww = x - y - z - t
        let t1 = self.x + self.y;
        let t2 = self.x - self.y;
        let t3 = self.z + self.t;
        let t4 = self.z - self.t;

        self.x = t1 + t4;
        self.y = t2 + t3;
        self.z = t1 - t4;
        self.t = t2 - t3;
    }

    /// Compute the squaring transform of the point in place
    /// Cost: 4s
    #[inline(always)]
    pub fn set_square(&mut self) {
        self.x.set_square();
        self.y.set_square();
        self.z.set_square();
        self.t.set_square();
    }

    /// Compute H(S(P)) of the point in place
    /// Cost: 4s + 8a
    pub fn set_square_hadamard(&mut self) {
        self.set_square();
        self.set_hadamard();
    }

    /// Return H(S(P)) of this point
    /// Cost: 4s + 8a
    pub fn square_hadamard(&self) -> Self {
        let mut P = *self;
        P.set_square_hadamard();
        P
    }

    /// Multiply each coordinate of self by the coordinates of other in place
    /// Cost: 4a
    pub fn set_pointwise_sub(&mut self, other: &Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
        self.t -= other.t;
    }

    /// Compute the coordinate wise multiplication of two points `P \star Q`
    /// Cost: 4a
    pub fn pointwise_sub(&self, other: &Self) -> Self {
        let mut P = *self;
        P.set_pointwise_sub(other);
        P
    }

    /// Multiply each coordinate of self by the coordinates of other in place
    /// Cost: 4m
    pub fn set_pointwise_mul(&mut self, other: &Self) {
        self.x *= other.x;
        self.y *= other.y;
        self.z *= other.z;
        self.t *= other.t;
    }

    /// Returns if any of the coordinates of `self` are zero
    pub fn has_zero_coordinate(&self) -> u32 {
        self.x.is_zero() | self.y.is_zero() | self.z.is_zero() | self.t.is_zero()
    }

    /// Conditionally swap two points
    #[inline]
    pub fn cond_swap(P: &mut Self, Q: &mut Self, ctl: u32) {
        Fq::cond_swap(&mut P.x, &mut Q.x, ctl);
        Fq::cond_swap(&mut P.y, &mut Q.y, ctl);
        Fq::cond_swap(&mut P.z, &mut Q.z, ctl);
        Fq::cond_swap(&mut P.t, &mut Q.t, ctl);
    }
}

/// Default element used for initialisation
impl<Fq: FqTrait> Default for ThetaPointDimTwo<Fq> {
    fn default() -> Self {
        Self {
            x: Fq::ZERO,
            y: Fq::ZERO,
            z: Fq::ZERO,
            t: Fq::ZERO,
        }
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for ThetaPointDimTwo<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        writeln!(f, "ThetaPointDimTwo: (")?;
        writeln!(f, "    {},", self.x)?;
        writeln!(f, "    {},", self.y)?;
        writeln!(f, "    {},", self.z)?;
        writeln!(f, "    {},", self.t)?;
        write!(f, ")")
    }
}

/// Theta Structure for dimension 2
#[derive(Clone, Copy, Debug)]
pub struct ThetaStructureDimTwo<Fp: FqTrait> {
    null_point: ThetaPointDimTwo<Fp>,
}

impl<Fq: FqTrait> ThetaStructureDimTwo<Fq> {
    pub fn new(null_point: ThetaPointDimTwo<Fq>) -> Self {
        Self { null_point }
    }

    /// Compute projective precomputed ratios used in differential addition.
    fn arithmetic_precomputation(&self) -> [Fq; 4] {
        // Compute projectively A^2/B^2 = A^4*C^2*D^2, etc.
        let (AA, BB, CC, DD) = self.null_point.square_hadamard().coords();
        let t1 = AA * BB;
        let t2 = CC * DD;
        let T0 = t1 * CC;
        let Z0 = t1 * DD;
        let Y0 = t2 * AA;
        let X0 = t2 * BB;

        [X0, Y0, Z0, T0]
    }

    /// Given P, Q and P - Q, compute P + Q
    pub fn diff_addition(
        &self,
        P: &ThetaPointDimTwo<Fq>,
        Q: &ThetaPointDimTwo<Fq>,
        PQ: &ThetaPointDimTwo<Fq>,
    ) -> ThetaPointDimTwo<Fq> {
        let [X0, Y0, Z0, T0] = self.arithmetic_precomputation();

        // H(S(P)) and H(S(Q)) — cost: 8S + 16a each
        let (p1, p2, p3, p4) = P.square_hadamard().coords();
        let (q1, q2, q3, q4) = Q.square_hadamard().coords();

        // Pointwise multiply with precomputed ratios — cost: 7M
        let xp = X0 * p1 * q1;
        let yp = Y0 * (p2 * q2);
        let zp = Z0 * (p3 * q3);
        let tp = T0 * (p4 * q4);

        // Hadamard — cost: 8a
        // TODO: make a hadamard on coordinates itself?
        let mut tmp = ThetaPointDimTwo::new(&xp, &yp, &zp, &tp);
        tmp.set_hadamard();
        let (X, Y, Z, T) = tmp.coords();

        // Replace four divisions with 10 multiplications
        let (PQx, PQy, PQz, PQt) = PQ.coords();
        let PQxy = PQx * PQy;
        let PQzt = PQz * PQt;

        let X = X * PQzt * PQy;
        let Y = Y * PQzt * PQx;
        let Z = Z * PQxy * PQt;
        let T = T * PQxy * PQz;

        ThetaPointDimTwo::new(&X, &Y, &Z, &T)
    }

    /// Given the 8-torsion above the kernel, compute the codomain of the
    /// (2,2)-isogeny and the image of all points in `image_points`.
    /// Cost:
    ///   Codomain: 8S + 9M
    ///   Image: 4S + 4M per point
    pub fn two_isogeny(
        T1: &ThetaPointDimTwo<Fq>,
        T2: &ThetaPointDimTwo<Fq>,
        image_points: &mut [ThetaPointDimTwo<Fq>],
        hadamard: [bool; 2],
    ) -> Self {
        // Conditionally hadamard to move to the dual
        let (T1, T2) = if hadamard[0] {
            (T1.hadamard(), T2.hadamard())
        } else {
            (*T1, *T2)
        };

        // coordinates of H(S(P))
        let (xA, xB, _, _) = T1.square_hadamard().coords();
        let (zA, tB, zC, tD) = T2.square_hadamard().coords();

        let xAtB = xA * tB;
        let zAxB = zA * xB;
        let zCtD = zC * tD;

        let A = zA * xAtB;
        let B = tB * zAxB;
        let C = zC * xAtB;
        let D = tD * zAxB;
        let mut codomain = ThetaPointDimTwo::new(&A, &B, &C, &D);

        let A_inv = xB * zCtD;
        let B_inv = xA * zCtD;
        let C_inv = D;
        let D_inv = C;
        let inverses = ThetaPointDimTwo::new(&A_inv, &B_inv, &C_inv, &D_inv);

        if hadamard[1] {
            codomain.set_hadamard();
        }

        for P in image_points.iter_mut() {
            if hadamard[0] {
                P.set_hadamard()
            }
            P.set_square_hadamard();
            P.set_pointwise_mul(&inverses);
            if hadamard[1] {
                P.set_hadamard();
            }
        }

        Self {
            null_point: codomain,
        }
    }
}
