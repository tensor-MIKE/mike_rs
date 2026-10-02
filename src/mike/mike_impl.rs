use crate::elliptic::curve::Curve;
use crate::elliptic::{basis::BasisX, jac_point::Point};
use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;
use nist_drbg_rs::{AesCtr256Drbg, Policy};
use rand_core::{TryCryptoRng, TryRngCore};
use sha3::{Digest, Sha3_256};

/// Constants associated to the field used in Mike
#[derive(Clone, Copy, Debug)]
pub struct MikeConstants {
    f: usize, // The base field has characteristic p = c*2^f - 1 for odd c
    c: usize, // The base field has characteristic p = c*2^f - 1 for odd c
    e: usize, // The length of the isogeny computed, e <= f - 2
}

impl MikeConstants {
    pub const fn new(f: usize, c: usize, e: usize) -> Self {
        Self { f, c, e }
    }
}

/// Precomputed parameters for Mike which stores the starting curve and a basis on this curve
#[derive(Clone, Copy, Debug)]
pub struct MikeParameters<Fq: Fp2Trait> {
    constants: MikeConstants,
    starting_curve: Curve<Fq>,
    basis: BasisX<Fq>,
}

impl<Fq: Fp2Trait> MikeParameters<Fq> {
    /// Create an instance of MikeParameters with (A, A24) for the curve and
    /// (x(P), x(Q) and x(P - Q)) for E[2^(e + 2)] = <P, Q>
    pub const fn new(
        constants: MikeConstants,
        starting_curve: Curve<Fq>,
        basis: BasisX<Fq>,
    ) -> Self {
        Self {
            constants,
            starting_curve,
            basis,
        }
    }
}

/// Errors for Mike runtime. Public keys are validated during secret generation
/// and RNG failure is captured during key generation
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MikeError {
    InvalidPublicKeyEncoding,
    UnsupportedPublicKey,
    NonSupersingularPublicKey,
    NonCanonicalPublicKey,
    RngFailure,
}

/// The public key of Mike is a starting curve from which to compute the
/// dimension four isogeny
pub struct MikePublicKey<Fq: Fp2Trait> {
    pub curve: Curve<Fq>,
}

impl<Fq: Fp2Trait> MikePublicKey<Fq> {
    pub fn new(curve: Curve<Fq>) -> Self {
        Self { curve }
    }

    /// An encoded Mike public key is represented by the Montgomery coefficient
    /// A for a curve E : y^2 = x^3 + Ax + x
    pub fn encode(&self) -> [u8; Fq::ENCODED_LENGTH] {
        self.curve.A.encode()
    }

    /// A public key is valid if the input bytes:
    ///
    /// 1. Represent a valid encoding (canonical within [0, p-1]) of an element Fp2
    /// 2. Represent a non-zero value of Fp2
    /// 3. Constructs an Elliptic curve in the Montgomery model, i.e. A^2 - 4 != 0
    ///
    /// NOTE: the constructed curve is *not* checked to be supersingular here,
    /// but this is instead proven at run-time by ensuring the canonical basis sampled
    /// is of the expected order
    pub fn decode(bytes: &[u8; Fq::ENCODED_LENGTH]) -> Result<Self, MikeError> {
        let (A, check) = Fq::decode(bytes);

        // If check is zero, then the value supplied was non-canonical and is rejected
        if check == 0 {
            return Err(MikeError::InvalidPublicKeyEncoding);
        }

        // An honestly generated public key from E0 has A in Fp with probability O(1/sqrt(p).
        // As the current implementation of gluing fails for these cases, we mark these
        // public keys as invalid.
        if A.x1().is_zero() == u32::MAX {
            return Err(MikeError::UnsupportedPublicKey);
        }

        // Decode the curve, here we only check whether the curve is non-singular
        // Supersingularity testing is performed during secret generation.
        let (curve, valid) = Curve::new(&A);
        valid
            .eq(&u32::MAX)
            .then_some(Self { curve })
            .ok_or(MikeError::InvalidPublicKeyEncoding)
    }
}

/// The secret key of Mike is an array of `N` bytes which represent the
/// secret scalar in little endian.
pub struct MikePrivateKey<Fp: FqTrait, Fp2: Fp2Trait, const N: usize> {
    _phantom_fp: core::marker::PhantomData<Fp>,
    _phantom_fp2: core::marker::PhantomData<Fp2>,
    constants: MikeConstants,
    secret: [u8; N],
}

impl<Fp: FqTrait, Fp2: Fp2Trait, const N: usize> MikePrivateKey<Fp, Fp2, N> {
    pub fn new(secret: [u8; N], constants: MikeConstants) -> Self {
        Self {
            _phantom_fp: core::marker::PhantomData,
            _phantom_fp2: core::marker::PhantomData,
            secret,
            constants,
        }
    }
    pub fn encode(&self) -> [u8; N] {
        self.secret
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Mike<Fp: FqTrait, Fp2: Fp2Trait, const N: usize> {
    _phantom_fp: core::marker::PhantomData<Fp>,
    params: MikeParameters<Fp2>,
}

impl<Fp: FqTrait, Fp2: Fp2Trait, const N: usize> Mike<Fp, Fp2, N> {
    /// Construct a `Mike` instance from a parameter set.
    pub const fn new(params: &MikeParameters<Fp2>) -> Self {
        Self {
            params: *params,
            _phantom_fp: core::marker::PhantomData,
        }
    }

    /// Clamp randomly sampled bytes to ensure the corresponding
    /// integer is congruent to 3 mod 8 with exactly e bits.
    fn clamp_scalar_bytes(&self, scalar: &mut [u8; N]) {
        // First we clamp the top-bits to ensure that the scalar
        // has exactly e bits
        let remainder = (self.params.constants.e) % 8;
        if remainder != 0 {
            scalar[N - 1] &= (1 << remainder) - 1;
        }

        // Then we clamp the bottom byte to ensure that x = 3 mod 8
        scalar[0] = (scalar[0] & !7) | 3;
    }

    /// Using a user supplied randomness source, sample and then clamp the secret key `x`
    fn generate_scalar<R: TryRngCore + TryCryptoRng>(
        &self,
        rng: &mut R,
    ) -> Result<[u8; N], MikeError> {
        // Populate the scalar with bytes from the rng
        let mut scalar = [0u8; N];
        rng.try_fill_bytes(&mut scalar)
            .map_err(|_| MikeError::RngFailure)?;
        self.clamp_scalar_bytes(&mut scalar);

        Ok(scalar)
    }

    /// Compute the codomain of the isogeny with kernel
    /// 4 (P + [x] Q) using the precomputed basis E0[2^e] = <P, Q>
    fn public_key(&self, scalar: &[u8]) -> MikePublicKey<Fp2>
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let E0 = self.params.starting_curve;

        // Compute P + [x]Q with x-only arithmetic
        let kernel = E0.three_point_ladder(&self.params.basis, scalar, self.params.constants.e);

        // Compute the codomain
        // NOTE: <P, Q> = E[2^(e + 2)] so this point is really the point of four torsion above the kernel
        // two_isogeny_chain is a specialised function which expects this and uses the four torsion at the
        // end of the chain to compute the randomised codomain without the requirement of an expensive sqrt
        let codomain = E0.two_isogeny_chain(&kernel, self.params.constants.e);

        MikePublicKey::new(codomain)
    }

    /// Compute the secret key given a randomly sampled scalar
    fn secret_key(&self, secret: &[u8; N]) -> MikePrivateKey<Fp, Fp2, N> {
        MikePrivateKey::new(*secret, self.params.constants)
    }

    fn keygen_impl(
        &self,
        scalar_bytes: &[u8; N],
    ) -> (MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>)
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let public_key = self.public_key(scalar_bytes);
        let secret_key = self.secret_key(scalar_bytes);

        (public_key, secret_key)
    }

    /// Generate a random keypair using user supplied entropy source
    pub fn keygen<R: TryRngCore + TryCryptoRng>(
        &self,
        rng: &mut R,
    ) -> Result<(MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>), MikeError>
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let scalar = self.generate_scalar(rng)?;
        Ok(self.keygen_impl(&scalar))
    }

    /// Generate a keypair using a AES-CTR-256 DRBG seeded with the provided value.
    /// Used for deterministic generation / KAT testing.
    pub fn keygen_seeded(
        &self,
        seed: &[u8],
    ) -> Result<(MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>), MikeError>
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let mut rng =
            AesCtr256Drbg::new(seed, &[], Policy::default()).map_err(|_| MikeError::RngFailure)?;
        let scalar = self.generate_scalar(&mut rng)?;
        Ok(self.keygen_impl(&scalar))
    }
}

/// Structure which holds what becomes the kernel for each of the three
/// gluing isogenies
pub struct KernelPoints<Fp2: Fp2Trait> {
    pub scholten_kernel: [Point<Fp2>; 2],
    pub diagonal_kernel: [Point<Fp2>; 2],
    pub gluing_kernel: [[Point<Fp2>; 2]; 3],
}

impl<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>, const N: usize> MikePrivateKey<Fp, Fp2, N> {
    /// Generate three kernels for the gluing isogeny: scholten kernel, diagonal kernel
    /// and the final gluing kernel, together with the chain kernel which generates the
    /// rest of the dim 4 isogeny chain.
    pub fn generate_kernel_elements(
        &self,
        curve: &Curve<Fp2>,
        a_div_three: &Fp2,
    ) -> (KernelPoints<Fp2>, [[Point<Fp2>; 2]; 2]) {
        // The user's secret key is used for a variety of scalar multiplications
        let s = self.secret;

        // Compute the x-only torsion basis of order 2^(e + 2)
        let e_diff: usize = self.constants.f - (self.constants.e + 2);
        let cofactor_bytes = self.constants.c.to_le_bytes();
        let cofactor_bits = (usize::BITS - self.constants.c.leading_zeros()) as usize;
        let basis = curve.torsion_basis_2e_vartime(e_diff, &cofactor_bytes, cofactor_bits);

        // Full points of order 2^(e + 2) with (X : Y : Z) coordinates
        let (P, Q) = curve.lift_to_jacobian(&basis);

        // Implementation Note:
        // We really want to perform scalar multiplication by (s - 1) / 2,
        // but instead of dividing by 2, we just ignore the first bit of s
        // when we compute the scalar multiplications below.
        // This is done by setting n_start = 1 in the `mul` methods

        // Writing n = (s - 1) / 2
        // Construct the kernel basis equal to (([(n]Q + Q, [n]Q), (-[n]P - P, -[n]P))
        let s_P = curve.mul(&P, &s, 1, self.constants.e);
        let s_Q = curve.mul(&Q, &s, 1, self.constants.e);
        let T1 = curve.add(&s_Q, &Q);
        let T3 = curve.add(&s_P, &P);
        let chain_kernel = [[T1, s_Q], [-T3, -s_P]];

        // Compute intermediate points of order 32 by doubling many times
        let P_32 = curve.double_iter_with_ws(&P, a_div_three, self.constants.e - 3);
        let Q_32 = curve.double_iter_with_ws(&Q, a_div_three, self.constants.e - 3);

        // Compute [n]P_32, [n]Q_32, for the final gluing kernel, we only need the bottom
        // 5 bits of the scalars because of the order of P, Q
        let s_P_32 = curve.mul(&P_32, &s, 1, 6);
        let s_Q_32 = curve.mul(&Q_32, &s, 1, 6);
        let T1 = curve.add(&s_Q_32, &Q_32);
        let T3 = curve.add(&s_P_32, &P_32);
        let T4 = curve.sub(&T1, &T3);
        let T5 = curve.sub(&s_Q_32, &s_P_32);
        let gluing_kernel = [[T1, s_Q_32], [-T3, -s_P_32], [T4, T5]];

        // Compute points of order 16 for the diagonal isogeny kernel
        let P_16 = curve.double(&P_32);
        let Q_16 = curve.double(&Q_32);
        let diagonal_kernel = [Q_16, P_16];

        // Compute points of order 8 for the Scholten isogeny kernel
        let P_8 = curve.double(&P_16);
        let Q_8 = curve.double(&Q_16);
        let scholten_kernel = [P_8, Q_8];

        // Package all data together for the three gluing kernels and the chain kernel
        let kernel_data = KernelPoints {
            scholten_kernel,
            diagonal_kernel,
            gluing_kernel,
        };

        (kernel_data, chain_kernel)
    }

    #[cfg(not(feature = "allow_non_canonical_pk"))]
    fn is_canonical_pk(domain: &Curve<Fp2>, p4: &Point<Fp2>) -> Result<(), MikeError>
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let E_normalized = Curve::normalize_curve(&p4.to_point_x());
        let is_canonical = domain.A.equals(&E_normalized.A) == u32::MAX;
        if !is_canonical {
            return Err(MikeError::NonCanonicalPublicKey);
        }

        Ok(())
    }

    /// A public key is considered supersingular providing that the points we sample
    /// with the torsion basis generation has exactly the order we expect.
    fn is_canonical_supersingular(
        domain: &Curve<Fp2>,
        p8_q8: &[Point<Fp2>],
    ) -> Result<(), MikeError>
    where
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let p4 = domain.double(&p8_q8[0]);
        let p2 = domain.double(&p4);
        let q2 = domain.double_iter(&p8_q8[1], 2);

        // P8 and Q8 have order exactly 8 if P2 and Q2 have y = 0, z != 0
        let p2_order_two = p2.Y.is_zero() & !p2.is_zero();
        let q2_order_two = q2.Y.is_zero() & !q2.is_zero();

        // Additionally, our gluing expects <P, Q> to be a basis with Q2 = (0 : 0 : 1)
        let q2_x_is_zero = q2.X.is_zero();

        let is_supersingular = p2_order_two & q2_order_two & q2_x_is_zero == u32::MAX;
        if !is_supersingular {
            return Err(MikeError::NonSupersingularPublicKey);
        }

        #[cfg(not(feature = "allow_non_canonical_pk"))]
        Self::is_canonical_pk(domain, &p4)?;

        Ok(())
    }

    /// We take the four absolute invariants computed from the codomain
    /// of the isogeny chain and hash them using Sha3_256 to compute a
    /// 32 byte secret used as a shared key
    fn hash_modular_invariants(&self, invariants: &[Fp]) -> [u8; 32]
    where
        [(); Fp::ENCODED_LENGTH]: Sized,
    {
        let mut hasher = Sha3_256::new();
        for inv in invariants.iter() {
            let encoded = inv.encode();
            hasher.update(encoded);
        }
        hasher.finalize().into()
    }

    /// Compute a shared secret between two parties
    pub fn shared_secret(&self, other_pk: &MikePublicKey<Fp2>) -> Result<[u8; 32], MikeError>
    where
        [(); Fp::ENCODED_LENGTH]: Sized,
        [(); Fp2::ENCODED_LENGTH]: Sized,
    {
        let domain = other_pk.curve;
        // TODO: we could precompute 1/3 for each field to avoid this inversion
        let a_div_three = domain.A / Fp2::THREE;

        // Compute the various bases needed for the isogeny chain kernel and the three
        // gluing isogenies
        let (gluing_kernel, chain_kernel) = self.generate_kernel_elements(&domain, &a_div_three);

        // Before we compute the chain, we check whether P8 and Q8 are points
        // of order exactly 8. This depends only on public data, so needs not
        // be constant time. We also ensure the public key is the correct
        // normalization as expected from the key generation.
        Self::is_canonical_supersingular(&domain, &gluing_kernel.scholten_kernel)?;

        // Compute the isogeny codomain
        let codomain = domain.mike_isogeny_chain(
            &a_div_three,
            &gluing_kernel,
            &chain_kernel,
            self.constants.e,
        );

        // Compute the absolute modular invariants
        let invariants = codomain.absolute_modular_invariants();

        // Compute the secret as the hash of these values
        Ok(self.hash_modular_invariants(&invariants))
    }
}

#[cfg(test)]
mod known_answer_test {
    use crate::mike::{MIKE_I, MIKE_III, MIKE_V};
    use fp2::traits::Fq as _;

    macro_rules! mike_kat {
        ($mod_name:ident, $instance:expr, $test_data:path, $seed_alice:expr, $seed_bob:expr) => {
            mod $mod_name {
                use super::*;
                use $test_data::{A_ALICE, A_BOB, SECRET, SK_ALICE, SK_BOB};

                const SEED_ALICE: [u8; 48] = $seed_alice;
                const SEED_BOB: [u8; 48] = $seed_bob;

                #[test]
                fn test_seeded_keygen_alice() {
                    let (pk, sk) = $instance.keygen_seeded(&SEED_ALICE).unwrap();
                    assert_eq!(pk.encode(), A_ALICE.encode());
                    assert_eq!(SK_ALICE, sk.encode());
                }

                #[test]
                fn test_seeded_keygen_bob() {
                    let (pk, sk) = $instance.keygen_seeded(&SEED_BOB).unwrap();
                    assert_eq!(pk.encode(), A_BOB.encode());
                    assert_eq!(SK_BOB, sk.encode());
                }

                #[test]
                fn test_seeded_shared_secret() {
                    let (pk_a, sk_a) = $instance.keygen_seeded(&SEED_ALICE).unwrap();
                    let (pk_b, sk_b) = $instance.keygen_seeded(&SEED_BOB).unwrap();

                    let s1 = sk_a.shared_secret(&pk_b).unwrap();
                    let s2 = sk_b.shared_secret(&pk_a).unwrap();
                    assert_eq!(s1, s2);
                    assert_eq!(s1, SECRET);
                }
            }
        };
    }

    mike_kat!(
        mike_i,
        MIKE_I,
        crate::mike::test_data::mike_i_test_data,
        [0u8; 48],
        [1u8; 48]
    );
    mike_kat!(
        mike_iii,
        MIKE_III,
        crate::mike::test_data::mike_iii_test_data,
        [0u8; 48],
        [1u8; 48]
    );
    mike_kat!(
        mike_v,
        MIKE_V,
        crate::mike::test_data::mike_v_test_data,
        [0u8; 48],
        [1u8; 48]
    );
}

#[cfg(test)]
mod test {
    use super::{MikeError, MikePublicKey};
    use crate::fields::MikeIFp2 as Fp2;
    use crate::mike::MIKE_I;
    use crate::mike::test_data::mike_i_test_data::A_ALICE;
    use fp2::traits::Fq as _;
    use rand_core::OsRng;

    #[test]
    fn test_decode_non_canonical_representation() {
        let bytes = [u8::MAX; Fp2::ENCODED_LENGTH];
        assert!(matches!(
            MikePublicKey::<Fp2>::decode(&bytes),
            Err(MikeError::InvalidPublicKeyEncoding)
        ));
    }

    #[test]
    fn test_pk_in_fp() {
        let bytes = Fp2::ONE.encode();
        assert!(matches!(
            MikePublicKey::<Fp2>::decode(&bytes),
            Err(MikeError::UnsupportedPublicKey)
        ));
    }

    #[test]
    fn test_decode_valid() {
        let bytes = A_ALICE.encode();
        let pk = MikePublicKey::<Fp2>::decode(&bytes);
        assert!(pk.is_ok());
        assert_eq!(pk.unwrap().encode(), bytes);
    }

    #[test]
    fn test_shared_secret() {
        for _ in 0..10 {
            let (pk_a, sk_a) = MIKE_I.keygen(&mut OsRng).unwrap();
            let (pk_b, sk_b) = MIKE_I.keygen(&mut OsRng).unwrap();

            // Ensure the keys are distinct
            assert_ne!(pk_a.encode(), pk_b.encode());
            assert_ne!(sk_a.encode(), sk_b.encode());

            // Ensure the shared secrets match
            let s1 = sk_a.shared_secret(&pk_b).unwrap();
            let s2 = sk_b.shared_secret(&pk_a).unwrap();
            assert_eq!(s1, s2);
        }
    }

    #[test]
    fn test_shared_secret_invalid_public_key() {
        let (_, sk) = MIKE_I.keygen(&mut OsRng).unwrap();
        let bad_a = Fp2::ZETA; // This is an ordinary curve for level one prime
        let bad_pk = MikePublicKey::<Fp2>::decode(&bad_a.encode()).unwrap();
        assert!(matches!(
            sk.shared_secret(&bad_pk),
            Err(MikeError::NonSupersingularPublicKey)
        ));
    }

    #[test]
    #[cfg(not(feature = "allow_non_canonical_pk"))]
    fn test_disallow_non_canonical_public_key() {
        let (pk_a, sk) = MIKE_I.keygen(&mut OsRng).unwrap();
        let non_canonical_a = -pk_a.curve.A; // This is a valid isomorphism, but not canonical
        let pk = MikePublicKey::<Fp2>::decode(&non_canonical_a.encode()).unwrap();
        assert!(matches!(
            sk.shared_secret(&pk),
            Err(MikeError::NonCanonicalPublicKey)
        ));
    }

    #[test]
    #[cfg(feature = "allow_non_canonical_pk")]
    fn test_allow_non_canonical_public_key() {
        let (pk_a, sk_a) = MIKE_I.keygen(&mut OsRng).unwrap();
        let (pk_b, sk_b) = MIKE_I.keygen(&mut OsRng).unwrap();

        let non_canonical_a = -pk_a.curve.A; // This is a valid isomorphism, but not canonical
        let non_canonical_pk = MikePublicKey::<Fp2>::decode(&non_canonical_a.encode()).unwrap();

        // Ensure the shared secrets match, non-canonical produces the same
        let s1 = sk_a.shared_secret(&pk_b).unwrap();
        let s2 = sk_b.shared_secret(&pk_a).unwrap();
        let s3 = sk_b.shared_secret(&non_canonical_pk).unwrap();
        assert_eq!(s1, s2);
        assert_eq!(s1, s3);
    }
}
