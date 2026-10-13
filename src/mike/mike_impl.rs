use crate::elliptic::curve::Curve;
use crate::elliptic::{basis::BasisX, jac_point::Point};
use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;
use nist_drbg_rs::{AesCtr256Drbg, Policy};
use rand_core::{TryCryptoRng, TryRngCore};
use sha3::{Digest, Sha3_256};

/// Constants associated to the field used in Mike
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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
    MismatchedParameters,
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
    pub fn encode(&self) -> Fq::Encoded {
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
    ///
    /// Inputs which are not exactly `Fq::ENCODED_LENGTH` bytes long are rejected.
    pub fn decode(bytes: &[u8]) -> Result<Self, MikeError> {
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

impl<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>, const N: usize> Mike<Fp, Fp2, N> {
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
    fn public_key(&self, scalar: &[u8]) -> MikePublicKey<Fp2> {
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
    ) -> (MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>) {
        let public_key = self.public_key(scalar_bytes);
        let secret_key = self.secret_key(scalar_bytes);

        (public_key, secret_key)
    }

    /// Generate a random keypair using user supplied entropy source
    pub fn keygen<R: TryRngCore + TryCryptoRng>(
        &self,
        rng: &mut R,
    ) -> Result<(MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>), MikeError> {
        let scalar = self.generate_scalar(rng)?;
        Ok(self.keygen_impl(&scalar))
    }

    /// Generate a keypair using a AES-CTR-256 DRBG seeded with the provided value.
    /// Used for deterministic generation / KAT testing.
    pub fn keygen_seeded(
        &self,
        seed: &[u8],
    ) -> Result<(MikePublicKey<Fp2>, MikePrivateKey<Fp, Fp2, N>), MikeError> {
        let mut rng =
            AesCtr256Drbg::new(seed, &[], Policy::default()).map_err(|_| MikeError::RngFailure)?;
        let scalar = self.generate_scalar(&mut rng)?;
        Ok(self.keygen_impl(&scalar))
    }

    /// Prepare and validate a peer public key once for agreement with multiple keys.
    ///
    /// This performs the same supersingularity and canonicality checks as
    /// `MikePrivateKey::shared_secret`, before any private-key computation.
    /// As with scalar agreement, the peer must first be parsed with
    /// `MikePublicKey::decode`; unchecked, malformed curves constructed with
    /// `MikePublicKey::new` are outside this method's input contract.
    pub fn prepare_peer(
        &self,
        other_pk: &MikePublicKey<Fp2>,
    ) -> Result<PreparedPeer<Fp, Fp2, N>, MikeError> {
        let domain = other_pk.curve;
        let constants = self.params.constants;

        let e_diff = constants.f - (constants.e + 2);
        let cofactor_bytes = constants.c.to_le_bytes();
        let cofactor_bits = (usize::BITS - constants.c.leading_zeros()) as usize;
        let basis = domain.torsion_basis_2e_vartime(e_diff, &cofactor_bytes, cofactor_bits);

        // Full points of order 2^(e + 2) with (X : Y : Z) coordinates
        let (P, Q) = domain.lift_to_jacobian(&basis);

        let points = self.prepare_lifted_peer_points(domain, P, Q)?;

        // Before we compute the chain, we check whether P8 and Q8 are points
        // of order exactly 8. This depends only on public data, so needs not
        // be constant time. We also ensure the public key is the correct
        // normalization as expected from the key generation.
        MikePrivateKey::<Fp, Fp2, N>::is_canonical_supersingular(&domain, &points.scholten_kernel)?;

        // The ladder uses the Jacobian source x-coordinate X/Z^2.
        let P_x = points.P.X / points.P.Z.square();
        let Q_x = points.Q.X / points.Q.Z.square();

        Ok(self.finish_prepared_peer(points, P_x, Q_x))
    }

    // Complete the public-point work after the (possibly batch-normalized) basis is lifted.
    fn prepare_lifted_peer_points(
        &self,
        curve: Curve<Fp2>,
        P: Point<Fp2>,
        Q: Point<Fp2>,
    ) -> Result<PeerPoints<Fp2>, MikeError> {
        let constants = self.params.constants;
        let a_div_three = curve.A / Fp2::THREE;

        let (chain_doublings, [P_32, Q_32]) =
            curve.prepare_chain_doublings(&a_div_three, &P, &Q, constants.e);

        // Compute points of order 16 for the diagonal isogeny kernel
        let P_16 = curve.double(&P_32);
        let Q_16 = curve.double(&Q_32);

        // Compute points of order 8 for the Scholten isogeny kernel
        let P_8 = curve.double(&P_16);
        let Q_8 = curve.double(&Q_16);

        // The batch path defers canonicalization until the coefficients of all
        // peers can share one inversion batch. Scalar preparation checks here.
        Ok(PeerPoints {
            domain: curve,
            a_div_three,
            chain_doublings,
            P,
            Q,
            P_32,
            Q_32,
            diagonal_kernel: [Q_16, P_16],
            scholten_kernel: [P_8, Q_8],
        })
    }

    fn finish_prepared_peer(
        &self,
        points: PeerPoints<Fp2>,
        P_x: Fp2,
        Q_x: Fp2,
    ) -> PreparedPeer<Fp, Fp2, N> {
        let PeerPoints {
            domain,
            a_div_three,
            chain_doublings,
            P,
            Q,
            P_32,
            Q_32,
            diagonal_kernel,
            scholten_kernel,
        } = points;
        PreparedPeer {
            _phantom_fp: core::marker::PhantomData,
            _phantom_n: core::marker::PhantomData,
            domain,
            a_div_three,
            chain_doublings,
            P,
            Q,
            P_x,
            Q_x,
            P_32,
            Q_32,
            diagonal_kernel,
            scholten_kernel,
            constants: self.params.constants,
        }
    }

    /// Prepare multiple peers, sharing the first lifted point's affine inversion.
    ///
    /// Each result remains at its input position. This retains the same validation
    /// checks as `prepare_peer`; variable-time torsion search still runs per peer.
    pub fn batch_prepare_peers(
        &self,
        peers: &[MikePublicKey<Fp2>],
    ) -> Vec<Result<PreparedPeer<Fp, Fp2, N>, MikeError>> {
        let constants = self.params.constants;
        let e_diff = constants.f - (constants.e + 2);
        let cofactor_bytes = constants.c.to_le_bytes();
        let cofactor_bits = (usize::BITS - constants.c.leading_zeros()) as usize;

        let bases: Vec<_> = peers
            .iter()
            .map(|peer| {
                peer.curve
                    .torsion_basis_2e_vartime(e_diff, &cofactor_bytes, cofactor_bits)
            })
            .collect();
        let mut p_z_inv: Vec<_> = bases.iter().map(|basis| basis.P.Z).collect();
        Fp2::batch_invert(&mut p_z_inv);
        let stages: Vec<_> = peers
            .iter()
            .zip(bases)
            .zip(p_z_inv)
            .map(|((peer, basis), p_z_inv)| {
                let domain = peer.curve;
                let p_x = basis.P.X * p_z_inv;
                let (P, Q) = domain.lift_to_jacobian_with_affine_p(&basis, p_x);
                self.prepare_lifted_peer_points(domain, P, Q)
            })
            .collect();

        // Validate point orders before canonicalization. Keep failed public
        // lanes in the same positions for subsequent inversion batches.
        let checked: Vec<_> = stages
            .into_iter()
            .map(|stage| {
                stage.and_then(|points| {
                    let p4 = MikePrivateKey::<Fp, Fp2, N>::check_supersingular_order(
                        &points.domain,
                        &points.scholten_kernel,
                    )?;
                    Ok((points, p4))
                })
            })
            .collect();
        #[cfg(not(feature = "allow_non_canonical_pk"))]
        let canonical_checks: Vec<_> = {
            let coefficients: Vec<_> = checked
                .iter()
                .map(|entry| match entry {
                    Ok((_, p4)) => Curve::normalization_coefficients(&p4.to_point_x()),
                    Err(_) => ([Fp2::ZERO; 3], [Fp2::ONE; 3]),
                })
                .collect();
            let mut inverses: Vec<_> = coefficients.iter().flat_map(|(_, c)| *c).collect();
            Fp2::batch_invert(&mut inverses);
            coefficients
                .into_iter()
                .zip(inverses.chunks_exact(3))
                .map(|((a, _), inv)| {
                    Curve::normalize_curve_with_inverses(a, [inv[0], inv[1], inv[2]])
                })
                .collect()
        };

        let stages: Vec<_> = checked
            .into_iter()
            .enumerate()
            .map(|(index, entry)| {
                #[cfg(not(feature = "allow_non_canonical_pk"))]
                {
                    entry.and_then(|(points, _p4)| {
                        if points.domain.A.equals(&canonical_checks[index].A) != u32::MAX {
                            return Err(MikeError::NonCanonicalPublicKey);
                        }
                        Ok(points)
                    })
                }
                #[cfg(feature = "allow_non_canonical_pk")]
                {
                    let _ = index;
                    entry.map(|(points, _p4)| points)
                }
            })
            .collect();

        // The ladder uses Jacobian source coordinates X/Z^2. Include failing
        // public lanes to keep their original positions and zero-inversion behavior.
        let mut z_inv: Vec<_> = stages
            .iter()
            .map(|stage| match stage {
                Ok(points) => points.P.Z.square(),
                Err(_) => Fp2::ONE,
            })
            .chain(stages.iter().map(|stage| match stage {
                Ok(points) => points.Q.Z.square(),
                Err(_) => Fp2::ONE,
            }))
            .collect();
        Fp2::batch_invert(&mut z_inv);
        let (p_z_inv, q_z_inv) = z_inv.split_at(peers.len());
        stages
            .into_iter()
            .enumerate()
            .map(|(index, stage)| {
                stage.map(|points| {
                    let P_x = points.P.X * p_z_inv[index];
                    let Q_x = points.Q.X * q_z_inv[index];
                    self.finish_prepared_peer(points, P_x, Q_x)
                })
            })
            .collect()
    }

    /// Agree with prepared peers, batching final invariant inversions across trials.
    ///
    /// Returns one result per peer/key pair in peer-major order. Prepare each
    /// peer with `prepare_peer` before calling this method; a failed preparation
    /// must be handled by the caller at its original input position. Errors for
    /// mismatched parameters occupy their original positions. All valid lanes
    /// reach the same inversion batch, including lanes with zero theta_0;
    /// no secret-dependent compaction is performed.
    pub fn batch_shared_secrets(
        &self,
        receiver_keys: &[MikePrivateKey<Fp, Fp2, N>],
        peers: &[PreparedPeer<Fp, Fp2, N>],
    ) -> Vec<Result<[u8; 32], MikeError>> {
        if receiver_keys.is_empty() || peers.is_empty() {
            return Vec::new();
        }

        let trials = peers.len() * receiver_keys.len();
        let mut theta_0_values = Vec::with_capacity(trials);
        let mut theta_n_values = Vec::with_capacity(trials);
        let mut results = vec![Err(MikeError::MismatchedParameters); trials];

        for peer in peers {
            for key in receiver_keys {
                if key.constants != peer.constants {
                    // Keep an inert lane so the batch shape is public and fixed.
                    theta_0_values.push(Fp::ONE);
                    theta_n_values.push([Fp::ZERO; 9]);
                    continue;
                }

                let (gluing_kernel, chain_kernel) = key.generate_kernels_from_prepared(
                    &peer.domain,
                    &peer.P,
                    Some(&peer.P_x),
                    &peer.Q,
                    Some(&peer.Q_x),
                    &peer.P_32,
                    &peer.Q_32,
                    peer.scholten_kernel,
                    peer.diagonal_kernel,
                );

                let codomain = peer.domain.mike_isogeny_chain_with_prepared_doublings(
                    &peer.a_div_three,
                    &gluing_kernel,
                    &chain_kernel,
                    key.constants.e,
                    Some(&peer.chain_doublings),
                );

                let (theta_0, theta_n) = codomain.special_modular_invariants_parts();
                theta_0_values.push(theta_0);
                theta_n_values.push(theta_n);
            }
        }

        Fp::batch_invert(&mut theta_0_values);

        for (position, (theta_0_inv, theta_n)) in
            theta_0_values.iter().zip(theta_n_values).enumerate()
        {
            let special =
                crate::theta::dim_four::ThetaStructure::<Fp>::normalize_special_modular_invariants(
                    theta_n,
                    theta_0_inv,
                );
            let invariants = crate::theta::dim_four::ThetaStructure::<Fp>::absolute_modular_invariants_from_special(special);
            let key_index = position % receiver_keys.len();
            let secret = receiver_keys[key_index].hash_modular_invariants(&invariants);

            if receiver_keys[key_index].constants == peers[position / receiver_keys.len()].constants
            {
                results[position] = Ok(secret);
            }
        }

        results
    }
}

// Public intermediate points shared by scalar and batched preparation.
struct PeerPoints<Fp2: Fp2Trait> {
    domain: Curve<Fp2>,
    a_div_three: Fp2,
    chain_doublings: super::isogeny_chain::PreparedDoublings<Fp2>,
    P: Point<Fp2>,
    Q: Point<Fp2>,
    P_32: Point<Fp2>,
    Q_32: Point<Fp2>,
    diagonal_kernel: [Point<Fp2>; 2],
    scholten_kernel: [Point<Fp2>; 2],
}

/// A validated peer public key with its public torsion points prepared for reuse.
///
/// Construct via [`Mike::prepare_peer`]. The parameter set is bound to this value;
/// use [`Mike::batch_shared_secrets`] to process multiple peers and receiver keys.
/// Preparation is independent of the receiver's private scalar. It retains public
/// basis multiples used by the isogeny strategy, so cloning also clones that cache.
///
/// This value contains no secret data, but its construction has variable runtime
/// dependent on the public key.
#[derive(Clone, Debug)]
pub struct PreparedPeer<Fp: FqTrait, Fp2: Fp2Trait, const N: usize> {
    _phantom_fp: core::marker::PhantomData<Fp>,
    _phantom_n: core::marker::PhantomData<[u8; N]>,
    domain: Curve<Fp2>,
    a_div_three: Fp2,
    chain_doublings: super::isogeny_chain::PreparedDoublings<Fp2>,
    P: Point<Fp2>,
    Q: Point<Fp2>,
    P_x: Fp2,
    Q_x: Fp2,
    P_32: Point<Fp2>,
    Q_32: Point<Fp2>,
    diagonal_kernel: [Point<Fp2>; 2],
    scholten_kernel: [Point<Fp2>; 2],
    constants: MikeConstants,
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
        // Compute the x-only torsion basis of order 2^(e + 2)
        let e_diff: usize = self.constants.f - (self.constants.e + 2);
        let cofactor_bytes = self.constants.c.to_le_bytes();
        let cofactor_bits = (usize::BITS - self.constants.c.leading_zeros()) as usize;
        let basis = curve.torsion_basis_2e_vartime(e_diff, &cofactor_bytes, cofactor_bits);

        // Full points of order 2^(e + 2) with (X : Y : Z) coordinates
        let (P, Q) = curve.lift_to_jacobian(&basis);

        // Compute intermediate points of order 32 by doubling many times
        let P_32 = curve.double_iter_with_ws(&P, a_div_three, self.constants.e - 3);
        let Q_32 = curve.double_iter_with_ws(&Q, a_div_three, self.constants.e - 3);

        // Compute points of order 16 for the diagonal isogeny kernel
        let P_16 = curve.double(&P_32);
        let Q_16 = curve.double(&Q_32);
        let diagonal_kernel = [Q_16, P_16];

        // Compute points of order 8 for the Scholten isogeny kernel
        let P_8 = curve.double(&P_16);
        let Q_8 = curve.double(&Q_16);
        let scholten_kernel = [P_8, Q_8];

        self.generate_kernels_from_prepared(
            curve,
            &P,
            None,
            &Q,
            None,
            &P_32,
            &Q_32,
            scholten_kernel,
            diagonal_kernel,
        )
    }

    /// Compute the private-key-dependent kernels from prepared public data.
    fn generate_kernels_from_prepared(
        &self,
        curve: &Curve<Fp2>,
        P: &Point<Fp2>,
        P_x: Option<&Fp2>,
        Q: &Point<Fp2>,
        Q_x: Option<&Fp2>,
        P_32: &Point<Fp2>,
        Q_32: &Point<Fp2>,
        scholten_kernel: [Point<Fp2>; 2],
        diagonal_kernel: [Point<Fp2>; 2],
    ) -> (KernelPoints<Fp2>, [[Point<Fp2>; 2]; 2]) {
        // The user's secret key is used for a variety of scalar multiplications
        let s = &self.secret;

        // Implementation Note:
        // We really want to perform scalar multiplication by (s - 1) / 2,
        // but instead of dividing by 2, we just ignore the first bit of s
        // when we compute the scalar multiplications below.
        // This is done by setting n_start = 1 in the `mul` methods

        // Writing n = (s - 1) / 2
        // Construct the kernel basis equal to (([(n]Q + Q, [n]Q), (-[n]P - P, -[n]P))
        let s_P = curve.mul_with_affine_x(P, P_x, s, 1, self.constants.e);
        let s_Q = curve.mul_with_affine_x(Q, Q_x, s, 1, self.constants.e);
        let T1 = curve.add(&s_Q, Q);
        let T3 = curve.add(&s_P, P);
        let chain_kernel = [[T1, s_Q], [-T3, -s_P]];

        // Compute [n]P_32, [n]Q_32, for the final gluing kernel, we only need the bottom
        // 5 bits of the scalars because of the order of P, Q
        let s_P_32 = curve.mul(P_32, s, 1, 6);
        let s_Q_32 = curve.mul(Q_32, s, 1, 6);
        let T1 = curve.add(&s_Q_32, Q_32);
        let T3 = curve.add(&s_P_32, P_32);
        let T4 = curve.sub(&T1, &T3);
        let T5 = curve.sub(&s_Q_32, &s_P_32);
        let gluing_kernel = [[T1, s_Q_32], [-T3, -s_P_32], [T4, T5]];

        // Package all data together for the three gluing kernels and the chain kernel
        let kernel_data = KernelPoints {
            scholten_kernel,
            diagonal_kernel,
            gluing_kernel,
        };

        (kernel_data, chain_kernel)
    }

    /// Agree with a prepared, validated peer public key.
    ///
    /// Reuse a prepared peer only with receiver keys of the same parameter set.
    /// This returns an error for a mismatched parameter set instead of using the
    /// wrong isogeny length. Validation errors are returned by `prepare_peer`.
    pub fn shared_secret_from_prepared(
        &self,
        peer: &PreparedPeer<Fp, Fp2, N>,
    ) -> Result<[u8; 32], MikeError> {
        if self.constants != peer.constants {
            return Err(MikeError::MismatchedParameters);
        }

        let (gluing_kernel, chain_kernel) = self.generate_kernels_from_prepared(
            &peer.domain,
            &peer.P,
            Some(&peer.P_x),
            &peer.Q,
            Some(&peer.Q_x),
            &peer.P_32,
            &peer.Q_32,
            peer.scholten_kernel,
            peer.diagonal_kernel,
        );

        let codomain = peer.domain.mike_isogeny_chain_with_prepared_doublings(
            &peer.a_div_three,
            &gluing_kernel,
            &chain_kernel,
            self.constants.e,
            Some(&peer.chain_doublings),
        );

        Ok(self.hash_modular_invariants(&codomain.absolute_modular_invariants()))
    }

    /// A public key is considered supersingular providing that the points we sample
    /// with the torsion basis generation has exactly the order we expect.
    fn is_canonical_supersingular(
        domain: &Curve<Fp2>,
        p8_q8: &[Point<Fp2>],
    ) -> Result<(), MikeError> {
        let p4 = Self::check_supersingular_order(domain, p8_q8)?;

        #[cfg(not(feature = "allow_non_canonical_pk"))]
        Self::is_canonical_pk(domain, &p4)?;
        #[cfg(feature = "allow_non_canonical_pk")]
        let _ = p4;

        Ok(())
    }

    fn check_supersingular_order(
        domain: &Curve<Fp2>,
        p8_q8: &[Point<Fp2>],
    ) -> Result<Point<Fp2>, MikeError> {
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
        Ok(p4)
    }

    #[cfg(not(feature = "allow_non_canonical_pk"))]
    fn is_canonical_pk(domain: &Curve<Fp2>, p4: &Point<Fp2>) -> Result<(), MikeError> {
        let E_normalized = Curve::normalize_curve(&p4.to_point_x());
        let is_canonical = domain.A.equals(&E_normalized.A) == u32::MAX;
        if !is_canonical {
            return Err(MikeError::NonCanonicalPublicKey);
        }

        Ok(())
    }

    /// We take the four absolute invariants computed from the codomain
    /// of the isogeny chain and hash them using Sha3_256 to compute a
    /// 32 byte secret used as a shared key
    fn hash_modular_invariants(&self, invariants: &[Fp]) -> [u8; 32] {
        let mut hasher = Sha3_256::new();
        for inv in invariants.iter() {
            let encoded = inv.encode();
            hasher.update(encoded);
        }
        hasher.finalize().into()
    }

    /// Compute a shared secret between two parties
    pub fn shared_secret(&self, other_pk: &MikePublicKey<Fp2>) -> Result<[u8; 32], MikeError> {
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
mod prepared_test {
    use super::*;

    fn check<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>, const N: usize>(mike: Mike<Fp, Fp2, N>) {
        for seed in [7u8, 19] {
            let (pk, _) = mike.keygen_seeded(&[seed; 48]).unwrap();
            let peer = mike.prepare_peer(&pk).unwrap();
            assert_eq!(
                peer.P_32.equals(&peer.domain.double_iter_with_ws(
                    &peer.P,
                    &peer.a_div_three,
                    peer.constants.e - 3
                )),
                u32::MAX
            );
            assert_eq!(
                peer.Q_32.equals(&peer.domain.double_iter_with_ws(
                    &peer.Q,
                    &peer.a_div_three,
                    peer.constants.e - 3
                )),
                u32::MAX
            );
            for index in 0..4 {
                let (_, key) = mike.keygen_seeded(&[index; 48]).unwrap();
                let (gluing, kernel) = key.generate_kernels_from_prepared(
                    &peer.domain,
                    &peer.P,
                    Some(&peer.P_x),
                    &peer.Q,
                    Some(&peer.Q_x),
                    &peer.P_32,
                    &peer.Q_32,
                    peer.scholten_kernel,
                    peer.diagonal_kernel,
                );
                let expected = peer.domain.mike_isogeny_chain(
                    &peer.a_div_three,
                    &gluing,
                    &kernel,
                    key.constants.e,
                );
                let actual = peer.domain.mike_isogeny_chain_with_prepared_doublings(
                    &peer.a_div_three,
                    &gluing,
                    &kernel,
                    key.constants.e,
                    Some(&peer.chain_doublings),
                );
                assert_eq!(actual.null_point().equals(&expected.null_point()), u32::MAX);
                assert_eq!(
                    key.shared_secret_from_prepared(&peer).unwrap(),
                    key.shared_secret(&pk).unwrap()
                );
            }
        }
    }

    #[test]
    fn prepared_doublings_match_unprepared() {
        check(crate::mike::MIKE_I);
        check(crate::mike::MIKE_III);
        check(crate::mike::MIKE_V);
    }
}

#[cfg(test)]
mod batch_test {
    use super::{Mike, MikeError, MikePrivateKey, MikePublicKey};
    use crate::mike::{MIKE_I, MIKE_III, MIKE_V};
    use fp2::traits::{Fp2 as Fp2Trait, Fq as FqTrait};

    fn check_batch<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>, const N: usize>(
        mike: Mike<Fp, Fp2, N>,
    ) {
        let (peer_a, key_a) = mike.keygen_seeded(&[0u8; 48]).unwrap();
        let (peer_b, key_b) = mike.keygen_seeded(&[1u8; 48]).unwrap();
        let (peer_c, key_c) = mike.keygen_seeded(&[2u8; 48]).unwrap();
        let keys = [key_a, key_b, key_c];
        let peers = [peer_a, peer_b, peer_c];
        let expected: Vec<_> = peers
            .iter()
            .flat_map(|peer| keys.iter().map(|key| key.shared_secret(peer)))
            .collect();
        for (peer_index, peer) in peers.iter().enumerate() {
            let prepared = mike.prepare_peer(peer).unwrap();
            for (key_index, key) in keys.iter().enumerate() {
                assert_eq!(
                    key.shared_secret_from_prepared(&prepared),
                    expected[peer_index * keys.len() + key_index]
                );
            }
        }
        // An ordinary but well-formed peer is rejected at preparation time.
        let invalid = MikePublicKey::<Fp2>::new(crate::elliptic::curve::Curve::new(&Fp2::ZETA).0);
        let invalid_error = mike.prepare_peer(&invalid).err().unwrap();
        let mixed = [
            MikePublicKey::new(peers[0].curve),
            invalid,
            MikePublicKey::new(peers[2].curve),
        ];
        let mixed_prepared = mike.batch_prepare_peers(&mixed);
        assert_eq!(mixed_prepared.len(), mixed.len());
        assert_eq!(mixed_prepared[1].as_ref().err(), Some(&invalid_error));
        assert!(mixed_prepared[0].is_ok() && mixed_prepared[2].is_ok());

        let mut mismatched = mike.params.constants;
        mismatched.e -= 1;
        let wrong_key = MikePrivateKey::new(keys[0].encode(), mismatched);
        let prepared = mike.prepare_peer(&peers[0]).unwrap();
        assert_eq!(
            wrong_key.shared_secret_from_prepared(&prepared),
            Err(MikeError::MismatchedParameters)
        );
        let prepared_peers: Vec<_> = peers
            .iter()
            .map(|peer| mike.prepare_peer(peer).unwrap())
            .collect();
        let batch_prepared: Vec<_> = mike
            .batch_prepare_peers(&peers)
            .into_iter()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(mike.batch_shared_secrets(&keys, &batch_prepared), expected);
        for (scalar, batched) in prepared_peers.iter().zip(&batch_prepared) {
            for key in &keys {
                assert_eq!(
                    key.shared_secret_from_prepared(scalar),
                    key.shared_secret_from_prepared(batched)
                );
            }
        }
        assert!(mike.batch_prepare_peers(&[]).is_empty());
        assert_eq!(mike.batch_shared_secrets(&keys, &prepared_peers), expected);
        assert!(mike.batch_shared_secrets(&keys, &[]).is_empty());
        assert!(mike.batch_shared_secrets(&[], &prepared_peers).is_empty());
        assert_eq!(
            mike.batch_shared_secrets(&[wrong_key], &prepared_peers)[0],
            Err(MikeError::MismatchedParameters)
        );
    }

    #[test]
    fn all_parameter_sets() {
        check_batch(MIKE_I);
        check_batch(MIKE_III);
        check_batch(MIKE_V);
    }

    #[test]
    fn mismatched_constants_are_rejected_in_batch() {
        let (peer, key) = MIKE_I.keygen_seeded(&[0u8; 48]).unwrap();
        let mut constants = MIKE_I.params.constants;
        constants.e -= 1;
        let wrong_key = MikePrivateKey::new(key.encode(), constants);
        assert_eq!(
            MIKE_I.batch_shared_secrets(&[wrong_key], &[MIKE_I.prepare_peer(&peer).unwrap()]),
            vec![Err(MikeError::MismatchedParameters)]
        );
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
