use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;

use crate::elliptic::curve::Curve;
use crate::elliptic::jac_point::Point;
use crate::mike::mike_impl::KernelPoints;
use crate::theta::dim_two::ThetaPointDimTwo;

pub struct ScholtenImages<Fp: FqTrait> {
    pub diagonal_kernel: [ThetaPointDimTwo<Fp>; 4],
    pub gluing_kernel: [[ThetaPointDimTwo<Fp>; 2]; 3],
    pub image_plus: Vec<[ThetaPointDimTwo<Fp>; 2]>,
    pub image_minus: Vec<[ThetaPointDimTwo<Fp>; 2]>,
}

pub struct ScholtenContext<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>> {
    pub shift: Point<Fp2>,
    pub precomp: [Fp; 5],
    pub J: ThetaPointDimTwo<Fp>,
}

/// Helper function for Fp2 arithmetic to compute the norm
/// n(a) for a = x + i*y = x^2 + y^2
#[inline(always)]
fn fp2_norm<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(a: &Fp2) -> Fp {
    a.x0().square() + a.x1().square()
}

/// Compute the matrix coefficients `mi` for the Scholten change of basis
fn symmetric_action_change_of_basis<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    E: &Curve<Fp2>,
    P8: &Point<Fp2>,
) -> [Fp; 4] {
    // Get point of order four from P8, we use that Q8 is above (0 : 1) to avoid
    // the need to double this, and can use x-only as we only need (X : Z)
    let P4 = E.xdbl(&P8.to_point_x());
    let (alpha, beta) = P4.coords();

    // compute values used to compute `mi`
    let n_alpha = fp2_norm(&alpha); // 2s + 1a on Fp
    let n_beta = fp2_norm(&beta); // 2s + 1a on Fp

    let m0 = n_alpha + n_beta;
    let m2 = n_beta - n_alpha;

    let t0 = (alpha * beta.conjugate()).mul2(); // 1M + 1A on Fp2
    let m1 = t0.x0();
    let m3 = t0.x1();

    [m0, m1, m2, m3]
}

// Given a point P on E, compute the square of the point after computing the
// basis change from (P, sigma(P)) to the theta point. The square is done within
// this function so we can cleanly return a point with coordinates over Fp rather
// than Fp2
fn point_to_squared_theta_point<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    mi: &[Fp; 4],
    P: &Point<Fp2>,
) -> ThetaPointDimTwo<Fp> {
    // This function tries to do a lot at once, but this is done so we can
    // algebraically track the real and imaginary components and work in Fp
    // as much as possible.
    //
    // The idea is that P = (X : Z) has coordinates in Fp2 but the output
    // will necessarily have coordinates in Fp due to the following.
    //
    // First we compute the theta point corresponding to (P, sigma(P))
    // which has coordinates (n(X) : X Z_bar : X_bar Z : n(Z)), with
    // the first and last coordinates in Fp but the second and third
    // in Fp2.
    //
    // We then want to compute the change of basis.
    // First we compute sums and differences to get coordinates
    // (t0, t1, t2, i * t3) where ti are all real
    //
    // The coordinates after the change of basis are then:
    //
    // X = m0 * t0 + m1 * t2;
    // Y = m1 * t0 + m0 * t2;
    // Z = m2 * t1 - m3 * t3;
    // T = i * (m3 * t1 + m2 * t3);
    //
    // Now X, Y, Z are seen to be purely real and T is purely imaginary
    // but as we never need these coordinates, only their squares, we
    // return (X^2 : Y^2 : Z^2 -T'^2) as our output where T' = m2 * t3 + m3 * t1

    let (x, z) = P.to_point_x().coords();
    let [m0, m1, m2, m3] = *mi;

    let n_x = fp2_norm(&x);
    let n_z = fp2_norm(&z);
    let tmp = (x * z.conjugate()).mul2();
    let t0 = n_x + n_z;
    let t1 = n_x - n_z;
    let t2 = tmp.x0();
    let t3 = tmp.x1();

    let mut X = m0 * t0 - m1 * t2;
    let mut Y = m0 * t2 - m1 * t0;
    let mut Z = m2 * t1 - m3 * t3;
    let mut T = m2 * t3 + m3 * t1;

    X.set_square();
    Y.set_square();
    Z.set_square();
    T.set_square();
    T.set_neg();

    ThetaPointDimTwo::new(&X, &Y, &Z, &T)
}

/// Compute the ThetaPoint J, used to compute images under the Scholten isogeny.
/// This is basically a codomain computation, but we don't need to compute the
/// codomain itself, as we use torsion above the kernel to directly compute the
/// following isogenies.
fn scholten_j_precomputation<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    mi: &[Fp; 4],
    P8: &Point<Fp2>,
    Q8: &Point<Fp2>,
) -> ThetaPointDimTwo<Fp> {
    // Move the kernel from E x E^sigma to the theta coordinates in the basis
    // dictated by the symmetric change of basis `mi`
    let mut T1 = point_to_squared_theta_point(mi, P8);
    let mut T2 = point_to_squared_theta_point(mi, Q8);
    T1.set_hadamard();
    T2.set_hadamard();

    // Extract out the coordinates from H(S(Ti)) to compute the codomain
    let (xA, _, yC, yD) = T1.coords();
    let (zA, _, zY, tD) = T2.coords();
    debug_assert_eq!(yD.is_zero(), u32::MAX);
    debug_assert_eq!(tD.is_zero(), u32::MAX);

    let x = xA * zY;
    let y = yC * zA;

    ThetaPointDimTwo::new(&y, &y, &x, &x)
}

/// Compute the five constants A, B, C, D, E required for the computation of images
/// under the Scholten isogeny using the fast "superglue" techniques. Requires the
/// matrix coefficients `mi` computed from the symmetric action of `P8`
fn scholten_constants_precomputation<Fp: FqTrait>(mi: &[Fp; 4]) -> [Fp; 5] {
    let [m0, m1, m2, m3] = *mi;

    // three squares
    let m0_sqr = m0.square();
    let m1_sqr = m1.square();
    let m2_sqr = m2.square();

    let m3_sqr = m3.square();
    debug_assert_eq!((m0_sqr - m3_sqr).equals(&(m1_sqr + m2_sqr)), u32::MAX);

    // three additions
    let a = m1_sqr + m2_sqr;
    let b = m0_sqr - m1_sqr;
    let c = m0_sqr - m2_sqr;

    // two final mul, with scaling
    let d = (m0 * m1).mul4();
    let e = (m2 * m3).mul4();

    [a, b, c, d, e]
}

/// Compute the image of a single elliptic curve point P under the Scholten
/// isogeny, using the precomputed context `ctx`.
fn scholten_image<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    E: &Curve<Fp2>,
    P: &Point<Fp2>,
    ctx: &ScholtenContext<Fp, Fp2>,
    apply_basis_change: bool,
) -> ThetaPointDimTwo<Fp> {
    // Compute two theta points from add components
    // TODO: things could be done over Fp here to
    // hyper optimise things.
    let (u, v, w) = E.add_components(P, &ctx.shift);

    let u_abs = fp2_norm(&u);
    let v_abs = fp2_norm(&v);
    let w_abs = fp2_norm(&w);

    // n(a) = a0^2 + a1^2
    // a^2 = a0^2 + i * 2 * a0 * a1 - a1^2
    // We could potentially make this faster by reusing the above
    // but for now we keep it simple
    let u_sqr = u.square();
    let v_sqr = v.square();
    let w_sqr = w.square();

    // Compute u^2 - v^2 ± w^2
    let t0 = u_sqr - v_sqr;
    let tp = t0 + w_sqr;
    let tm = t0 - w_sqr;

    // Compute n(u^2 - v^2 ± w^2)
    let tp_abs = fp2_norm(&tp);
    let tm_abs = fp2_norm(&tm);

    // Compute U + iV = (u^2 - v^2 + w^2) * bar(u * w)
    let t0 = tp * (u * w).conjugate();
    let U = t0.x0();
    let V = t0.x1();

    let [a, b, c, d, e] = ctx.precomp;
    let t0 = (u_abs * w_abs).mul4();
    let t1 = d * U;
    let t2 = e * V;

    // alpha = A * n(u^2 - v^2 + w^2) + 4C * n(u) * n(w) - DU - EV
    let alpha = a * tp_abs + c * t0 - t1 - t2;

    // gamma = C * n(u^2 - v^2 + w^2) + 4A * n(u) * n(w) - DU + EV
    let gamma = c * tp_abs + a * t0 - t1 + t2;

    // beta = B * n(u^2 - v^2 - w^2)
    let beta = b * tm_abs;

    // beta = (-1)^b * B * n(u^2 - v^2 - w^2)
    // For us, b = 1 for all cases as we only eval (P, -sigma(P))
    let delta = -b * (v_abs * w_abs).mul4();

    let mut image = ThetaPointDimTwo::new(&alpha, &beta, &gamma, &delta);
    image.set_pointwise_mul(&ctx.J);
    image.set_hadamard();

    if apply_basis_change {
        image.set_scholten_change_of_basis();
    }

    image
}

/// Helper function which computes images of points structured into arrays of the form
/// [[P; M]; N], just hides some ugliness of how to apply this map.
fn apply_scholten_isogeny<
    Fp: FqTrait,
    Fp2: Fp2Trait<BaseField = Fp>,
    const N: usize,
    const M: usize,
>(
    E: &Curve<Fp2>,
    basis: &[[Point<Fp2>; M]; N],
    ctx: &ScholtenContext<Fp, Fp2>,
    apply_basis_change: bool,
) -> [[ThetaPointDimTwo<Fp>; M]; N] {
    basis.map(|row| row.map(|point| scholten_image(E, &point, ctx, apply_basis_change)))
}

/// Compute the Ui points directly from the Vi points with 8a
fn U_from_V<Fp: FqTrait>(V: &ThetaPointDimTwo<Fp>) -> ThetaPointDimTwo<Fp> {
    let (x, y, z, t) = V.coords();
    // H . negate t . H
    let (t0, t1) = (x + t, x - t);
    let (t2, t3) = (y + z, y - z);
    let a = t1 + t2;
    let b = t0 + t3;
    let c = t0 - t3;
    let d = t2 - t1;
    ThetaPointDimTwo::new(&a, &b, &c, &d)
}

/// Compute the first gluing Scholten isogeny and push through all kernel data
pub fn scholten_isogeny<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    E: &Curve<Fp2>,
    kernel_points: &KernelPoints<Fp2>,
    image_points: &[[Point<Fp2>; 2]],
) -> ScholtenImages<Fp> {
    // Extract the kernel of the Scholten isogeny
    let [P8, Q8] = kernel_points.scholten_kernel;

    // Compute the matrix coefficents `mi` from `P8` using that `Q8` is above (1 : 0)
    let mi = symmetric_action_change_of_basis(E, &P8);

    // Compute the value J (y : y : x : x) used for computing generic images
    // which requires the data from `mi` together with the kernel `P8` and `Q8`
    let J = scholten_j_precomputation(&mi, &P8, &Q8);

    // Precompute the gluing constants and store P8 as the shift point
    let precomp = scholten_constants_precomputation(&mi);

    // Create a ctx object which holds the data required for all isogeny eval computation
    let ctx: ScholtenContext<Fp, Fp2> = ScholtenContext {
        shift: P8,
        precomp,
        J,
    };

    // Compute the kernel of the diagonal isogeny, we compute Ui directly from Vi without a Scholten image
    let [[V1, V2]] = apply_scholten_isogeny(E, &[kernel_points.diagonal_kernel], &ctx, false);
    let [U1, U2] = [V2, V1].map(|v| U_from_V(&v));

    // Apply change of basis to all elements of the diagonal kernel, we cannot do this before computing Ui without extra cost?
    let mut diagonal_kernel: [ThetaPointDimTwo<Fp>; 4] = [U1, U2, V1, V2];
    for P in diagonal_kernel.iter_mut() {
        P.set_scholten_change_of_basis();
    }

    // Push the final gluing kernel through the Scholten isogeny
    let gluing_kernel: [[ThetaPointDimTwo<Fp>; 2]; 3] =
        apply_scholten_isogeny(E, &kernel_points.gluing_kernel, &ctx, true);

    // Push images through. For each pair of points (P, Q) we compute (P + T, Q + T) and (P - T, Q - T)
    // and then each pair of these points are kept within image_plus and image_minus
    let T_trans = kernel_points.gluing_kernel[0];
    let image_plus: Vec<[ThetaPointDimTwo<Fp>; 2]> = image_points
        .iter()
        .map(|&[P, Q]| {
            [
                scholten_image(E, &E.add(&P, &T_trans[0]), &ctx, true),
                scholten_image(E, &E.add(&Q, &T_trans[1]), &ctx, true),
            ]
        })
        .collect();

    let image_minus: Vec<[ThetaPointDimTwo<Fp>; 2]> = image_points
        .iter()
        .map(|&[P, Q]| {
            [
                scholten_image(E, &E.sub(&P, &T_trans[0]), &ctx, true),
                scholten_image(E, &E.sub(&Q, &T_trans[1]), &ctx, true),
            ]
        })
        .collect();

    ScholtenImages {
        diagonal_kernel,
        gluing_kernel,
        image_plus,
        image_minus,
    }
}
