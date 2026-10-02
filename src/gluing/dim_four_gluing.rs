use fp2::traits::Fq as FqTrait;

use crate::gluing::diagonal_gluing::DiagonalImages;
use crate::theta::dim_four::{ThetaPoint, ThetaPointCompressed, ThetaStructure};
use crate::theta::dim_two::ThetaPointDimTwo;

/// Convert from two dim two points to one dim four point
/// with the correct change of basis applied
fn dim_two_to_dim_four<Fp: FqTrait>(
    P: &ThetaPointDimTwo<Fp>,
    Q: &ThetaPointDimTwo<Fp>,
) -> ThetaPoint<Fp> {
    let mut coords = [Fp::ZERO; 16];
    let (px, py, pz, pt) = P.coords();
    let (qx, qy, qz, qt) = Q.coords();

    coords[0] = px * qx;
    coords[8] = px * qz;
    coords[4] = px * qy;
    coords[12] = px * qt;
    coords[2] = pz * qx;
    coords[10] = pz * qz;
    coords[6] = pz * qy;
    coords[14] = pz * qt;
    coords[1] = py * qx;
    coords[9] = py * qz;
    coords[5] = py * qy;
    coords[13] = py * qt;
    coords[3] = pt * qx;
    coords[11] = pt * qz;
    coords[7] = pt * qy;
    coords[15] = pt * qt;
    let mut R = ThetaPoint::new_from_array(coords);
    R.set_apply_gluing_basis_change();
    R
}

/// helper for solve_HIIP_gluing where we need this on 8 elements, rather
/// than on a point.
fn proj_batch_inversion<const N: usize, Fq: FqTrait>(v: &mut [Fq; N]) {
    // TODO: a very similar thing is needed for theta points
    // I should refactor this and reuse it in both places.
    let mut multiples = [Fq::ONE; N];
    multiples[0] = v[0];
    for i in 1..N {
        multiples[i] = multiples[i - 1] * v[i];
    }

    let mut inverses = [Fq::ONE; N];
    for i in 1..N {
        inverses[i] = inverses[i - 1] * v[N - i];
    }

    v[0] = inverses[N - 1];
    for i in 1..N {
        v[i] = inverses[N - 1 - i] * multiples[i - 1];
    }
}

/// Compute the inverse dual null point of the codomain
fn solve_HIIP_gluing<Fq: FqTrait>(
    kernel: &[ThetaPoint<Fq>; 4],
) -> (ThetaPoint<Fq>, ThetaPoint<Fq>) {
    // leg_4: uses HSK_8[0] i.e. kernel[0], indices 0, 2, 4, 6
    // 0--(4,2)--6
    let leg_4 = [
        kernel[0][4] * kernel[0][6],
        kernel[0][0] * kernel[0][6],
        kernel[0][0] * kernel[0][2],
    ];

    // leg_8: uses HSK_8[1] i.e. kernel[1], indices 0, 1, 8, 9
    // 0--(8,1)--9
    let leg_8 = [
        kernel[1][8] * kernel[1][9],
        kernel[1][0] * kernel[1][9],
        kernel[1][0] * kernel[1][1],
    ];

    // leg_12: uses HSK_8[2] i.e. kernel[2], indices 0, 3, 12, 15
    // 0--(12,3)--15
    let leg_12 = [
        kernel[2][12] * kernel[2][15],
        kernel[2][0] * kernel[2][15],
        kernel[2][0] * kernel[2][3],
    ];

    // Build the 10 compressed coordinates of inv_theta_null.
    let mut inv_null_comp = ThetaPointCompressed::default();

    let tmp_8_12 = leg_8[0] * leg_12[0];
    inv_null_comp[0] = tmp_8_12 * leg_4[0];
    inv_null_comp[2] = tmp_8_12 * leg_4[1];
    inv_null_comp[5] = tmp_8_12 * leg_4[2];

    let tmp_4_12 = leg_4[0] * leg_12[0];
    inv_null_comp[1] = tmp_4_12 * leg_8[1];
    inv_null_comp[7] = tmp_4_12 * leg_8[2];

    let tmp_4_8 = leg_4[0] * leg_8[0];
    inv_null_comp[3] = tmp_4_8 * leg_12[1];
    inv_null_comp[9] = tmp_4_8 * leg_12[2];

    // Expand compressed -> full 16 coords via the duplication pattern
    let inv_theta_null: ThetaPoint<Fq> = inv_null_comp.to_theta_point();

    // Build T3: 8 values using kernel[0] (= H(S(T3))) and inv_theta_null
    // Note T3[6] = kernel[0][10] * inv_theta_null[10] will be zero
    // since inv_theta_null[10] = inv_null_comp[4] = 0
    let mut t3: [Fq; 8] = [
        kernel[0][0] * inv_theta_null[0],
        kernel[0][1] * inv_theta_null[1],
        kernel[0][2] * inv_theta_null[2],
        kernel[0][3] * inv_theta_null[3],
        kernel[0][8] * inv_theta_null[8],
        kernel[0][9] * inv_theta_null[9],
        kernel[0][10] * inv_theta_null[10],
        kernel[0][15] * inv_theta_null[15],
    ];

    // We want to correct the value of t3[6] which requires an inversion
    //
    // t3[6] =  t3[0] * kernel[3][2] * inv_theta_null[2] / (kernel[3][8] * inv_theta_null[8])
    //
    // But we can instead clear the denominator correction_den = kernel[3][8] * inv_theta_null[8]
    // with 7M instead of 1I
    t3[6] = t3[0] * kernel[3][2] * inv_theta_null[2];

    let correction_den = kernel[3][8] * inv_theta_null[8];
    t3[0] *= correction_den;
    t3[1] *= correction_den;
    t3[2] *= correction_den;
    t3[3] *= correction_den;
    t3[4] *= correction_den;
    t3[5] *= correction_den;
    // Do not multiply t3[6] b y the factor here here
    t3[7] *= correction_den;

    // Batch invert the 8 T3 values
    proj_batch_inversion(&mut t3);

    // Expand to 16 with duplication: inv_T3[8i+j] = inv_T3[8i+j+4] = inv_small_t3[4i+j]
    let mut inv_t3_coords = [Fq::ZERO; 16];
    for i in 0..2 {
        for j in 0..4 {
            inv_t3_coords[8 * i + j] = t3[4 * i + j];
            inv_t3_coords[8 * i + j + 4] = t3[4 * i + j];
        }
    }
    let inv_t3 = ThetaPoint::new_from_array(inv_t3_coords);

    (inv_theta_null, inv_t3)
}

/// Compute the final gluing 2-isogeny from A_2 x A_2 -> A_4
pub fn dim_four_gluing_isogeny<Fp: FqTrait>(
    diagonal_images: &DiagonalImages<Fp>,
) -> (ThetaStructure<Fp>, Vec<ThetaPoint<Fp>>) {
    // Build and transform the gluing kernel
    let gluing_kernel: [ThetaPoint<Fp>; 4] = diagonal_images.gluing_kernel.map(|[u, v]| {
        let mut p = dim_two_to_dim_four(&u, &v);
        p.set_square();
        p.set_hadamard();
        p
    });

    let (inv_null_point_dual, inv_T3) = solve_HIIP_gluing(&gluing_kernel);

    const NON_ZERO: [u8; 10] = [0, 1, 2, 3, 4, 6, 8, 9, 12, 15];
    let mut non_zero_coordinates: [Fp; 10] = NON_ZERO.map(|i| inv_null_point_dual[i as usize]);
    proj_batch_inversion(&mut non_zero_coordinates);
    let mut codomain_null_point = ThetaPoint::default();
    for i in 0..10 {
        codomain_null_point[NON_ZERO[i] as usize] = non_zero_coordinates[i];
    }
    codomain_null_point.set_hadamard();
    let codomain = ThetaStructure::new_from_null_point(codomain_null_point);

    // Compute the images the idea is we take t1 made from the pairs of points in image_plus
    // and t2 from pairs of points in image_minus. Then the final image is t3 = t1 \star t2
    // then image = H(inv_T3 \star H(t1 \star t2))
    let chain_kernel: Vec<ThetaPoint<Fp>> = diagonal_images
        .image_plus
        .iter()
        .zip(diagonal_images.image_minus.iter())
        .map(|([plus_u, plus_v], [minus_u, minus_v])| {
            let t1 = dim_two_to_dim_four(plus_u, plus_v);
            let t2 = dim_two_to_dim_four(minus_u, minus_v);
            let mut t3 = t1.coordinate_multiply(&t2);
            t3.set_hadamard();
            t3.set_coordinate_multiply(&inv_T3);
            t3.set_hadamard();
            t3
        })
        .collect();

    (codomain, chain_kernel)
}
