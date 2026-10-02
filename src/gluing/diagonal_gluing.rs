use crate::gluing::scholten_gluing::ScholtenImages;
use crate::theta::dim_two::{ThetaPointDimTwo, ThetaStructureDimTwo};
use fp2::traits::Fq as FqTrait;

pub struct DiagonalImages<Fp: FqTrait> {
    pub gluing_kernel: [[ThetaPointDimTwo<Fp>; 2]; 4],
    pub image_plus: Vec<[ThetaPointDimTwo<Fp>; 2]>,
    pub image_minus: Vec<[ThetaPointDimTwo<Fp>; 2]>,
}

/// Compute the second gluing step, a diagonal isogeny A x A -> A x A, for A a dimension two theta structure
pub fn diagonal_isogeny<Fp: FqTrait>(scholten_images: &ScholtenImages<Fp>) -> DiagonalImages<Fp> {
    // We package the first of each pair of points into U_images, and the second into
    // V_images. However, when x = 5 mod 8, we need to swap this decision in constant
    // time.

    // For each pair of theta points, put the first into U and the second into V, then, if x = 5 mod 8
    // we swap everything
    let n = scholten_images.image_plus.len();
    let mut U_images: Vec<ThetaPointDimTwo<Fp>> = Vec::with_capacity(3 + 2 * n);
    let mut V_images: Vec<ThetaPointDimTwo<Fp>> = Vec::with_capacity(3 + 2 * n);
    for &[u, v] in &scholten_images.gluing_kernel {
        U_images.push(u);
        V_images.push(v);
    }
    for &[u, v] in &scholten_images.image_plus {
        U_images.push(u);
        V_images.push(v);
    }
    for &[u, v] in &scholten_images.image_minus {
        U_images.push(u);
        V_images.push(v);
    }

    // Compute the isogeny and push all the points through. Note that we treat phi_U as a dual by reversing the Hadamard bools
    let [U1, U2, V1, V2] = scholten_images.diagonal_kernel;
    let codomain_u = ThetaStructureDimTwo::two_isogeny(&U1, &U2, &mut U_images, [true, false]);
    let codomain_v = ThetaStructureDimTwo::two_isogeny(&V1, &V2, &mut V_images, [false, true]);

    // Lastly, we extend the gluing basis with a final differential addition
    let u_21 = codomain_u.diff_addition(&U_images[2], &U_images[1], &U_images[0]);
    let v_21 = codomain_v.diff_addition(&V_images[2], &V_images[1], &V_images[0]);

    // Reconstruct all the various images into the right parts for the final gluing
    DiagonalImages {
        gluing_kernel: [
            [U_images[0], V_images[0]],
            [U_images[1], V_images[1]],
            [U_images[2], V_images[2]],
            [u_21, v_21],
        ],
        image_plus: (0..n).map(|i| [U_images[3 + i], V_images[3 + i]]).collect(),
        image_minus: (0..n)
            .map(|i| [U_images[3 + n + i], V_images[3 + n + i]])
            .collect(),
    }
}
