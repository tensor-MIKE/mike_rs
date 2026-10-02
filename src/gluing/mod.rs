use fp2::traits::Fp2 as Fp2Trait;
use fp2::traits::Fq as FqTrait;

pub mod diagonal_gluing;
pub mod dim_four_gluing;
pub mod scholten_gluing;

use crate::elliptic::curve::Curve;
use crate::elliptic::jac_point::Point;
use crate::mike::mike_impl::KernelPoints;
use crate::theta::dim_four::{ThetaPoint, ThetaStructure};

use crate::gluing::diagonal_gluing::diagonal_isogeny;
use crate::gluing::dim_four_gluing::dim_four_gluing_isogeny;
use crate::gluing::scholten_gluing::scholten_isogeny;

/// Compute three gluing isogenies to get a dimension four codomain and two points
/// which are the kernel of a 2^(e - 3) isogeny
pub fn gluing_isogeny<Fp: FqTrait, Fp2: Fp2Trait<BaseField = Fp>>(
    E: &Curve<Fp2>,
    kernel_points: &KernelPoints<Fp2>,
    image_points: &[[Point<Fp2>; 2]],
) -> (ThetaStructure<Fp>, Vec<ThetaPoint<Fp>>) {
    // Scholten isogeny from E x E^sigma -> A_2 x A_2
    let scholten_images = scholten_isogeny(E, kernel_points, image_points);

    // Diagonal isogeny from A_2 x A_2 -> A_2 x A_2
    let diagonal_images = diagonal_isogeny(&scholten_images);

    // Final gluing isogeny from A_2 x A_2 -> A_4
    let (gluing_codomain, chain_kernel) = dim_four_gluing_isogeny(&diagonal_images);
    (gluing_codomain, chain_kernel)
}
