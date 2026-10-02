pub(crate) mod isogeny_chain;
pub(crate) mod mike_impl;
pub(crate) mod public_parameters;
pub(crate) mod test_data;

use crate::fields::{MikeIFp, MikeIFp2, MikeIIIFp, MikeIIIFp2, MikeVFp, MikeVFp2};
use mike_impl::Mike;

// Generics required are the Fp, and Fp2 fields as well as the length of the secret key in bytes
type MikeI = Mike<MikeIFp, MikeIFp2, { public_parameters::N_I }>;
type MikeIII = Mike<MikeIIIFp, MikeIIIFp2, { public_parameters::N_III }>;
type MikeV = Mike<MikeVFp, MikeVFp2, { public_parameters::N_V }>;

pub const MIKE_I: MikeI = MikeI::new(&public_parameters::MIKE_I_PARAMS);
pub const MIKE_III: MikeIII = MikeIII::new(&public_parameters::MIKE_III_PARAMS);
pub const MIKE_V: MikeV = MikeV::new(&public_parameters::MIKE_V_PARAMS);
