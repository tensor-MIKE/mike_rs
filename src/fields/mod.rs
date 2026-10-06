#[cfg(feature = "conservative_parameters")]
mod constants {
    // p374.117
    pub const MIKE_I_MODULUS: [u64; 6] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x1D3F_FFFF_FFFF_FFFF,
    ];

    // p566.77
    pub const MIKE_III_MODULUS: [u64; 9] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x133F_FFFF_FFFF_FFFF,
    ];

    // p758.41
    pub const MIKE_V_MODULUS: [u64; 12] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x0A3F_FFFF_FFFF_FFFF,
    ];
}

#[cfg(not(feature = "conservative_parameters"))]
mod constants {
    // p308.644
    pub const MIKE_I_MODULUS: [u64; 5] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x278F_FFFF_FFFF_FFFF,
    ];

    // p474.593
    pub const MIKE_III_MODULUS: [u64; 8] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x0000_0009_43FF_FFFF,
    ];

    // p628.317
    pub const MIKE_V_MODULUS: [u64; 10] = [
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0xFFFF_FFFF_FFFF_FFFF,
        0x13CF_FFFF_FFFF_FFFF,
    ];
}

// Re-export so they are directly accessible at the module level
pub use constants::*;

#[cfg(test)]
pub const MIKE_TEST_MODULUS: [u64; 4] = [
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0x04FF_FFFF_FFFF_FFFF,
];

#[cfg(test)]
fp2::define_fp_core!(typename = MikeTestFp, modulus = MIKE_TEST_MODULUS,);

#[cfg(test)]
fp2::define_fp2_from_type!(
    typename = MikeTestFp2,
    base_field = MikeTestFp,
    use_sum_of_products = true,
);

fp2::define_fp_core!(typename = MikeIFp, modulus = MIKE_I_MODULUS,);

fp2::define_fp_core!(typename = MikeIIIFp, modulus = MIKE_III_MODULUS,);

fp2::define_fp_core!(typename = MikeVFp, modulus = MIKE_V_MODULUS,);

fp2::define_fp2_from_type!(
    typename = MikeIFp2,
    base_field = MikeIFp,
    use_sum_of_products = true,
);

fp2::define_fp2_from_type!(
    typename = MikeIIIFp2,
    base_field = MikeIIIFp,
    use_sum_of_products = true,
);

fp2::define_fp2_from_type!(
    typename = MikeVFp2,
    base_field = MikeVFp,
    use_sum_of_products = true,
);

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! mike_tests {
        (fp2: $mod:ident, $fp2:ident, $modulus:ident, $nqr:expr) => {
            mod $mod {
                use super::{$fp2, $modulus};
                fp2::define_fp2_tests!($fp2, $modulus, $nqr);
            }
        };
        (fp: $mod:ident, $fp:ident $(, cfg($($attr:tt)*))?) => {
            $(#[cfg($($attr)*)])?
            mod $mod {
                use super::$fp;
                fp2::define_fp_tests!($fp);
            }
        };
    }

    mike_tests!(fp: mike_test_fp,   MikeTestFp);
    mike_tests!(fp: mike_i_fp,   MikeIFp);
    mike_tests!(fp: mike_iii_fp, MikeIIIFp);
    mike_tests!(fp: mike_v_fp,   MikeVFp);
    mike_tests!(fp2: mike_test_fp2,   MikeTestFp2,   MIKE_TEST_MODULUS,   5);

    #[cfg(feature = "conservative_parameters")]
    mike_tests!(fp2: mike_i_fp2,   MikeIFp2,   MIKE_I_MODULUS,   2);
    #[cfg(not(feature = "conservative_parameters"))]
    mike_tests!(fp2: mike_i_fp2,   MikeIFp2,   MIKE_I_MODULUS,   2);

    #[cfg(feature = "conservative_parameters")]
    mike_tests!(fp2: mike_iii_fp2, MikeIIIFp2, MIKE_III_MODULUS, 2);
    #[cfg(not(feature = "conservative_parameters"))]
    mike_tests!(fp2: mike_iii_fp2, MikeIIIFp2, MIKE_III_MODULUS, 15);

    #[cfg(feature = "conservative_parameters")]
    mike_tests!(fp2: mike_v_fp2,   MikeVFp2,   MIKE_V_MODULUS,   2);
    #[cfg(not(feature = "conservative_parameters"))]
    mike_tests!(fp2: mike_v_fp2,   MikeVFp2,   MIKE_V_MODULUS,   4);
}
