#[cfg(test)]
pub mod test_data {
    use crate::fields::MikeTestFp2 as Fp2;

    // A1 is a QR in the field
    const A1_RE_BYTES: [u8; 32] = [
        152, 137, 139, 24, 148, 181, 15, 114, 226, 232, 62, 52, 255, 78, 160, 252, 198, 152, 194,
        208, 236, 58, 127, 192, 252, 125, 128, 226, 56, 67, 46, 0,
    ];
    const A1_IM_BYTES: [u8; 32] = [
        34, 58, 173, 76, 208, 170, 127, 6, 145, 18, 209, 105, 136, 253, 212, 227, 117, 246, 231,
        59, 134, 209, 52, 42, 76, 18, 96, 76, 128, 114, 51, 4,
    ];
    pub const A1: Fp2 = Fp2::const_decode_no_check(&A1_RE_BYTES, &A1_IM_BYTES);

    // A2 is a NQR in the field
    const A2_RE_BYTES: [u8; 32] = [
        198, 232, 251, 8, 104, 243, 12, 83, 239, 92, 144, 156, 82, 232, 69, 217, 79, 60, 217, 7,
        120, 94, 249, 178, 204, 66, 178, 246, 69, 197, 227, 3,
    ];
    const A2_IM_BYTES: [u8; 32] = [
        88, 201, 129, 250, 30, 103, 239, 169, 169, 134, 195, 94, 168, 10, 73, 96, 52, 214, 17, 215,
        215, 212, 51, 153, 1, 236, 52, 84, 43, 174, 118, 0,
    ];
    pub const A2: Fp2 = Fp2::const_decode_no_check(&A2_RE_BYTES, &A2_IM_BYTES);
}

#[cfg(test)]
pub use test_data::*;
