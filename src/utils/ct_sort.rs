// ==============================================================================
// u8 constant time operations
// ==============================================================================

/// Returns `u32::MAX` if `a < b` (unsigned), `0` otherwise.
#[inline(always)]
fn ct_lt_u8(a: u8, b: u8) -> u32 {
    let borrow_bit = (a as u32).wrapping_sub(b as u32) >> 8 & 1;
    borrow_bit.wrapping_neg()
}

/// Returns `u32::MAX` if `a == b`, `0` otherwise.
#[inline(always)]
fn ct_eq_u8(a: u8, b: u8) -> u32 {
    let x = (a ^ b) as u32;
    let nonzero = (x | x.wrapping_neg()) >> 31;
    !(nonzero).wrapping_neg()
}

// ==============================================================================
// Byte slice constant time operations
// ==============================================================================

/// Returns `(less, equal)` masks where each is `u32::MAX` or `0`.
/// `less`  = u32::MAX iff a < b (as little-endian integers)
/// `equal` = u32::MAX iff a == b
/// Requires `a` and `b` to have the same length.
#[inline(always)]
fn ct_lt_and_eq_le_bytes(a: &[u8], b: &[u8]) -> (u32, u32) {
    debug_assert_eq!(a.len(), b.len());
    let mut result: u32 = 0;
    let mut all_equal_so_far: u32 = u32::MAX;

    for (&ai, &bi) in a.iter().zip(b.iter()).rev() {
        let less = ct_lt_u8(ai, bi);
        let equal = ct_eq_u8(ai, bi);
        result |= all_equal_so_far & less;
        all_equal_so_far &= equal;
    }

    (result, all_equal_so_far)
}

/// Returns `u32::MAX` if the little-endian integer `a < b`, `0` otherwise.
#[inline(always)]
fn ct_lt_le_bytes(a: &[u8], b: &[u8]) -> u32 {
    ct_lt_and_eq_le_bytes(a, b).0
}

/// Sets each byte of `a` to either its current value or the corresponding
/// byte of `b`, depending on whether `ctl` is `0u32` or `u32::MAX`.
/// Requires `a` and `b` to have the same length.
#[inline(always)]
fn ct_select_le_bytes(ctl: u32, a: &mut [u8], b: &[u8]) {
    debug_assert_eq!(a.len(), b.len());
    let mask = ctl as u8;
    for (ai, &bi) in a.iter_mut().zip(b.iter()) {
        *ai ^= mask & (*ai ^ bi);
    }
}

/// Finds the lexicographically smallest value of an array of encoded Fp2 values
pub fn ct_find_smallest_in_array<T: Copy + AsRef<[u8]> + AsMut<[u8]>>(input: &[T]) -> T {
    let mut output = input[0];
    for value in &input[1..] {
        let ctl = ct_lt_le_bytes(value.as_ref(), output.as_ref());
        ct_select_le_bytes(ctl, output.as_mut(), value.as_ref())
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    fn le8(v: u64) -> [u8; 8] {
        v.to_le_bytes()
    }

    #[test]
    fn test_ct_lt_u8() {
        assert_eq!(ct_lt_u8(0, 1), u32::MAX);
        assert_eq!(ct_lt_u8(1, 1), 0);
        assert_eq!(ct_lt_u8(2, 1), 0);
        assert_eq!(ct_lt_u8(0, 255), u32::MAX);
        assert_eq!(ct_lt_u8(255, 0), 0);
    }

    #[test]
    fn test_ct_eq_u8() {
        assert_eq!(ct_eq_u8(0, 0), u32::MAX);
        assert_eq!(ct_eq_u8(1, 1), u32::MAX);
        assert_eq!(ct_eq_u8(0, 1), 0);
        assert_eq!(ct_eq_u8(255, 255), u32::MAX);
    }

    #[test]
    fn test_ct_lt_equal() {
        assert_eq!(ct_lt_le_bytes(&le8(1), &le8(1)), 0);
        assert_eq!(ct_lt_le_bytes(&le8(0), &le8(0)), 0);
    }

    #[test]
    fn test_ct_lt_simple() {
        assert_eq!(ct_lt_le_bytes(&le8(0), &le8(1)), u32::MAX);
        assert_eq!(ct_lt_le_bytes(&le8(2), &le8(1)), 0);
    }

    #[test]
    fn test_ct_lt_msb_differs() {
        assert_eq!(
            ct_lt_le_bytes(&le8(0x00FF_0000), &le8(0x0100_0000)),
            u32::MAX
        );
        assert_eq!(ct_lt_le_bytes(&le8(0x0100_0000), &le8(0x00FF_0000)), 0);
    }

    #[test]
    fn test_ct_lt_large_values() {
        assert_eq!(
            ct_lt_le_bytes(&le8(0x00FF_0000), &le8(0x0100_0000)),
            u32::MAX
        );
        assert_eq!(ct_lt_le_bytes(&le8(0x0100_0000), &le8(0x00FF_0000)), 0);
        assert_eq!(
            ct_lt_le_bytes(&le8(0xFFFF_FFFE), &le8(0xFFFF_FFFF)),
            u32::MAX
        );
        assert_eq!(ct_lt_le_bytes(&le8(0xFFFF_FFFF), &le8(0xFFFF_FFFE)), 0);
        assert_eq!(ct_lt_le_bytes(&le8(0x0000_0200), &le8(0x0000_0101)), 0);
        assert_eq!(
            ct_lt_le_bytes(&le8(0x0000_0101), &le8(0x0000_0200)),
            u32::MAX
        );
    }

    #[test]
    fn test_ct_lt_and_eq_le_bytes() {
        assert_eq!(ct_lt_and_eq_le_bytes(&le8(1), &le8(2)), (u32::MAX, 0));
        assert_eq!(ct_lt_and_eq_le_bytes(&le8(2), &le8(1)), (0, 0));
        assert_eq!(ct_lt_and_eq_le_bytes(&le8(2), &le8(2)), (0, u32::MAX));
    }

    #[test]
    fn test_ct_select_le_bytes() {
        let (a, b) = (le8(111), le8(222));

        let mut out = a;
        ct_select_le_bytes(0, &mut out, &b);
        assert_eq!(out, a);

        let mut out = a;
        ct_select_le_bytes(u32::MAX, &mut out, &b);
        assert_eq!(out, b);
    }

    #[test]
    fn test_ct_find_smallest_in_array() {
        assert_eq!(ct_find_smallest_in_array(&[le8(3), le8(1), le8(2)]), le8(1));
        assert_eq!(ct_find_smallest_in_array(&[le8(5)]), le8(5));
    }
}
