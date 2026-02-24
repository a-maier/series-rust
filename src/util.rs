use num_traits::Zero;

pub(crate) fn trim_zero<C: Zero>(coeffs: &mut Vec<C>) -> usize {
    trim_end_zero(coeffs);
    trim_start_zero(coeffs)
}

pub(crate) fn trim_start_zero<C: Zero>(coeffs: &mut Vec<C>) -> usize {
    let leading_zeros_end = coeffs
        .iter()
        .position(|c| !c.is_zero())
        .unwrap_or(coeffs.len());
    coeffs.drain(..leading_zeros_end);
    leading_zeros_end
}

fn trim_end_zero<C: Zero>(coeffs: &mut Vec<C>) -> usize {
    let old_len = coeffs.len();
    let trailing_zeros_start = coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .map(|p| p + 1)
        .unwrap_or(0);
    coeffs.truncate(trailing_zeros_start);
    old_len - trailing_zeros_start
}

pub(crate) fn trim_slice_zero<C: Zero>(coeffs: &mut &[C]) -> usize {
    trim_slice_end_zero(coeffs);
    trim_slice_start_zero(coeffs)
}

pub(crate) fn trim_slice_start_zero<C: Zero>(coeffs: &mut &[C]) -> usize {
    let leading_zeros_end = coeffs
        .iter()
        .position(|c| !c.is_zero())
        .unwrap_or(coeffs.len());
    *coeffs = &coeffs[leading_zeros_end..];
    leading_zeros_end
}

fn trim_slice_end_zero<C: Zero>(coeffs: &mut &[C]) -> usize {
    let old_len = coeffs.len();
    let trailing_zeros_start = coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .map(|p| p + 1)
        .unwrap_or(0);
    *coeffs = &coeffs[..trailing_zeros_start];
    old_len - trailing_zeros_start
}

pub(crate) trait NumDisplay {
    fn starts_with_minus(&self) -> bool;
    fn abs(self) -> Self;
}

macro_rules! impl_display_traits_signed {
    ($($t:ty), *) => {
        $(
            impl NumDisplay for $t {
                fn starts_with_minus(&self) -> bool {
                    *self < Zero::zero()
                }

                fn abs(self) -> $t {
                    -self
                }
            }
        )*
    }
}

impl_display_traits_signed!(i8, i16, i32, i64, i128, isize, f32, f64);

macro_rules! impl_display_traits_unsigned {
    ($($t:ty), *) => {
        $(
            impl NumDisplay for $t {
                fn starts_with_minus(&self) -> bool {
                    false
                }

                fn abs(self) -> $t {
                    self
                }
            }
        )*
    }
}

impl_display_traits_unsigned!(u8, u16, u32, u64, u128, usize);
