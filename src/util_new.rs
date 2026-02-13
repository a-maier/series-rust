use num_traits::Zero;

pub(crate) fn trim_zero<C: Zero>(coeffs: &mut Vec<C>) -> usize {
    trim_end_zero(coeffs);
    trim_start_zero(coeffs)
}

fn trim_start_zero<C: Zero>(coeffs: &mut Vec<C>) -> usize {
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
