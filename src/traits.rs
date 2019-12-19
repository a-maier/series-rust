// Helper traits
use crate::Coeff;
use crate::slice::SeriesSlice;

/// Multiplicative inverse
pub trait MulInverse {
    type Output;

    fn mul_inverse(self) -> Self::Output;
}

/// View a slice of the original object
pub trait AsSlice<'a, T> {
    type Output;

    fn as_slice(&'a self, t: T) -> Self::Output;
}

// Logarithm for series starting with var^0
pub(crate) trait LnVarFree {
    type Output;
    fn ln_var_free(self) -> Self::Output;
}

// Spurious traits, needed for rust >= 1.28
// direct implementation used to work in 1.24
pub(crate) trait AddAssignHelper<Var, C: Coeff> {
    fn truncate_cutoff_pow<'a>(&mut self, other: SeriesSlice<'a, Var, C>);
    fn add_overlap<'a>(&mut self, other: SeriesSlice<'a, Var, C>);
    fn num_leading<'a>(&mut self, other: SeriesSlice<'a, Var, C>) -> usize;
}

pub(crate) trait ExpCoeff {
    type Output;
    fn exp_coeff(&self) -> Self::Output;
}
