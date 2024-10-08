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

/// Karatsuba multiplication
pub trait KaratsubaMul<Rhs> {
    type Output;

    /// Calculate `self * rhs` using Karatsuba multiplication.
    ///
    /// `min_size` marks the threshold size below which naive
    /// multiplication should be used.
    fn karatsuba_mul(self, rhs: Rhs, min_size: usize) -> Self::Output;
}

// Logarithm for series starting with var^0
pub(crate) trait LnVarFree {
    type Output;
    fn ln_var_free(self) -> Self::Output;
}

pub(crate) trait ExpCoeff {
    type Output;
    fn exp_coeff(&self) -> Self::Output;
}
