/// Traits for common functions

/// Trait for the natural logarithm
pub trait Ln {
    type Output;

    fn ln(self) -> Self::Output;
}

impl Ln for f64 {
    type Output = f64;

    fn ln(self) -> Self::Output {
        <f64>::ln(self)
    }
}

impl<'a> Ln for &'a f64 {
    type Output = f64;

    fn ln(self) -> Self::Output {
        (*self).ln()
    }
}

impl Ln for f32 {
    type Output = f32;

    fn ln(self) -> Self::Output {
        <f32>::ln(self)
    }
}

impl<'a> Ln for &'a f32 {
    type Output = f32;

    fn ln(self) -> Self::Output {
        (*self).ln()
    }
}

/// Trait for the exponential function
pub trait Exp {
    type Output;

    fn exp(self) -> Self::Output;
}

impl Exp for f64 {
    type Output = Self;

    fn exp(self) -> Self::Output {
        <f64>::exp(self)
    }
}

impl<'a> Exp for &'a f64 {
    type Output = f64;

    fn exp(self) -> Self::Output {
        (*self).exp()
    }
}

impl Exp for f32 {
    type Output = Self;

    fn exp(self) -> Self::Output {
        <f32>::exp(self)
    }
}

impl<'a> Exp for &'a f32 {
    type Output = f32;

    fn exp(self) -> Self::Output {
        (*self).exp()
    }
}

pub use num_traits::Pow;
