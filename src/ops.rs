
/// Trait for the natural logarithm
pub trait Ln {
    fn ln(self) -> Self;
}

impl Ln for f64 {
    fn ln(self) -> Self {
        <f64>::ln(self)
    }
}

impl Ln for f32 {
    fn ln(self) -> Self {
        <f32>::ln(self)
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

impl Exp for f32 {
    type Output = Self;

    fn exp(self) -> Self::Output {
        <f32>::exp(self)
    }
}
