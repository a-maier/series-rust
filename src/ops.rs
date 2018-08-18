
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
    fn exp(self) -> Self;
}

impl Exp for f64 {
    fn exp(self) -> Self {
        <f64>::exp(self)
    }
}

impl Exp for f32 {
    fn exp(self) -> Self {
        <f32>::exp(self)
    }
}
