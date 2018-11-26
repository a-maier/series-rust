
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

/// Trait for powers
pub trait Pow<T> {
    type Output;

    fn pow(self, T) -> Self::Output;
}

impl Pow<f64> for f64 {
    type Output = Self;

    fn pow(self, exponent: f64) -> Self::Output {
        <f64>::powf(self, exponent)
    }
}

impl<'a> Pow<f64> for &'a f64 {
    type Output = f64;

    fn pow(self, exponent: f64) -> Self::Output {
        (*self).pow(exponent)
    }
}

impl<'a> Pow<&'a f64> for f64 {
    type Output = f64;

    fn pow(self, exponent: &'a f64) -> Self::Output {
        self.pow(*exponent)
    }
}

impl<'a, 'b> Pow<&'a f64> for &'b f64 {
    type Output = f64;

    fn pow(self, exponent: &'a f64) -> Self::Output {
        (*self).pow(*exponent)
    }
}

impl Pow<i32> for f64 {
    type Output = Self;

    fn pow(self, exponent: i32) -> Self::Output {
        <f64>::powi(self, exponent)
    }
}

impl<'a> Pow<i32> for &'a f64 {
    type Output = f64;

    fn pow(self, exponent: i32) -> Self::Output {
        (*self).pow(exponent)
    }
}

impl<'a> Pow<&'a i32> for f64 {
    type Output = f64;

    fn pow(self, exponent: &'a i32) -> Self::Output {
        self.pow(*exponent)
    }
}

impl<'a, 'b> Pow<&'a i32> for &'b f64 {
    type Output = f64;

    fn pow(self, exponent: &'a i32) -> Self::Output {
        (*self).pow(*exponent)
    }
}

impl Pow<f32> for f32 {
    type Output = Self;

    fn pow(self, exponent: f32) -> Self::Output {
        <f32>::powf(self, exponent)
    }
}

impl<'a> Pow<f32> for &'a f32 {
    type Output = f32;

    fn pow(self, exponent: f32) -> Self::Output {
        (*self).pow(exponent)
    }
}

impl<'a> Pow<&'a f32> for f32 {
    type Output = f32;

    fn pow(self, exponent: &'a f32) -> Self::Output {
        self.pow(*exponent)
    }
}

impl<'a, 'b> Pow<&'a f32> for &'b f32 {
    type Output = f32;

    fn pow(self, exponent: &'a f32) -> Self::Output {
        (*self).pow(*exponent)
    }
}

impl Pow<i32> for f32 {
    type Output = Self;

    fn pow(self, exponent: i32) -> Self::Output {
        <f32>::powi(self, exponent)
    }
}

impl<'a> Pow<i32> for &'a f32 {
    type Output = f32;

    fn pow(self, exponent: i32) -> Self::Output {
        (*self).pow(exponent)
    }
}

impl<'a> Pow<&'a i32> for f32 {
    type Output = f32;

    fn pow(self, exponent: &'a i32) -> Self::Output {
        self.pow(*exponent)
    }
}

impl<'a, 'b> Pow<&'a i32> for &'b f32 {
    type Output = f32;

    fn pow(self, exponent: &'a i32) -> Self::Output {
        (*self).pow(*exponent)
    }
}

// TODO: Pow for all integer types
