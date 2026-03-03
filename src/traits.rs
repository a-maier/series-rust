use std::fmt::Display;

/// Multiplicative inverse
pub trait MulInverse {
    type Output;

    fn mul_inverse(self) -> Self::Output;
}

/// View a slice of the original object
pub trait AsSlice<T> {
    type Output;

    fn as_slice(self, t: T) -> Self::Output;
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

/// Split an object into the leading sign and a remainder
///
/// This trait is needed for the [Display] implementations of [Series]
/// and [Polynomials](Polynomial). Note that it should usually be
/// implemented on _reference types_. `split_sign` does not take a
/// reference to facilitate an efficient implementation on complex
/// types: it makes it easier to define `Signless` as a lightweight reference type:
/// ```
/// # use series::{SplitSign, Sign};
/// # struct SomeComplexType();
/// # struct LightweightReference<'a>(std::marker::PhantomData<&'a ()>);
/// impl<'a> SplitSign for &'a SomeComplexType {
///     type Signless = LightweightReference<'a>; // can reuse lifetime 'a here
///
///     fn split_sign(self) -> (Sign, Self::Signless) {
///         todo!()
///     }
/// }
/// ```
/// # Example
///
/// ```
/// # use std::cmp::Ordering;
/// # use series::{Sign, SplitSign};
/// // Define a custom complex number type
/// // with integer real and imaginary parts
/// #[derive(Copy, Clone)]
/// struct Complex {
///     re: i32,
///     im: i32,
/// }
///
/// impl SplitSign for &Complex {
///     type Signless = Complex;
///
///     fn split_sign(self) -> (Sign, Self::Signless) {
///         let Complex{re, im} = self;
///         // split off the sign of the real part, unless it's zero
///         match re.cmp(&0) {
///             Ordering::Greater => (Sign::Plus, *self),
///             Ordering::Less => (Sign::Minus, Complex{re: -re, im: *im}),
///             Ordering::Equal => if *im >= 0 {
///                 (Sign::Plus, *self)
///             } else {
///                 (Sign::Minus, Complex{re: *re, im: -im})
///             }
///         }
///     }
/// }
/// ```
pub trait SplitSign {
    type Signless;

    fn split_sign(self) -> (Sign, Self::Signless);
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Sign {
    Plus,
    Minus,
}

impl Display for Sign {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Sign::Plus => '+',
            Sign::Minus => '-',
        }.fmt(f)
    }
}

macro_rules! impl_split_sign_signed_int {
    ($($t:ty), *) => {
        $(
            impl SplitSign for $t {
                type Signless = $t;

                fn split_sign(self) -> (Sign, Self::Signless) {
                    if self > 0 {
                        (Sign::Plus, self)
                    } else {
                        (Sign::Minus, -self)
                    }
                }
            }
        )*
    };
}

impl_split_sign_signed_int!(i8, i16, i32, i64, i128, isize);

macro_rules! impl_split_sign_unsigned_int {
    ($($t:ty), *) => {
        $(
            impl SplitSign for $t {
                type Signless = $t;

                fn split_sign(self) -> (Sign, Self::Signless) {
                    (Sign::Plus, self)
                }
            }
        )*
    };
}

impl_split_sign_unsigned_int!(u8, u16, u32, u64, u128, usize);

macro_rules! impl_split_sign_float {
    ($($t:ty), *) => {
        $(
            impl SplitSign for $t {
                type Signless = $t;

                fn split_sign(self) -> (Sign, Self::Signless) {
                    // use `<` for comparison so that NaN ends up with a Plus sign
                    if self < 0.0 {
                        (Sign::Minus, -self)
                    } else {
                        (Sign::Plus, self)
                    }
                }
            }
        )*
    };
}

impl_split_sign_float!(f32, f64);

/// Check if a bracket is needed when used as a coefficient in a [Series] or [Polynomial]
///
/// This trait is needed for the [Display] implementations to make
/// sure that brackets are added where necessary. Usually that happens
/// when the object uses an operator with a lower precedence than
/// multiplication, for instance if it is some kind of sum.
///
/// # Example
///
/// ```
/// # use series::NeedsCoeffBracket;
/// // Define a custom complex number type
/// // with integer real and imaginary parts
/// #[derive(Copy, Clone)]
/// struct Complex {
///     re: i32,
///     im: i32,
/// }
///
/// impl NeedsCoeffBracket for Complex {
///     fn needs_coeff_bracket(&self) -> bool {
///         let Self{re, im} = self;
///         // assuming the usual way to display complex numbers as
///         // re + im*i and omitting any vanishing terms
///         // we need brackets in expressions like (re + im*i)*x,
///         // if the following condition is fulfilled:
///         (*re != 0) && (*im != 0)
///     }
/// }
/// ```
pub trait NeedsCoeffBracket {
    fn needs_coeff_bracket(&self) -> bool;
}

impl<T: NeedsCoeffBracket> NeedsCoeffBracket for &T {
    fn needs_coeff_bracket(&self) -> bool {
        (*self).needs_coeff_bracket()
    }
}

macro_rules! impl_never_needs_coeff_bracket {
    ($($t:ty), *) => {
        $(
            impl NeedsCoeffBracket for $t {
                fn needs_coeff_bracket(&self) -> bool {
                    false
                }
            }
        )*
    };
}

impl_never_needs_coeff_bracket!(
    i8, i16, i32, i64, i128, isize,
    u8, u16, u32, u64, u128, usize,
    f32, f64
);
