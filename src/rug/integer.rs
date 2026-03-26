use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use derive_more::{
    Add, AddAssign, AsMut, AsRef, Debug, Display, From, FromStr, Into,
    MulAssign, Neg, Sub, SubAssign,
};
use num_traits::{One, Zero};
use rug::{Complete, integer::BorrowInteger};

use crate::{
    Laurent, NeedsCoeffBracket, Polynomial, PolynomialSlice, Series,
    SeriesSlice, Sign, SplitSign, anon_series::AnonSeries,
    anon_series_slice::AnonSeriesSlice,
};

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(
    Clone,
    Default,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Add,
    AddAssign,
    AsMut,
    AsRef,
    Debug,
    Display,
    From,
    FromStr,
    Into,
    MulAssign,
    Neg,
    Sub,
    SubAssign,
)]
pub struct Integer(pub rug::Integer);

impl Integer {
    pub const ZERO: Integer = Integer(rug::Integer::ZERO);
}

impl<'a> AddAssign<&'a Integer> for Integer {
    fn add_assign(&mut self, rhs: &'a Integer) {
        self.0.add_assign(&rhs.0);
    }
}

impl<'a> SubAssign<&'a Integer> for Integer {
    fn sub_assign(&mut self, rhs: &'a Integer) {
        self.0.sub_assign(&rhs.0);
    }
}

impl<'a> MulAssign<&'a Integer> for Integer {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.mul_assign(&rhs.0);
    }
}

impl MulAssign for Integer {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.mul_assign(&rhs.0);
    }
}

impl Add for &Integer {
    type Output = Integer;

    fn add(self, rhs: Self) -> Self::Output {
        Integer((&self.0).add(&rhs.0).complete())
    }
}

impl Add<&Integer> for Integer {
    type Output = Integer;

    fn add(self, rhs: &Self) -> Self::Output {
        Integer(self.0.add(&rhs.0))
    }
}

impl Add<Integer> for &Integer {
    type Output = Integer;

    fn add(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).add(rhs.0))
    }
}

impl Sub for &Integer {
    type Output = Integer;

    fn sub(self, rhs: Self) -> Self::Output {
        Integer((&self.0).sub(&rhs.0).complete())
    }
}

impl Sub<&Integer> for Integer {
    type Output = Integer;

    fn sub(self, rhs: &Self) -> Self::Output {
        Integer(self.0.sub(&rhs.0))
    }
}

impl Sub<Integer> for &Integer {
    type Output = Integer;

    fn sub(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).sub(rhs.0))
    }
}

impl Mul for Integer {
    type Output = Integer;

    fn mul(self, rhs: Self) -> Self::Output {
        Integer(self.0.mul(rhs.0))
    }
}

impl Mul for &Integer {
    type Output = Integer;

    fn mul(self, rhs: Self) -> Self::Output {
        Integer((&self.0).mul(&rhs.0).complete())
    }
}

impl Mul<&Integer> for Integer {
    type Output = Integer;

    fn mul(self, rhs: &Self) -> Self::Output {
        Integer(self.0.mul(&rhs.0))
    }
}

impl Mul<Integer> for &Integer {
    type Output = Integer;

    fn mul(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).mul(rhs.0))
    }
}

impl std::ops::Neg for &Integer {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer((-&self.0).complete())
    }
}

impl Zero for Integer {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for Integer {
    fn one() -> Self {
        Self(rug::Integer::one())
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

macro_rules! impl_mul {
    ($($t:ty), *) => {
        $(
            impl<'a, Var> Mul<$t> for Integer
            where $t: Mul<Integer>
            {
                type Output = <$t as Mul<Integer>>::Output;

                fn mul(self, other: $t) -> Self::Output {
                    other.mul(self)
                }
            }
        )*
    };
}

impl_mul!(
    Polynomial<Var, Integer>, &'a Polynomial<Var, Integer>, PolynomialSlice<'a, Var, Integer>,
    Series<Var, Integer>, &'a Series<Var, Integer>, SeriesSlice<'a, Var, Integer>,
    Laurent<Var, Integer>, &'a Laurent<Var, Integer>
);

macro_rules! impl_mul_anon {
    ($($t:ty), *) => {
        $(
            impl<'a> Mul<$t> for Integer {
                type Output = <$t as Mul<Integer>>::Output;

                fn mul(self, other: $t) -> Self::Output {
                    other.mul(self)
                }
            }
        )*
    };
}
impl_mul_anon!(
    AnonSeries<Integer>,
    &'a AnonSeries<Integer>,
    AnonSeriesSlice<'a, Integer>
);

macro_rules! impl_from {
    ($($t:ty), *) => {
        $(
            impl From<$t> for Integer {
                fn from(val: $t) -> Self {
                    Self(rug::Integer::from(val))
                }
            }
        )*
    };
}

impl_from!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

macro_rules! impl_try_from {
    ($($t:ty), *) => {
        $(
            impl TryFrom<Integer> for $t {
                type Error = <$t as TryFrom<rug::Integer>>::Error;

                fn try_from(val: Integer) -> Result<Self, Self::Error> {
                    Self::try_from(val.0)
                }
            }

            impl<'a> TryFrom<&'a Integer> for $t {
                type Error = <$t as TryFrom<&'a rug::Integer>>::Error;

                fn try_from(val: &'a Integer) -> Result<Self, Self::Error> {
                    Self::try_from(&val.0)
                }
            }
        )*
    };
}

impl_try_from!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl NeedsCoeffBracket for Integer {
    fn needs_coeff_bracket(&self) -> bool {
        false
    }
}

impl<'a> SplitSign<'a> for Integer {
    type Signless = SignlessInteger<'a>;

    fn split_sign(&'a self) -> (Sign, Self::Signless) {
        let sign = if self.0.is_negative() {
            Sign::Minus
        } else {
            Sign::Plus
        };
        (sign, SignlessInteger(self.0.as_abs()))
    }
}

pub struct SignlessInteger<'a>(BorrowInteger<'a>);

impl<'a> Mul for SignlessInteger<'a> {
    type Output = Self;

    fn mul(self, _: Self) -> Self::Output {
        unimplemented!(
            "`Mul` is only implemented to satisfy the trait bounds for `One`."
        )
    }
}

impl<'a> One for SignlessInteger<'a> {
    fn one() -> Self {
        unimplemented!(
            "`One` is only implemented to satisfy trait bounds. Only the `is_one` function should be used"
        )
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

impl<'a> PartialEq for SignlessInteger<'a> {
    fn eq(&self, _: &Self) -> bool {
        unimplemented!(
            "`PartialEq` is only implemented to satisfy the trait bounds for `One::is_one`."
        )
    }
}

impl<'a> std::fmt::Display for SignlessInteger<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}
