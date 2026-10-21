use std::{
    cmp::Ordering,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign},
};

use derive_more::{
    Add, AddAssign, AsMut, AsRef, Debug, Display, From, FromStr, Into,
    MulAssign, Neg, Sub, SubAssign,
};
use malachite::{
    Natural,
    base::num::{
        arithmetic::traits::Sign as MalachiteSign,
        basic::traits::{One, Zero},
    },
};

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
pub struct Rational(pub malachite::Rational);

impl<'a> AddAssign<&'a Rational> for Rational {
    fn add_assign(&mut self, rhs: &'a Rational) {
        self.0.add_assign(&rhs.0);
    }
}

impl<'a> SubAssign<&'a Rational> for Rational {
    fn sub_assign(&mut self, rhs: &'a Rational) {
        self.0.sub_assign(&rhs.0);
    }
}

impl<'a> MulAssign<&'a Rational> for Rational {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.mul_assign(&rhs.0);
    }
}

impl<'a> DivAssign<&'a Rational> for Rational {
    fn div_assign(&mut self, rhs: &'a Self) {
        self.0.div_assign(&rhs.0);
    }
}

impl MulAssign for Rational {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.mul_assign(rhs.0);
    }
}

impl DivAssign for Rational {
    fn div_assign(&mut self, rhs: Self) {
        self.0.div_assign(rhs.0);
    }
}

impl Add for &Rational {
    type Output = Rational;

    fn add(self, rhs: Self) -> Self::Output {
        Rational((&self.0).add(&rhs.0))
    }
}

impl Add<&Rational> for Rational {
    type Output = Rational;

    fn add(self, rhs: &Self) -> Self::Output {
        Rational(self.0.add(&rhs.0))
    }
}

impl Add<Rational> for &Rational {
    type Output = Rational;

    fn add(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).add(rhs.0))
    }
}

impl Sub for &Rational {
    type Output = Rational;

    fn sub(self, rhs: Self) -> Self::Output {
        Rational((&self.0).sub(&rhs.0))
    }
}

impl Sub<&Rational> for Rational {
    type Output = Rational;

    fn sub(self, rhs: &Self) -> Self::Output {
        Rational(self.0.sub(&rhs.0))
    }
}

impl Sub<Rational> for &Rational {
    type Output = Rational;

    fn sub(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).sub(rhs.0))
    }
}

impl Mul for Rational {
    type Output = Rational;

    fn mul(self, rhs: Self) -> Self::Output {
        Rational(self.0.mul(rhs.0))
    }
}

impl Mul for &Rational {
    type Output = Rational;

    fn mul(self, rhs: Self) -> Self::Output {
        Rational((&self.0).mul(&rhs.0))
    }
}

impl Mul<&Rational> for Rational {
    type Output = Rational;

    fn mul(self, rhs: &Self) -> Self::Output {
        Rational(self.0.mul(&rhs.0))
    }
}

impl Mul<Rational> for &Rational {
    type Output = Rational;

    fn mul(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).mul(rhs.0))
    }
}

impl Div for Rational {
    type Output = Rational;

    fn div(self, rhs: Self) -> Self::Output {
        Rational(self.0.div(rhs.0))
    }
}

impl Div for &Rational {
    type Output = Rational;

    fn div(self, rhs: Self) -> Self::Output {
        Rational((&self.0).div(&rhs.0))
    }
}

impl Div<&Rational> for Rational {
    type Output = Rational;

    fn div(self, rhs: &Self) -> Self::Output {
        Rational(self.0.div(&rhs.0))
    }
}

impl Div<Rational> for &Rational {
    type Output = Rational;

    fn div(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).div(rhs.0))
    }
}

impl std::ops::Neg for &Rational {
    type Output = Rational;

    fn neg(self) -> Self::Output {
        Rational(-&self.0)
    }
}

impl num_traits::Zero for Rational {
    fn zero() -> Self {
        Self(malachite::Rational::ZERO)
    }

    fn is_zero(&self) -> bool {
        self.0 == malachite::Rational::ZERO
    }
}

impl num_traits::One for Rational {
    fn one() -> Self {
        Self(malachite::Rational::ONE)
    }

    fn is_one(&self) -> bool {
        self.0 == malachite::Rational::ONE
    }
}

macro_rules! impl_mul {
    ($($t:ty), *) => {
        $(
            impl<'a, Var> Mul<$t> for Rational
            where $t: Mul<Rational>
            {
                type Output = <$t as Mul<Rational>>::Output;

                fn mul(self, other: $t) -> Self::Output {
                    other.mul(self)
                }
            }
        )*
    };
}

impl_mul!(
    Polynomial<Var, Rational>, &'a Polynomial<Var, Rational>, PolynomialSlice<'a, Var, Rational>,
    Series<Var, Rational>, &'a Series<Var, Rational>, SeriesSlice<'a, Var, Rational>,
    Laurent<Var, Rational>, &'a Laurent<Var, Rational>
);

macro_rules! impl_mul_anon {
    ($($t:ty), *) => {
        $(
            impl<'a> Mul<$t> for Rational {
                type Output = <$t as Mul<Rational>>::Output;

                fn mul(self, other: $t) -> Self::Output {
                    other.mul(self)
                }
            }
        )*
    };
}
impl_mul_anon!(
    AnonSeries<Rational>,
    &'a AnonSeries<Rational>,
    AnonSeriesSlice<'a, Rational>
);

macro_rules! impl_from {
    ($($t:ty), *) => {
        $(
            impl From<$t> for Rational {
                fn from(val: $t) -> Self {
                    Self(malachite::Rational::from(val))
                }
            }
        )*
    };
}

impl_from!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

#[cfg(feature = "malachite-integer")]
impl_from!(malachite::Integer);

#[cfg(feature = "malachite-integer")]
impl From<(malachite::Integer, malachite::Integer)> for Rational {
    fn from(val: (malachite::Integer, malachite::Integer)) -> Self {
        Self(malachite::Rational::from_integers(val.0, val.1))
    }
}

#[cfg(feature = "malachite-integer")]
impl From<super::Integer> for Rational {
    fn from(value: super::Integer) -> Self {
        Self::from(value.0)
    }
}

#[cfg(feature = "malachite-integer")]
impl From<(super::Integer, super::Integer)> for Rational {
    fn from(value: (super::Integer, super::Integer)) -> Self {
        Self::from((value.0.0, value.1.0))
    }
}

impl NeedsCoeffBracket for Rational {
    fn needs_coeff_bracket(&self) -> bool {
        false
    }
}

impl<'a> SplitSign<'a> for Rational {
    type Signless = SignlessRational<'a>;

    fn split_sign(&'a self) -> (Sign, Self::Signless) {
        let sign = match self.0.sign() {
            Ordering::Less => Sign::Minus,
            Ordering::Equal | Ordering::Greater => Sign::Plus,
        };
        let (num, den) = self.0.numerator_and_denominator_ref();
        (sign, SignlessRational { num, den })
    }
}

pub struct SignlessRational<'a> {
    num: &'a Natural,
    den: &'a Natural,
}

impl<'a> Mul for SignlessRational<'a> {
    type Output = Self;

    fn mul(self, _: Self) -> Self::Output {
        unimplemented!(
            "`Mul` is only implemented to satisfy the trait bounds for `One`."
        )
    }
}

impl<'a> num_traits::One for SignlessRational<'a> {
    fn one() -> Self {
        unimplemented!(
            "`One` is only implemented to satisfy trait bounds. Only the `is_one` function should be used"
        )
    }

    fn is_one(&self) -> bool {
        self.num == &Natural::ONE && self.den == &Natural::ONE
    }
}

impl<'a> PartialEq for SignlessRational<'a> {
    fn eq(&self, _: &Self) -> bool {
        unimplemented!(
            "`PartialEq` is only implemented to satisfy the trait bounds for `One::is_one`."
        )
    }
}

impl<'a> std::fmt::Display for SignlessRational<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self { num, den } = self;
        write!(f, "{num}/{den}")
    }
}
