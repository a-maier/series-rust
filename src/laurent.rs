use std::{
    fmt::{Debug, Display},
    ops::{
        Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

use derive_more::{Display, From, IsVariant};
use num_traits::{One, Zero};

use crate::{
    Coeff, NeedsCoeffBracket, Polynomial, PolynomialSlice, Series, SeriesSlice,
    Sign, SplitSign, poly::SignlessPoly, series::SignlessSeries,
};

/// A Laurent polynomial or series in a single variable
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(
    Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Display, From, IsVariant,
)]
pub enum Laurent<Var, C: Coeff> {
    Polynomial(Polynomial<Var, C>),
    Series(Series<Var, C>),
}

impl<Var, C: Coeff> Laurent<Var, C> {
    pub fn zero() -> Self {
        Polynomial::zero().into()
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Self::Polynomial(p) if p.is_zero() => true,
            _ => false,
        }
    }

    pub fn one() -> Self {
        Polynomial::one().into()
    }

    pub fn is_one(&self) -> bool {
        match self {
            Self::Polynomial(p) if p.is_one() => true,
            _ => false,
        }
    }
}

impl<Var, C: Coeff> Neg for Laurent<Var, C>
where
    C: Neg<Output = C>,
{
    type Output = Laurent<Var, C>;

    fn neg(self) -> Self::Output {
        match self {
            Laurent::Polynomial(p) => (-p).into(),
            Laurent::Series(s) => (-s).into(),
        }
    }
}

impl<'a, Var, C: Coeff> Neg for &'a Laurent<Var, C>
where
    Laurent<Var, C>: Clone + Neg<Output = Laurent<Var, C>>,
{
    type Output = Laurent<Var, C>;

    fn neg(self) -> Self::Output {
        self.clone().neg()
    }
}

impl<Var, C: Coeff> Zero for Laurent<Var, C>
where
    Self: Add<Output = Self>,
{
    fn zero() -> Self {
        Laurent::zero()
    }

    fn is_zero(&self) -> bool {
        Laurent::is_zero(&self)
    }
}

impl<Var, C: Coeff> One for Laurent<Var, C>
where
    Self: Mul<Output = Self>,
{
    fn one() -> Self {
        Laurent::one()
    }

    fn is_one(&self) -> bool {
        Laurent::is_one(&self)
    }
}

impl<Var, C: Coeff> Default for Laurent<Var, C>
where
    Self: Zero,
{
    fn default() -> Self {
        Self::zero()
    }
}

macro_rules! impl_add {
    ($lhs:ty, $rhs:ty) => {
        impl<'a, Var, C: Coeff> Add<$rhs> for $lhs
        where
            for<'c> C: AddAssign<&'c C>,
            C: Clone + AddAssign,
            Var: Clone + Debug + PartialEq,
        {
            type Output = Laurent<Var, C>;

            fn add(self, rhs: $rhs) -> Self::Output {
                match (self, rhs) {
                    (Laurent::Polynomial(p), Laurent::Polynomial(q)) => {
                        p.add(q).into()
                    }
                    (Laurent::Polynomial(p), Laurent::Series(s)) => {
                        s.add(p).into()
                    }
                    (Laurent::Series(s), Laurent::Polynomial(p)) => {
                        p.add(s).into()
                    }
                    (Laurent::Series(s), Laurent::Series(t)) => s.add(t).into(),
                }
            }
        }
    };
}

impl_add!(Laurent<Var, C>, Laurent<Var, C>);
impl_add!(Laurent<Var, C>, &'a Laurent<Var, C>);
impl_add!(&'a Laurent<Var, C>, Laurent<Var, C>);
impl_add!(&'a Laurent<Var, C>, &'a Laurent<Var, C>);

macro_rules! impl_add_foreign {
    ($lhs:ty, $rhs:ty) => {
        impl<'a, Var, C: Coeff> Add<$rhs> for $lhs
        where
            for<'c> C: AddAssign<&'c C>,
            C: Clone + AddAssign,
            Var: Clone + Debug + PartialEq,
        {
            type Output = Laurent<Var, C>;

            fn add(self, rhs: $rhs) -> Self::Output {
                match self {
                    Laurent::Polynomial(p) => p.add(rhs).into(),
                    Laurent::Series(s) => s.add(rhs).into(),
                }
            }
        }
    };
}

impl_add_foreign!(Laurent<Var, C>, Polynomial<Var, C>);
impl_add_foreign!(Laurent<Var, C>, &'a Polynomial<Var, C>);
impl_add_foreign!(Laurent<Var, C>, PolynomialSlice<'a, Var, C>);
impl_add_foreign!(Laurent<Var, C>, Series<Var, C>);
impl_add_foreign!(Laurent<Var, C>, &'a Series<Var, C>);
impl_add_foreign!(Laurent<Var, C>, SeriesSlice<'a, Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, Polynomial<Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, &'a Polynomial<Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, PolynomialSlice<'a, Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, Series<Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, &'a Series<Var, C>);
impl_add_foreign!(&'a Laurent<Var, C>, SeriesSlice<'a, Var, C>);

macro_rules! impl_add_to_foreign {
    ($lhs:ty, $rhs:ty) => {
        impl<'a, Var, C: Coeff> Add<$rhs> for $lhs
        where
            for<'c> C: AddAssign<&'c C>,
            C: Clone + AddAssign,
            Var: Clone + Debug + PartialEq,
        {
            type Output = Laurent<Var, C>;

            fn add(self, rhs: $rhs) -> Self::Output {
                rhs.add(self)
            }
        }
    };
}

impl_add_to_foreign!(Polynomial<Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(Polynomial<Var, C>, &'a Laurent<Var, C>);
impl_add_to_foreign!(&'a Polynomial<Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(&'a Polynomial<Var, C>, &'a Laurent<Var, C>);
impl_add_to_foreign!(PolynomialSlice<'a, Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(PolynomialSlice<'a, Var, C>, &'a Laurent<Var, C>);
impl_add_to_foreign!(Series<Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(Series<Var, C>, &'a Laurent<Var, C>);
impl_add_to_foreign!(&'a Series<Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(&'a Series<Var, C>, &'a Laurent<Var, C>);
impl_add_to_foreign!(SeriesSlice<'a, Var, C>, Laurent<Var, C>);
impl_add_to_foreign!(SeriesSlice<'a, Var, C>, &'a Laurent<Var, C>);

macro_rules! impl_add_assign {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> AddAssign<$rhs> for Laurent<Var, C>
            where Self: Add<$rhs, Output = Self> + Default
            {
                fn add_assign(&mut self, rhs: $rhs) {
                    *self = std::mem::take(self).add(rhs);
                }
            }
        )*
    }
}

impl_add_assign!(
    Laurent<Var, C>,
    &'a Laurent<Var, C>,
    Polynomial<Var, C>,
    &'a Polynomial<Var, C>,
    PolynomialSlice<'a, Var, C>,
    Series<Var, C>,
    &'a Series<Var, C>,
    SeriesSlice<'a, Var, C>
);

macro_rules! impl_sub {
    ($lhs:ty, $rhs:ty) => {
        impl<'a, Var, C: Coeff> Sub<$rhs> for $lhs
        where
            for<'c> C: AddAssign<&'c C>,
            for<'c> &'c C: Neg<Output = C>,
            C: Clone + AddAssign + Neg<Output = C>,
            Var: Clone + Debug + PartialEq,
        {
            type Output = Laurent<Var, C>;

            fn sub(self, rhs: $rhs) -> Self::Output {
                match (self, rhs) {
                    (Laurent::Polynomial(p), Laurent::Polynomial(q)) => {
                        p.sub(q).into()
                    }
                    (Laurent::Polynomial(p), Laurent::Series(s)) => {
                        s.sub(p).into()
                    }
                    (Laurent::Series(s), Laurent::Polynomial(p)) => {
                        p.sub(s).into()
                    }
                    (Laurent::Series(s), Laurent::Series(t)) => s.sub(t).into(),
                }
            }
        }
    };
}

impl_sub!(Laurent<Var, C>, Laurent<Var, C>);
impl_sub!(Laurent<Var, C>, &'a Laurent<Var, C>);
impl_sub!(&'a Laurent<Var, C>, Laurent<Var, C>);
impl_sub!(&'a Laurent<Var, C>, &'a Laurent<Var, C>);

macro_rules! impl_sub_foreign {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<$rhs> for Laurent<Var, C>
            where
                Polynomial<Var, C>: Sub<$rhs>,
            <Polynomial<Var, C> as Sub<$rhs>>::Output: Into<Self>,
            Series<Var, C>: Sub<$rhs>,
            <Series<Var, C> as Sub<$rhs>>::Output: Into<Self>,
            {
                type Output = Laurent<Var, C>;

                fn sub(self, rhs: $rhs) -> Self::Output {
                    match self {
                        Laurent::Polynomial(p) => p.sub(rhs).into(),
                        Laurent::Series(s) => s.sub(rhs).into(),
                    }
                }
            }
        )*
    };
}

impl_sub_foreign!(Polynomial<Var, C>, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, Series<Var, C>, &'a Series<Var, C>, SeriesSlice<'a, Var, C>);

macro_rules! impl_sub_ref_foreign {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<$rhs> for &'a Laurent<Var, C>
            where
                &'a Polynomial<Var, C>: Sub<$rhs>,
            <&'a Polynomial<Var, C> as Sub<$rhs>>::Output: Into<Laurent<Var, C>>,
            &'a Series<Var, C>: Sub<$rhs>,
            <&'a Series<Var, C> as Sub<$rhs>>::Output: Into<Laurent<Var, C>>,
            {
                type Output = Laurent<Var, C>;

                fn sub(self, rhs: $rhs) -> Self::Output {
                    match self {
                        Laurent::Polynomial(p) => p.sub(rhs).into(),
                        Laurent::Series(s) => s.sub(rhs).into(),
                    }
                }
            }
        )*
    };
}

impl_sub_ref_foreign!(Polynomial<Var, C>, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, Series<Var, C>, &'a Series<Var, C>, SeriesSlice<'a, Var, C>);

macro_rules! impl_sub_from_foreign {
    ($($lhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<Laurent<Var, C>> for $lhs
            where
                Self: Sub<Polynomial<Var, C>> + Sub<Series<Var, C>>,
            <Self as Sub<Polynomial<Var, C>>>::Output: Into<Laurent<Var, C>>,
            <Self as Sub<Series<Var, C>>>::Output: Into<Laurent<Var, C>>,
            {
                type Output = Laurent<Var, C>;

                fn sub(self, rhs: Laurent<Var, C>) -> Self::Output {
                    match rhs {
                        Laurent::Polynomial(p) => self.sub(p).into(),
                        Laurent::Series(s) => self.sub(s).into(),
                    }
                }
            }
        )*
    };
}

impl_sub_from_foreign!(Polynomial<Var, C>, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, Series<Var, C>, &'a Series<Var, C>, SeriesSlice<'a, Var, C>);

macro_rules! impl_sub_ref_from_foreign {
    ($($lhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<&'a Laurent<Var, C>> for $lhs
            where
                Self: Sub<&'a Polynomial<Var, C>> + Sub<&'a Series<Var, C>>,
            <Self as Sub<&'a Polynomial<Var, C>>>::Output: Into<Laurent<Var, C>>,
            <Self as Sub<&'a Series<Var, C>>>::Output: Into<Laurent<Var, C>>,
            {
                type Output = Laurent<Var, C>;

                fn sub(self, rhs: &'a Laurent<Var, C>) -> Self::Output {
                    match rhs {
                        Laurent::Polynomial(p) => self.sub(p).into(),
                        Laurent::Series(s) => self.sub(s).into(),
                    }
                }
            }
        )*
    };
}

impl_sub_ref_from_foreign!(Polynomial<Var, C>, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, Series<Var, C>, &'a Series<Var, C>, SeriesSlice<'a, Var, C>);

macro_rules! impl_sub_assign {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> SubAssign<$rhs> for Laurent<Var, C>
            where Self: Sub<$rhs, Output = Self> + Default
            {
                fn sub_assign(&mut self, rhs: $rhs) {
                    *self = std::mem::take(self).sub(rhs);
                }
            }
        )*
    }
}

impl_sub_assign!(
    Laurent<Var, C>,
    &'a Laurent<Var, C>,
    Polynomial<Var, C>,
    &'a Polynomial<Var, C>,
    PolynomialSlice<'a, Var, C>,
    Series<Var, C>,
    &'a Series<Var, C>,
    SeriesSlice<'a, Var, C>
);

impl<Var, C: Coeff> Mul for Laurent<Var, C>
where
    Polynomial<Var, C>: Mul<Output = Polynomial<Var, C>>
        + Mul<Series<Var, C>, Output = Series<Var, C>>,
    Series<Var, C>: Mul<Output = Series<Var, C>>
        + Mul<Polynomial<Var, C>, Output = Series<Var, C>>,
    Self: Zero,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Laurent::Polynomial(p1), Laurent::Polynomial(p2)) => {
                p1.mul(p2).into()
            }
            (Laurent::Polynomial(p), Laurent::Series(s))
            | (Laurent::Series(s), Laurent::Polynomial(p)) => {
                if p.is_zero() {
                    Self::zero()
                } else {
                    s.mul(p).into()
                }
            }
            (Laurent::Series(s1), Laurent::Series(s2)) => s1.mul(s2).into(),
        }
    }
}

impl<'a, Var, C: Coeff> Mul<&'a Laurent<Var, C>> for Laurent<Var, C>
where
    Polynomial<Var, C>: Mul<&'a Polynomial<Var, C>, Output = Polynomial<Var, C>>
        + Mul<&'a Series<Var, C>, Output = Series<Var, C>>,
    Series<Var, C>: Mul<&'a Series<Var, C>, Output = Series<Var, C>>
        + Mul<&'a Polynomial<Var, C>, Output = Series<Var, C>>,
    Self: Zero,
{
    type Output = Self;

    fn mul(self, rhs: &'a Laurent<Var, C>) -> Self::Output {
        match (self, rhs) {
            (Laurent::Polynomial(p1), Laurent::Polynomial(p2)) => {
                p1.mul(p2).into()
            }
            (Laurent::Polynomial(p), Laurent::Series(s)) => {
                if p.is_zero() {
                    Self::zero()
                } else {
                    p.mul(s).into()
                }
            }
            (Laurent::Series(s), Laurent::Polynomial(p)) => {
                if p.is_zero() {
                    Self::zero()
                } else {
                    s.mul(p).into()
                }
            }
            (Laurent::Series(s1), Laurent::Series(s2)) => s1.mul(s2).into(),
        }
    }
}

impl<'a, Var, C: Coeff> Mul<Laurent<Var, C>> for &'a Laurent<Var, C>
where
    Laurent<Var, C>: Mul<&'a Laurent<Var, C>, Output = Laurent<Var, C>>,
{
    type Output = Laurent<Var, C>;

    fn mul(self, rhs: Laurent<Var, C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, Var, C: Coeff> Mul for &'a Laurent<Var, C>
where
    &'a Polynomial<Var, C>: Mul<Output = Polynomial<Var, C>>
        + Mul<&'a Series<Var, C>, Output = Series<Var, C>>,
    &'a Series<Var, C>: Mul<Output = Series<Var, C>>
        + Mul<&'a Polynomial<Var, C>, Output = Series<Var, C>>,
    Laurent<Var, C>: Zero,
{
    type Output = Laurent<Var, C>;

    fn mul(self, rhs: &'a Laurent<Var, C>) -> Self::Output {
        match (self, rhs) {
            (Laurent::Polynomial(p1), Laurent::Polynomial(p2)) => {
                p1.mul(p2).into()
            }
            (Laurent::Polynomial(p), Laurent::Series(s))
            | (Laurent::Series(s), Laurent::Polynomial(p)) => {
                if p.is_zero() {
                    Self::Output::zero()
                } else {
                    s.mul(p).into()
                }
            }
            (Laurent::Series(s1), Laurent::Series(s2)) => s1.mul(s2).into(),
        }
    }
}

impl<Var, C: Coeff> Div for Laurent<Var, C>
where
    Polynomial<Var, C>: Div<C, Output = Polynomial<Var, C>>,
    Series<Var, C>: Div<Output = Series<Var, C>>
        + Div<Polynomial<Var, C>, Output = Series<Var, C>>
        + Mul<C, Output = Series<Var, C>>,
    Laurent<Var, C>: Zero,
    C: Div<Output = C>,
{
    type Output = Self;

    /// Divide two Laurent series or polynomials
    ///
    /// # Panics
    ///
    /// Panics if the dividend is a polynomial and the divisor is a
    /// non-constant polynomial.
    fn div(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (
                Laurent::Polynomial(p),
                Laurent::Polynomial(Polynomial::Const(c)),
            ) => p.div(c).into(),
            (
                Laurent::Polynomial(_),
                Laurent::Polynomial(Polynomial::Poly(_)),
            ) => unimplemented!(
                "Dividing a polynomial by a non-constant polynomial"
            ),
            (Laurent::Polynomial(Polynomial::Const(c)), Laurent::Series(s)) => {
                if c.is_zero() {
                    Self::zero()
                } else {
                    // we cannot implement Div<Series> for C due to orphan rules
                    let c_inv = C::one().div(c);
                    s.mul(c_inv).into()
                }
            }
            (Laurent::Polynomial(Polynomial::Poly(p)), Laurent::Series(s)) => {
                let cutoff_pow = s.len() as isize + p.min_pow();
                p.cutoff_at(cutoff_pow).div(s).into()
            }
            (Laurent::Series(s), Laurent::Polynomial(p)) => s.div(p).into(),
            (Laurent::Series(s1), Laurent::Series(s2)) => s1.div(s2).into(),
        }
    }
}

macro_rules! impl_add_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> AddAssign<$rhs> for Laurent<Var, C>
            where C: AddAssign<$rhs>
            {
                fn add_assign(&mut self, rhs: $rhs) {
                    match self {
                        Laurent::Polynomial(p) => p.add_assign(rhs),
                        Laurent::Series(s) => s.add_assign(rhs),
                    }
                }
            }
        )*
    }
}

impl_add_assign_const!(C, &'a C);

macro_rules! impl_sub_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> SubAssign<$rhs> for Laurent<Var, C>
            where C: SubAssign<$rhs>
            {
                fn sub_assign(&mut self, rhs: $rhs) {
                    match self {
                        Laurent::Polynomial(p) => p.sub_assign(rhs),
                        Laurent::Series(s) => s.sub_assign(rhs),
                    }
                }
            }
        )*
    }
}

impl_sub_assign_const!(C, &'a C);

macro_rules! impl_mul_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> MulAssign<$rhs> for Laurent<Var, C>
            where for<'c> C: MulAssign<&'c C>
            {
                fn mul_assign(&mut self, rhs: $rhs) {
                    match self {
                        Laurent::Polynomial(p) => p.mul_assign(rhs),
                        Laurent::Series(s) => if rhs.is_zero() {
                            *self = Self::zero();
                        } else {
                            s.mul_assign(rhs)
                        },
                    }
                }
            }
        )*
    }
}

impl_mul_assign_const!(C, &'a C);

macro_rules! impl_div_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> DivAssign<$rhs> for Laurent<Var, C>
            where for<'c> C: DivAssign<&'c C>
            {
                fn div_assign(&mut self, rhs: $rhs) {
                    match self {
                        Laurent::Polynomial(p) => p.div_assign(rhs),
                        Laurent::Series(s) => s.div_assign(rhs),
                    }
                }
            }
        )*
    }
}

impl_div_assign_const!(C, &'a C);

macro_rules! impl_add_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Add<$rhs> for Laurent<Var, C>
            where Self: AddAssign<$rhs>
            {
                type Output = Self;

                fn add(mut self, rhs: $rhs) -> Self::Output{
                    self.add_assign(rhs);
                    self
                }
            }

            impl<'a, Var, C: Coeff> Add<$rhs> for &'a Laurent<Var, C>
            where Laurent<Var, C>: Clone + AddAssign<$rhs>
            {
                type Output = Laurent<Var, C>;

                fn add(self, rhs: $rhs) -> Self::Output{
                    self.clone().add(rhs)
                }
            }
        )*
    }
}

impl_add_const!(C, &'a C);

macro_rules! impl_sub_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<$rhs> for Laurent<Var, C>
            where Self: SubAssign<$rhs>
            {
                type Output = Self;

                fn sub(mut self, rhs: $rhs) -> Self::Output{
                    self.sub_assign(rhs);
                    self
                }
            }

            impl<'a, Var, C: Coeff> Sub<$rhs> for &'a Laurent<Var, C>
            where Laurent<Var, C>: Clone + SubAssign<$rhs>
            {
                type Output = Laurent<Var, C>;

                fn sub(self, rhs: $rhs) -> Self::Output{
                    self.clone().sub(rhs)
                }
            }
        )*
    }
}

impl_sub_const!(C, &'a C);

macro_rules! impl_mul_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Mul<$rhs> for Laurent<Var, C>
            where Self: MulAssign<$rhs>
            {
                type Output = Self;

                fn mul(mut self, rhs: $rhs) -> Self::Output{
                    self.mul_assign(rhs);
                    self
                }
            }

            impl<'a, Var, C: Coeff> Mul<$rhs> for &'a Laurent<Var, C>
            where Laurent<Var, C>: Clone + MulAssign<$rhs>
            {
                type Output = Laurent<Var, C>;

                fn mul(self, rhs: $rhs) -> Self::Output{
                    self.clone().mul(rhs)
                }
            }
        )*
    }
}

impl_mul_const!(C, &'a C);

macro_rules! impl_div_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Div<$rhs> for Laurent<Var, C>
            where Self: DivAssign<$rhs>
            {
                type Output = Self;

                fn div(mut self, rhs: $rhs) -> Self::Output{
                    self.div_assign(rhs);
                    self
                }
            }

            impl<'a, Var, C: Coeff> Div<$rhs> for &'a Laurent<Var, C>
            where Laurent<Var, C>: Clone + DivAssign<$rhs>
            {
                type Output = Laurent<Var, C>;

                fn div(self, rhs: $rhs) -> Self::Output{
                    self.clone().div(rhs)
                }
            }
        )*
    }
}

impl_div_const!(C, &'a C);

impl<Var, C: Coeff> NeedsCoeffBracket for Laurent<Var, C>
where
    Polynomial<Var, C>: NeedsCoeffBracket,
    Series<Var, C>: NeedsCoeffBracket,
{
    fn needs_coeff_bracket(&self) -> bool {
        match self {
            Laurent::Polynomial(p) => p.needs_coeff_bracket(),
            Laurent::Series(s) => s.needs_coeff_bracket(),
        }
    }
}

impl<'a, Var: 'a, C: Coeff + 'a> SplitSign<'a> for Laurent<Var, C>
where
    Polynomial<Var, C>: SplitSign<'a>,
    Series<Var, C>: SplitSign<'a>,
{
    type Signless = SignlessLaurent<'a, Var, C>;

    fn split_sign(&'a self) -> (Sign, Self::Signless) {
        let sign = match self {
            Laurent::Polynomial(p) => p.split_sign().0,
            Laurent::Series(s) => s.split_sign().0,
        };
        (sign, SignlessLaurent(self))
    }
}

#[derive(PartialEq)]
pub struct SignlessLaurent<'a, Var, C: Coeff>(&'a Laurent<Var, C>);

impl<'a, C: Coeff, Var: Display> Display for SignlessLaurent<'a, Var, C>
where
    SignlessPoly<'a, Var, C>: Display,
    SignlessSeries<'a, Var, C>: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Laurent::Polynomial(p) => SignlessPoly(p).fmt(f),
            Laurent::Series(s) => SignlessSeries(s).fmt(f),
        }
    }
}

impl<'a, Var, C: Coeff> Mul for SignlessLaurent<'a, Var, C> {
    type Output = Self;

    fn mul(self, _: Self) -> Self::Output {
        unimplemented!(
            "`Mul` is only implemented to satisfy the trait bounds for `One`."
        )
    }
}

impl<'a, Var: PartialEq, C: Coeff + SplitSign<'a>> One
    for SignlessLaurent<'a, Var, C>
where
    <C as SplitSign<'a>>::Signless: One + PartialEq,
{
    fn one() -> Self {
        unimplemented!(
            "`One` is only implemented to satisfy trait bounds. Only the `is_one` function should be used"
        )
    }

    fn is_one(&self) -> bool {
        match self.0 {
            Laurent::Polynomial(p) => p.split_sign().1.is_one(),
            Laurent::Series(_) => false,
        }
    }
}
