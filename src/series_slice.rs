use crate::ops::{Exp, Ln, Pow};
use crate::traits::{AsSlice, ExpCoeff, MulInverse};
use crate::{
    Coeff, Iter, PolynomialSlice, Series, anon_series_slice::AnonSeriesSlice,
};

use num_traits::{One, Zero};
use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign,
};

// TODO: lots of code duplication with SeriesSlice

/// View into a Laurent series
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct SeriesSlice<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) series: AnonSeriesSlice<'a, C>,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<Var, C: Coeff> std::marker::Copy for SeriesSlice<'_, Var, C> {}

impl<Var, C: Coeff> std::clone::Clone for SeriesSlice<'_, Var, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, Var, C: Coeff> SeriesSlice<'a, Var, C> {
    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec![1,2,3]);
    /// assert_eq!(s.as_slice(..).min_pow(), -1);
    /// assert_eq!(s.as_slice(0..).min_pow(), 0);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.series.min_pow()
    }

    /// Get the power of the expansion variable where the slice is
    /// truncated.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.as_slice(..).cutoff_pow(), 2);
    /// assert_eq!(s.as_slice(..1).cutoff_pow(), 1);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.series.cutoff_pow()
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..);
    /// let mut iter = slice.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<'_, C> {
        self.series.iter()
    }

    /// Split a series slice into two at the given power
    /// of the expansion variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let (lower, upper) = s.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), -1);
    /// assert_eq!(upper.min_pow(), 0);
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let Self { series, var } = self;
        let (lower, upper) = series.split_at(pos);
        let lower = Self { series: lower, var };
        let upper = Self { series: upper, var };
        (lower, upper)
    }

    /// View as polynomial slice
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..).as_poly();
    /// let p = series::Polynomial::from(s.clone());
    /// assert_eq!(slice, p.as_slice(..));
    /// ```
    pub fn as_poly(self) -> PolynomialSlice<'a, Var, C> {
        let Self { var, series } = self;
        let AnonSeriesSlice { min_pow, coeffs } = series;
        PolynomialSlice::Poly {
            min_pow,
            coeffs,
            var,
        }
    }

    /// Get the expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.as_slice(..).var(), &"x");
    /// ```
    pub fn var(&self) -> &'a Var {
        self.var
    }

    /// Try to get the series coefficient of the expansion variable to
    /// the given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power or below the leading power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..);
    /// assert_eq!(slice.try_coeff(-5), None);
    /// assert_eq!(slice.try_coeff(-2), None);
    /// assert_eq!(slice.try_coeff(-1), Some(&1));
    /// assert_eq!(slice.try_coeff(0), Some(&2));
    /// assert_eq!(slice.try_coeff(1), Some(&3));
    /// assert_eq!(slice.try_coeff(2), None);
    /// assert_eq!(slice.try_coeff(5), None);
    /// ```
    pub fn try_coeff(&self, pow: isize) -> Option<&'a C> {
        self.series.try_coeff(pow)
    }
}

impl<'a, Var, C: 'static + Coeff + Send + Sync> SeriesSlice<'a, Var, C> {
    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..);
    /// assert_eq!(slice.coeff(-5), Some(&0));
    /// assert_eq!(slice.coeff(-2), Some(&0));
    /// assert_eq!(slice.coeff(-1), Some(&1));
    /// assert_eq!(slice.coeff(0), Some(&2));
    /// assert_eq!(slice.coeff(1), Some(&3));
    /// assert_eq!(slice.coeff(2), None);
    /// assert_eq!(slice.coeff(5), None);
    /// ```
    pub fn coeff(&self, pow: isize) -> Option<&'a C> {
        self.series.coeff(pow)
    }
}

impl<Var, C: Coeff> Index<isize> for SeriesSlice<'_, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.series[index]
    }
}

impl<Var, C> MulInverse for SeriesSlice<'_, Var, C>
where
    Var: Clone,
    C: Coeff + SubAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
{
    type Output = Series<Var, C>;

    fn mul_inverse(self) -> Self::Output {
        let inv = self.series.mul_inverse();
        inv.in_var(self.var.clone())
    }
}

impl<Var: Clone, C: Coeff> Neg for SeriesSlice<'_, Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Series<Var, C>;

    fn neg(self) -> Self::Output {
        self.series.neg().in_var(self.var.clone())
    }
}

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for SeriesSlice<'_, Var, C>
where
    Series<Var, C>: AddAssign<Rhs>,
{
    type Output = Series<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = Series::from(self);
        res += other;
        res
    }
}

impl<Var, C: Coeff, T> Sub<T> for SeriesSlice<'_, Var, C>
where
    C: Clone,
    Var: Clone,
    Series<Var, C>: SubAssign<T>,
{
    type Output = Series<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = Series::from(self);
        res -= other;
        res
    }
}

impl<'a, Var, C: Coeff> Mul for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: Mul<SeriesSlice<'a, Var, C>, Output = Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: SeriesSlice<'a, Var, C>) -> Self::Output {
        Series::from(self) * other
    }
}

impl<Var, C: Coeff> Mul<Series<Var, C>> for SeriesSlice<'_, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: MulAssign<Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: Series<Var, C>) -> Self::Output {
        Series::from(self) * other
    }
}

impl<'b, Var, C: Coeff> Mul<&'b Series<Var, C>> for SeriesSlice<'_, Var, C>
where
    C: Clone,
    Var: Clone,
    for<'c> Series<Var, C>:
        Mul<SeriesSlice<'c, Var, C>, Output = Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: &'b Series<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<Var, C: Coeff> Mul<C> for SeriesSlice<'_, Var, C>
where
    Var: Clone,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: C) -> Self::Output {
        (self.series * other).in_var(self.var.clone())
    }
}

impl<'b, Var, C: Coeff> Mul<&'b C> for SeriesSlice<'_, Var, C>
where
    Var: Clone,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: &'b C) -> Self::Output {
        (self.series * other).in_var(self.var.clone())
    }
}

impl<Var, C: Coeff, T> Div<T> for SeriesSlice<'_, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: DivAssign<T>,
{
    type Output = Series<Var, C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = Series::from(self);
        res /= other;
        res
    }
}

impl<Var, C: Coeff> Exp for SeriesSlice<'_, Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
    for<'b> C: MulAssign<&'b C>,
    Var: Clone,
    C: Clone
        + Div<Output = C>
        + Mul<Output = C>
        + AddAssign
        + Exp<Output = C>
        + From<i32>,
{
    type Output = Series<Var, C>;

    fn exp(self) -> Self::Output {
        Series::new(self.var.clone(), 0, self.exp_coeff())
    }
}

impl<Var, C: Coeff> Ln for SeriesSlice<'_, Var, C>
where
    for<'b> C: Div<&'b C, Output = C>,
    for<'b> &'b C: Mul<Output = C> + Ln<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
        + From<i32>
        + From<Var>,
    Var: Clone,
{
    type Output = Series<Var, C>;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(self) -> Self::Output {
        assert!(!self.series.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = &self.series.coeffs[0];
        // self.coeffs[0] = C::one();
        // for i in 1..self.coeffs.len() {
        //     self.coeffs[i] /= &c_k0;
        // }
        let a = &self.series.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = if k0 != 0 {
            let var = self.var.clone();
            c_k0.ln() + C::from(k0 as i32) * C::from(var).ln()
        } else {
            c_k0.ln()
        };
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone() / c_k0);
            for i in 1..n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                let tmp = num_factor * (&a[n - i] * &b[i]) / c_k0;
                b[n] -= tmp;
            }
        }
        Series::new(self.var.clone(), 0, b)
    }
}

impl<Var, C: Coeff, T> Pow<T> for SeriesSlice<'_, Var, C>
where
    for<'b> SeriesSlice<'b, Var, C>: Ln<Output = Series<Var, C>>,
    Series<Var, C>: Mul<T>,
    <Series<Var, C> as Mul<T>>::Output: Exp,
{
    type Output = <<Series<Var, C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

macro_rules! impl_num_display {
    ($($t:ty), *) => {
        $(
            impl<'a, Var: Display> Display for SeriesSlice<'a, Var, $t> {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    let terms = self.iter().filter(|(_, c)| !c.is_zero());
                    let mut first = true;
                    for (pow, c) in terms {
                        let mut c = *c;
                        if !first {
                            if c < <$t>::zero() {
                                c = -c;
                                write!(f, " - ")?;
                            } else {
                                write!(f, " + ")?;
                            }
                        }
                        first = false;
                        if pow == 0 {
                            write!(f, "{c}")?;
                        } else {
                            if !c.is_one() {
                                write!(f, "{c}*")?;
                            }
                            write!(f, "{}", self.var())?;
                            if pow != 1 {
                                write!(f, "^{pow}")?;
                            }
                        }
                    }
                    if !first {
                        write!(f, " + ")?;
                    }
                    if self.cutoff_pow() == 1 {
                        write!(f, "O({})", self.var())
                    } else {
                        write!(f, "O({}^{})", self.var(), self.cutoff_pow())
                    }
                }
            }
        )*
    };
}

impl_num_display!(i8, i16, i32, i64, i128, isize, f32, f64);

macro_rules! impl_unsigned_display {
    ($($t:ty), *) => {
        $(
            impl<'a, Var: Display> Display for SeriesSlice<'a, Var, $t> {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    let terms = self.iter().filter(|(_, c)| !c.is_zero());
                    let mut first = true;
                    for (pow, c) in terms {
                        if !first {
                            write!(f, " + ")?;
                        }
                        first = false;
                        if pow == 0 {
                            write!(f, "{c}")?;
                        } else {
                            if !c.is_one() {
                                write!(f, "{c}*")?;
                            }
                            write!(f, "{}", self.var())?;
                            if pow != 1 {
                                write!(f, "^{pow}")?;
                            }
                        }
                    }
                    if !first {
                        write!(f, " + ")?;
                    }
                    if self.cutoff_pow() == 1 {
                        write!(f, "O({})", self.var())
                    } else {
                        write!(f, "O({}^{})", self.var(), self.cutoff_pow())
                    }
                }
            }
        )*
    };
}

impl_unsigned_display!(u8, u16, u32, u64, u128, usize);
