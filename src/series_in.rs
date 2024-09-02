use crate::ops::{Exp, Ln, Pow};
use crate::{series_slice_in::*, Series};
use crate::traits::*;
use crate::{Coeff, IntoIter, Iter};

use std::convert::From;
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};

/// Laurent series in a single variable up to some power
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct SeriesIn<Var, C: Coeff> {
    pub(crate) series: Series<C>,
    pub(crate) var: Var,
}

impl<Var, C: Coeff> SeriesIn<Var, C> {
    /// Create a new series
    ///
    /// # Example
    ///
    /// This creates a series in the variable "x", starting at "x"^-1
    /// with coefficients 1, 2, 3. In other words, the series x^-1 + 2 +
    /// 3*x + O(x^2).
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// ```
    pub fn new(var: Var, min_pow: isize, coeffs: Vec<C>) -> SeriesIn<Var, C> {
        let series = Series::new(min_pow, coeffs);
        SeriesIn {
            series,
            var,
        }
    }

    /// Create a new series with a given cutoff power
    ///
    /// # Example
    ///
    /// This creates a series in the variable "x", starting at "x"^-1
    /// with coefficients 1, 2, 3. and vanishing coefficients up to
    /// "x"^5 .In other words, the series
    /// x^-1 + 2 + 3*x + O(x^5).
    /// ```rust
    /// let s = series::SeriesIn::with_cutoff("x", -1, 5, vec!(1,2,3));
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the cutoff power is lower than the starting power
    ///
    pub fn with_cutoff(
        var: Var,
        min_pow: isize,
        cutoff_pow: isize,
        coeffs: Vec<C>,
    ) -> SeriesIn<Var, C> {
        let series = Series::with_cutoff(min_pow..cutoff_pow, coeffs);
        SeriesIn {
            series,
            var,
        }
    }

    /// Get the expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        &self.var
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.min_pow(), -1);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.series.min_pow()
    }

    /// Get the power of the expansion variable where the series is
    /// truncated.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.cutoff_pow(), 2);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.series.cutoff_pow()
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// let mut iter = s.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        self.series.iter()
    }

    /// Try to get the series coefficient of the expansion variable to
    /// the given power.
    ///
    /// Returns [None] if the requested power is above the highest known
    /// power or below the leading power.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.try_coeff(-5), None);
    /// assert_eq!(s.try_coeff(-2), None);
    /// assert_eq!(s.try_coeff(-1), Some(&1));
    /// assert_eq!(s.try_coeff(0), Some(&2));
    /// assert_eq!(s.try_coeff(1), Some(&3));
    /// assert_eq!(s.try_coeff(2), None);
    /// assert_eq!(s.try_coeff(5), None);
    /// ```
    pub fn try_coeff(&self, pow: isize) -> Option<&C> {
        self.series.try_coeff(pow)
    }

    /// Apply a function to a specific coefficient
    ///
    /// `f(c)` is applied to the coefficient `c` of the variable to
    /// the power `pow`.
    ///
    /// # Panics
    ///
    /// Panics if `pow` is equal to or larger than
    /// [cutoff_pow](Self::cutoff_pow).
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut s = series::SeriesIn::new("x", -1, vec![1,2,3]);
    /// s.apply_at(0, |c| *c = 0);
    /// assert_eq!(s.coeff(0), Some(&0));
    ///
    /// // We can remove existing terms and add new ones, provided the
    /// // variable power is less than `s.cutoff_pow()`!
    /// s.apply_at(-1, |c| *c = 0);
    /// assert_eq!(s.min_pow(), 1);
    /// s.apply_at(-3, |c| *c = 1);
    /// assert_eq!(s.min_pow(), -3);
    /// assert_eq!(s.coeff(-3), Some(&1));
    /// ```
    pub fn apply_at<F: FnOnce(&mut C)>(&mut self, pow: isize, f: F) {
        self.series.apply_at(pow, f)
    }

    /// Transform all coefficients
    ///
    /// `f(p, c)` is applied to each term, where `p` is the power of
    /// the variable and `c` a mutable reference to the
    /// coefficient. `p` takes all values in the range
    /// `min_pow()..cutoff_pow()`.
    ///
    /// # Example
    ///
    /// Replace each coefficient by its square
    /// ```rust
    /// let mut s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// s.for_each(|_, c| *c *= *c);
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&4));
    /// assert_eq!(s.coeff(1), Some(&9));
    /// assert_eq!(s.coeff(2), Some(&16));
    /// ```
    pub fn for_each<F>(&mut self, f: F)
    where
        F: FnMut(isize, &mut C)
    {
        self.series.for_each(f)
    }
}

impl<Var, C: 'static + Coeff + Send + Sync> SeriesIn<Var, C> {
    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.coeff(-5), Some(&0));
    /// assert_eq!(s.coeff(-2), Some(&0));
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&2));
    /// assert_eq!(s.coeff(1), Some(&3));
    /// assert_eq!(s.coeff(2), None);
    /// assert_eq!(s.coeff(5), None);
    /// ```
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        self.series.coeff(pow)
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, Range<isize>>
    for SeriesIn<Var, C>
{
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the lower bound is smaller than the leading power
    /// or the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..2);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeInclusive<isize>>
    for SeriesIn<Var, C>
{
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the lower bound is smaller than the leading power
    /// or the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..=1);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>>
    for SeriesIn<Var, C>
{
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..=1);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>>
    for SeriesIn<Var, C>
{
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the lower bound is smaller than the leading power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert!(std::panic::catch_unwind(|| t[-1]).is_err());
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>>
    for SeriesIn<Var, C>
{
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..2);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFull> for SeriesIn<Var, C> {
    type Output = SeriesSliceIn<'a, Var, C>;

    /// A slice containing the complete series.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert_eq!(t[-1], s[-1]);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        self.series.as_slice(r).in_var(self.var())
    }
}

impl<Var, C: Coeff> Index<isize> for SeriesIn<Var, C> {
    type Output = C;

    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// # Panics
    ///
    /// Panics if the index is smaller than the leading power or
    /// at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s[-1], 1);
    /// assert_eq!(s[0], 2);
    /// assert_eq!(s[1], 3);
    /// assert!(std::panic::catch_unwind(|| s[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| s[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        &self.series[index]
    }
}

impl<Var, C: Coeff> std::iter::IntoIterator for SeriesIn<Var, C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the series coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// let mut iter = s.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        self.series.into_iter()
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> MulInverse for &'a SeriesIn<Var, C>
where
    Var: Clone,
    C: Coeff + SubAssign,
    for<'c> &'c C: Div<Output = C> + Mul<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::SeriesIn::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = (&s).mul_inverse();
    /// let one = series::SeriesIn::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self.series).mul_inverse().in_var(self.var.clone())
    }
}

impl<Var: Clone, C: Coeff + SubAssign> MulInverse for SeriesIn<Var, C>
where
    for<'a> &'a SeriesIn<Var, C>: MulInverse<Output = SeriesIn<Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::SeriesIn::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = s.clone().mul_inverse();
    /// let one = series::SeriesIn::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self).mul_inverse()
    }
}

impl<Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for SeriesIn<Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.as_slice(..).fmt(f)
    }
}

impl<Var, C: Coeff + Neg<Output = C>> Neg for SeriesIn<Var, C> {
    type Output = SeriesIn<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::SeriesIn::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let Self {
            series,
            var,
        } = self;
        let series = -series;
        Self {
            series,
            var,
        }
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a SeriesIn<Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::SeriesIn::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-&s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<&'a SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let t = SeriesIn::new("x", -1, vec!(3., 4., 5.));
    /// let res = SeriesIn::new("x", -3, vec!(1.,0.,0.));
    /// s += &t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(&mut self, other: &'a SeriesIn<Var, C>) {
        assert_eq!(self.var(), other.var());
        self.series.add_assign(&other.series)
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<SeriesSliceIn<'a, Var, C>> for SeriesIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: SeriesSliceIn<'a, Var, C>) {
        assert_eq!(self.var(), other.var());
        self.series.add_assign(other.series)
    }
}

impl<Var, C: Coeff> AddAssign<SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
    Var: PartialEq + fmt::Debug,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let t = SeriesIn::new("x", -1, vec!(3., 4., 5.));
    /// let res = SeriesIn::new("x", -3, vec!(1.,0.,0.));
    /// s += t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(&mut self, other: SeriesIn<Var, C>) {
        assert_eq!(self.var(), other.var());
        self.series.add_assign(other.series)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for &'a SeriesIn<Var, C>
where
    SeriesIn<Var, C>: AddAssign<Rhs>,
{
    type Output = SeriesIn<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: AddAssign<Rhs>,
{
    type Output = SeriesIn<Var, C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, Var, C: Coeff> SubAssign<&'a SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    for<'c> &'c SeriesIn<Var, C>: Neg<Output = SeriesIn<Var, C>>,
    SeriesIn<Var, C>: AddAssign<SeriesIn<Var, C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let res = SeriesIn::new("x", 0, vec!());
    /// s -= &s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(&mut self, other: &'a SeriesIn<Var, C>) {
        *self += -other;
    }
}

impl<'a, Var, C: Coeff> SubAssign<SeriesSliceIn<'a, Var, C>>
    for SeriesIn<Var, C>
where
    for<'c> SeriesSliceIn<'c, Var, C>: Neg<Output = SeriesIn<Var, C>>,
    SeriesIn<Var, C>: AddAssign<SeriesIn<Var, C>>,
{
    fn sub_assign(&mut self, other: SeriesSliceIn<'a, Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff> SubAssign<SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: AddAssign + Neg<Output = SeriesIn<Var, C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// let res = SeriesIn::new("x", 0, vec!());
    /// s -= s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(&mut self, other: SeriesIn<Var, C>) {
        *self += -other;
    }
}

// TODO: somehow make addition symmetric?
impl<'a, Var, C: Coeff, T> Sub<T> for &'a SeriesIn<Var, C>
where
    SeriesIn<Var, C>: Clone + SubAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res -= other;
        res
    }
}

impl<Var, C: Coeff, T> Sub<T> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: SubAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone + AddAssign>
    MulAssign<&'a SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = SeriesIn::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: &'a SeriesIn<Var, C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, Var, C> MulAssign<SeriesSliceIn<'a, Var, C>> for SeriesIn<Var, C>
where
    Var: PartialEq + fmt::Debug,
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C> + Coeff + Clone + AddAssign,
{
    fn mul_assign(&mut self, other: SeriesSliceIn<'a, Var, C>) {
        assert_eq!(self.var(), other.var());
        self.series.mul_assign(other.series);
    }
}

impl<Var, C: Coeff> MulAssign for SeriesIn<Var, C>
where
    for<'a> SeriesIn<Var, C>: MulAssign<&'a SeriesIn<Var, C>>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = SeriesIn::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: SeriesIn<Var, C>) {
        *self *= &other
    }
}

// TODO: somehow make multiplication symmetric?
impl<Var, C: Coeff> Mul for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign,
{
    type Output = SeriesIn<Var, C>;

    fn mul(mut self, other: SeriesIn<Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a SeriesIn<Var, C>> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign<SeriesSliceIn<'a, Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: &'a SeriesIn<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var, C: Coeff> Mul<SeriesSliceIn<'a, Var, C>> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign<SeriesSliceIn<'a, Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(mut self, other: SeriesSliceIn<'a, Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<Var, C: Coeff> Mul<C> for SeriesIn<Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a C> for SeriesIn<Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff, T> Mul<T> for &'a SeriesIn<Var, C>
where
    SeriesSliceIn<'a, Var, C>: Mul<T, Output = SeriesIn<Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> DivAssign<&'a SeriesIn<Var, C>>
    for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c SeriesIn<Var, C>: MulInverse<Output = SeriesIn<Var, C>>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// s /= &s.clone();
    /// let res = SeriesIn::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: &'a SeriesIn<Var, C>) {
        self.div_assign(other.as_slice(..));
    }
}

impl<Var: Clone, C: Coeff + SubAssign> DivAssign for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign + MulInverse<Output = SeriesIn<Var, C>>,
    for<'a> &'a C: Div<Output = C> + Mul<Output = C>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::SeriesIn;
    /// let mut s = SeriesIn::new("x", -3, vec!(1.,0.,-3.));
    /// s /= s.clone();
    /// let res = SeriesIn::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: SeriesIn<Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> DivAssign<SeriesSliceIn<'a, Var, C>>
    for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c SeriesIn<Var, C>: MulInverse<Output = SeriesIn<Var, C>>,
{
    fn div_assign(&mut self, other: SeriesSliceIn<'a, Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, Var, C: Coeff, T> Div<T> for &'a SeriesIn<Var, C>
where
    SeriesIn<Var, C>: Clone + DivAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res /= other;
        res
    }
}

impl<Var, C: Coeff, T> Div<T> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: DivAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn div(mut self, other: T) -> Self::Output {
        self /= other;
        self
    }
}

impl<Var, C: Coeff> Exp for SeriesIn<Var, C>
where
    for<'a> &'a C: Mul<Output = C>,
    for<'a> C: MulAssign<&'a C>,
    C: Clone
        + Div<Output = C>
        + Mul<Output = C>
        + AddAssign
        + Exp<Output = C>
        + From<i32>,
{
    type Output = Self;

    /// Computes the exponential of a series
    ///
    /// # Panics
    ///
    /// Panics if the series contains negative powers of the expansion
    /// variable
    fn exp(self) -> Self::Output {
        let coeff = self.exp_coeff();
        SeriesIn::new(self.var, 0, coeff)
    }
}

impl<'a, Var, C: Coeff> Exp for &'a SeriesIn<Var, C>
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
    type Output = SeriesIn<Var, C>;

    /// Computes the exponential of a series
    ///
    /// # Panics
    ///
    /// Panics if the series contains negative powers of the expansion
    /// variable
    fn exp(self) -> Self::Output {
        self.as_slice(..).exp()
    }
}

impl<Var, C: Coeff + From<i32>> Ln for SeriesIn<Var, C>
where
    for<'a> C: DivAssign<&'a C>,
    for<'a> &'a C: Mul<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
        + Ln<Output = C>
        + From<Var>,
    Var: Clone,
{
    type Output = Self;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(mut self) -> Self {
        assert!(!self.series.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = self.series.coeffs[0].clone();
        self.series.coeffs[0] = C::one();
        for i in 1..self.series.coeffs.len() {
            self.series.coeffs[i] /= &c_k0;
        }
        let a = self.series.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = if k0 != 0 {
            let var = self.var.clone();
            c_k0.ln() + C::from(k0 as i32) * C::from(var).ln()
        } else {
            c_k0.ln()
        };
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone());
            for i in 1..n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                let tmp = num_factor * (&a[n - i] * &b[i]);
                b[n] -= tmp;
            }
        }
        SeriesIn::new(self.var, 0, b)
    }
}

impl<'a, Var, C: Coeff> Ln for &'a SeriesIn<Var, C>
where
    for<'b> C: Div<&'b C, Output = C>,
    for<'b> &'b C: Mul<Output = C> + Ln<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
        + From<Var>
        + From<i32>,
    Var: Clone,
{
    type Output = SeriesIn<Var, C>;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(self) -> Self::Output {
        self.as_slice(..).ln()
    }
}

impl<Var, C: Coeff, T> Pow<T> for SeriesIn<Var, C>
where
    SeriesIn<Var, C>: Ln,
    <SeriesIn<Var, C> as Ln>::Output: Mul<T>,
    <<SeriesIn<Var, C> as Ln>::Output as std::ops::Mul<T>>::Output: Exp,
{
    type Output = <<<SeriesIn<Var, C> as Ln>::Output as std::ops::Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

impl<'a, Var, C: Coeff, T> Pow<T> for &'a SeriesIn<Var, C>
where
    for<'b> SeriesSliceIn<'b, Var, C>: Ln<Output = SeriesIn<Var, C>>,
    SeriesIn<Var, C>: Mul<T>,
    <SeriesIn<Var, C> as Mul<T>>::Output: Exp,
{
    type Output = <<SeriesIn<Var, C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        self.as_slice(..).pow(exponent)
    }
}

impl<Var, C: Coeff + Clone> AddAssign<C> for SeriesIn<Var, C>
where
    C: AddAssign<C>,
{
    fn add_assign(&mut self, rhs: C) {
        self.series.add_assign(rhs)
    }
}

impl<Var, C: Coeff + Clone> SubAssign<C> for SeriesIn<Var, C>
where
    C: Neg<Output = C> + SubAssign<C>,
{
    fn sub_assign(&mut self, rhs: C) {
        self.series.sub_assign(rhs)
    }
}

impl<'a, Var, C: Coeff> MulAssign<&'a C> for SeriesIn<Var, C>
where
    C: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a C) {
        self.series.mul_assign(rhs)
    }
}

impl<Var, C: Coeff> MulAssign<C> for SeriesIn<Var, C>
where
    for<'a> Series<C>: MulAssign<C>,
{
    fn mul_assign(&mut self, rhs: C) {
        self.series.mul_assign(rhs)
    }
}

impl<'a, Var, C: Coeff> DivAssign<&'a C> for SeriesIn<Var, C>
where
    C: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a C) {
        self.series.div_assign(rhs)
    }
}

impl<Var, C: Coeff> DivAssign<C> for SeriesIn<Var, C>
where
    for<'a> Series<C>: DivAssign<C>,
{
    fn div_assign(&mut self, rhs: C) {
        self.series.div_assign(rhs)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> From<SeriesSliceIn<'a, Var, C>>
    for SeriesIn<Var, C>
{
    fn from(s: SeriesSliceIn<'a, Var, C>) -> Self {
        SeriesIn::new(s.var.clone(), s.series.min_pow, s.series.coeffs.to_vec())
    }
}

/// Data parts of a series
///
/// # Example
///
/// ```rust
/// // destructure a series
/// let s = series::SeriesIn::new("x", -1, vec![1,2,3]);
/// let series::SeriesInParts{var, min_pow, coeffs} = s.into();
/// assert_eq!(var, "x");
/// assert_eq!(min_pow, -1);
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct SeriesInParts<Var, C> {
    pub var: Var,
    pub min_pow: isize,
    pub coeffs: Vec<C>,
}

impl<Var, C: Coeff> From<SeriesIn<Var, C>> for SeriesInParts<Var, C> {
    fn from(s: SeriesIn<Var, C>) -> Self {
        let SeriesIn {
            series,
            var,
        } = s;
        SeriesInParts {
            var,
            min_pow: series.min_pow,
            coeffs: series.coeffs,
        }
    }
}

impl<Var, C: Coeff> From<SeriesInParts<Var, C>> for SeriesIn<Var, C> {
    fn from(parts: SeriesInParts<Var, C>) -> Self {
        SeriesIn::new(parts.var, parts.min_pow, parts.coeffs)
    }
}

impl<Var, C: Coeff + From<i32>> LnVarFree for SeriesIn<Var, C>
where
    for<'a> C: DivAssign<&'a C>,
    for<'a> &'a C: Mul<Output = C>,
    C: Clone + SubAssign + Ln<Output = C> + Div<Output = C> + Mul<Output = C>,
{
    type Output = Self;

    fn ln_var_free(self) -> Self {
        let Self {
            series,
            var,
        } = self;
        let series = series.ln_var_free();
        Self {
            series,
            var
        }
    }
}

impl<Var, C: Coeff + From<i32>> SeriesIn<Var, C> {
    /// Calculate series to some integer power
    ///
    /// In contrast to the more general `pow` method this does _not_ require
    /// the variable type to be convertible to the coefficient type and gives
    /// an overall nicer output
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::SeriesIn::new("x", -1, vec![1.,3.,7.]);
    /// let s_to_minus_5 = series::SeriesIn::new("x", 5, vec![1.,-15.,100.]);
    /// assert_eq!(s.powi(-5), s_to_minus_5);
    /// ```
    pub fn powi(self, exp: i32) -> Self
    where
        for<'a> C: DivAssign<&'a C>,
        for<'a> &'a C: Mul<Output = C>,
        C: Clone
            + SubAssign
            + Ln<Output = C>
            + Div<Output = C>
            + Mul<Output = C>,
        Series<C>: Mul<C, Output = Series<C>>
            + Exp<Output = Series<C>>
            + MulInverse<Output = Series<C>>,
    {
        let Self {
            series,
            var,
        } = self;
        let series = series.powi(exp);
        Self {
            series,
            var
        }
    }
}

impl<Var, C: Coeff + From<i32>> ExpCoeff for SeriesIn<Var, C>
where
    for<'a> &'a C: Mul<Output = C>,
    for<'a> C: MulAssign<&'a C>,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
{
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        self.as_slice(..).exp_coeff()
    }
}

impl<'a, Var, C: Coeff + From<i32>> ExpCoeff for SeriesSliceIn<'a, Var, C>
where
    for<'c> &'c C: Mul<Output = C>,
    for<'c> C: MulAssign<&'c C>,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
{
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        self.series.exp_coeff()
    }
}
