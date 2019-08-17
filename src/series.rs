use crate::ops::{Exp, Ln, Pow};
use crate::slice::*;
use crate::traits::*;
use crate::util::trim_start;
use crate::{Coeff, Iter, IterMut, IntoIter};

use std::cmp::min;
use std::convert::From;
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Range,
    RangeInclusive, RangeToInclusive, RangeFrom, RangeFull, RangeTo,
    Sub, SubAssign, Index, IndexMut
};

/// Laurent series in a single variable up to some power
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct Series<Var, C: Coeff> {
    pub(crate) var: Var,
    pub(crate) min_pow: isize,
    pub(crate) coeffs: Vec<C>,
    pub(crate) zero: C, // TODO: evil hack, learn how to do this better
}

impl<Var, C: Coeff> Series<Var, C> {
    /// Create a new series
    ///
    /// # Example
    ///
    /// This creates a series in the variable "x", starting at "x"^-1
    /// with coefficients 1, 2, 3. In other words, the series x^-1 + 2 +
    /// 3*x + O(x^2).
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// ```
    pub fn new(var: Var, min_pow: isize, coeffs: Vec<C>) -> Series<Var, C> {
        let mut res = Series {
            var,
            min_pow,
            coeffs,
            zero: C::from(0),
        };
        res.trim();
        res
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
    /// let s = series::Series::with_cutoff("x", -1, 5, vec!(1,2,3));
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
        mut coeffs: Vec<C>,
    ) -> Series<Var, C> {
        assert!(cutoff_pow >= min_pow);
        let len = (cutoff_pow - min_pow) as usize;
        // can't use resize here, because C is not Clone
        if len < coeffs.len() {
            coeffs.truncate(len)
        } else {
            let num_missing = len - coeffs.len();
            coeffs.reserve(num_missing);
            for _ in 0..num_missing {
                coeffs.push(C::from(0));
            }
        }
        Series::new(var, min_pow, coeffs)
    }

    /// Get the expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        &self.var
    }

    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.coeff(-5), Some(&0));
    /// assert_eq!(s.coeff(-2), Some(&0));
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&2));
    /// assert_eq!(s.coeff(1), Some(&3));
    /// assert_eq!(s.coeff(2), None);
    /// assert_eq!(s.coeff(5), None);
    /// ```
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        if pow < self.min_pow() {
            return Some(&self.zero); // TODO this is a bad hack
        }
        if pow >= self.cutoff_pow() {
            return None;
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(&self.coeffs[idx])
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.min_pow(), -1);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    /// Get the power of the expansion variable where the series is
    /// truncated, that is \\( N \\) in a series of the form \\(
    /// \sum_{n=n_0}^{N-1} a_n x^n + O(x^N), \\)
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.cutoff_pow(), 2);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.as_slice(..).cutoff_pow()
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let mut iter = s.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        // TODO: code duplication with SeriesSlice
        (self.min_pow..).zip(self.coeffs.iter())
    }

    /// An iterator that allows modifying each series coefficient.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut s = series::Series::new("x", -1, vec!(1,2,3));
    /// for (_pow, coeff) in s.iter_mut() {
    ///     *coeff += 1
    /// }
    /// let inc = series::Series::new("x", -1, vec!(2,3,4));
    /// assert_eq!(s, inc);
    /// ```
    pub fn iter_mut(&mut self) -> IterMut<C> {
        (self.min_pow..).zip(self.coeffs.iter_mut())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, Range<isize>> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

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
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..2);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        let start = (r.start - self.min_pow()) as usize;
        let end = (r.end - self.min_pow()) as usize;
        SeriesSlice::new(
            &self.var,
            r.start,
            &self.coeffs[start..end],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeInclusive<isize>> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

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
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..=1);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        let (start, end) = r.into_inner();
        let ustart = (start - self.min_pow()) as usize;
        let end = (end - self.min_pow()) as usize;
        SeriesSlice::new(
            &self.var,
            start,
            &self.coeffs[ustart..=end],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeToInclusive<isize>> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..=1);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        SeriesSlice::new(
            &self.var,
            self.min_pow,
            &self.coeffs[..=end],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeFrom<isize>> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the lower bound is smaller than the leading power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(0..);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert!(std::panic::catch_unwind(|| t[-1]).is_err());
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        let start = (r.start - self.min_pow()) as usize;
        SeriesSlice::new(
            &self.var,
            r.start,
            &self.coeffs[start..],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeTo<isize>> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..2);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        SeriesSlice::new(
            &self.var,
            self.min_pow,
            &self.coeffs[..end],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeFull> for Series<Var, C> {
    type Output = SeriesSlice<'a, Var, C>;

    /// A slice containing the complete series.
    ///
    /// # Example
    ///
    /// ```rust
    /// use crate::series::AsSlice;
    /// let s = series::Series::new("x", -1, vec!(1,2,3,4));
    /// let t = s.as_slice(..);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert_eq!(t[-1], s[-1]);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        SeriesSlice::new(
            &self.var,
            self.min_pow,
            &self.coeffs[r],
            &self.zero,
        )
    }
}

impl<Var, C: Coeff> Index<isize> for Series<Var, C> {
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
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s[-1], 1);
    /// assert_eq!(s[0], 2);
    /// assert_eq!(s[1], 3);
    /// assert!(std::panic::catch_unwind(|| s[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| s[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index-self.min_pow) as usize]
    }
}

impl<Var, C: Coeff> IndexMut<isize> for Series<Var, C> {
    /// Access the (mutable) series coefficient of the expansion
    /// variable to the given power.
    ///
    /// # Panics
    ///
    /// Panics if the index is smaller than the leading power or
    /// at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut s = series::Series::new("x", -1, vec!(1,2,3));
    /// s[-1] = 0;
    /// assert_eq!(s[-1], 0);
    /// assert_eq!(s[0], 2);
    /// assert_eq!(s[1], 3);
    /// assert!(std::panic::catch_unwind(|| s[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| s[2]).is_err());
    /// ```
    fn index_mut(&mut self, index: isize) -> &mut Self::Output {
        &mut self.coeffs[(index-self.min_pow) as usize]
    }
}

impl<Var, C: Coeff> std::iter::IntoIterator for Series<Var, C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the series coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// let mut iter = s.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        (self.min_pow..).zip(self.coeffs.into_iter())
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> MulInverse for &'a Series<Var, C>
where
    Var: Clone,
    C: Coeff + SubAssign,
    for<'c> &'c C: Div<Output = C> + Mul<Output = C>,
{
    type Output = Series<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = (&s).mul_inverse();
    /// let one = series::Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        self.as_slice(..).mul_inverse()
    }
}

impl<Var: Clone, C: Coeff + SubAssign> MulInverse for Series<Var, C>
where
    for<'a> &'a Series<Var, C>: MulInverse<Output = Series<Var, C>>,
{
    type Output = Series<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = s.clone().mul_inverse();
    /// let one = series::Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self).mul_inverse()
    }
}

impl<Var, C: Coeff> Series<Var, C> {
    fn trim(&mut self) {
        self.min_pow += trim_start(&mut self.coeffs, &C::from(0)) as isize;
    }
}

impl<Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for Series<Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.as_slice(..).fmt(f)
    }
}

impl<Var, C: Coeff + Neg<Output = C>> Neg for Series<Var, C> {
    type Output = Series<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        Series::new(self.var, self.min_pow, neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a Series<Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Series<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-&s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<&'a Series<Var, C>> for Series<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = Series::new("x", -1, vec!(3., 4., 5.));
    /// let res = Series::new("x", -3, vec!(1.,0.,0.));
    /// s += &t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(&mut self, other: &'a Series<Var, C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<SeriesSlice<'a, Var, C>> for Series<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: SeriesSlice<'a, Var, C>) {
        assert_eq!(self.var, *other.var);
        self.truncate_cutoff_pow(other);
        self.add_overlap(other);
        if other.min_pow() < self.min_pow() {
            let num_leading = self.num_leading(other);
            let leading_coeff = other.coeffs[0..num_leading].iter().cloned();
            self.coeffs.splice(0..0, leading_coeff);
            self.min_pow = other.min_pow;
        }
        debug_assert!(other.cutoff_pow() >= self.cutoff_pow());
        self.trim();
    }
}

impl<Var, C: Coeff> AddAssign<Series<Var, C>> for Series<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
    Var: PartialEq + fmt::Debug,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = Series::new("x", -1, vec!(3., 4., 5.));
    /// let res = Series::new("x", -3, vec!(1.,0.,0.));
    /// s += t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(&mut self, mut other: Series<Var, C>) {
        assert_eq!(self.var, other.var);
        self.truncate_cutoff_pow(other.as_slice(..));
        self.add_overlap(other.as_slice(..));
        if other.min_pow() < self.min_pow() {
            let num_leading = self.num_leading(other.as_slice(..));
            let leading_coeff = other.coeffs.drain(0..num_leading);
            self.coeffs.splice(0..0, leading_coeff);
            self.min_pow = other.min_pow;
        }
        self.trim();
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for &'a Series<Var, C>
where
    Series<Var, C>: AddAssign<Rhs>,
{
    type Output = Series<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for Series<Var, C>
where
    Series<Var, C>: AddAssign<Rhs>,
{
    type Output = Series<Var, C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, Var, C: Coeff> SubAssign<&'a Series<Var, C>> for Series<Var, C>
where
    for<'c> &'c Series<Var, C>: Neg<Output = Series<Var, C>>,
    Series<Var, C>: AddAssign<Series<Var, C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!());
    /// s -= &s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(&mut self, other: &'a Series<Var, C>) {
        *self += -other;
    }
}

impl<'a, Var, C: Coeff> SubAssign<SeriesSlice<'a, Var, C>> for Series<Var, C>
where
    for<'c> SeriesSlice<'c, Var, C>: Neg<Output = Series<Var, C>>,
    Series<Var, C>: AddAssign<Series<Var, C>>,
{
    fn sub_assign(&mut self, other: SeriesSlice<'a, Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff> SubAssign<Series<Var, C>> for Series<Var, C>
where
    Series<Var, C>: AddAssign + Neg<Output = Series<Var, C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!());
    /// s -= s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(&mut self, other: Series<Var, C>) {
        *self += -other;
    }
}

// TODO: somehow make addition symmetric?
impl<'a, Var, C: Coeff, T> Sub<T> for &'a Series<Var, C>
where
    Series<Var, C>: Clone + SubAssign<T>,
{
    type Output = Series<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res -= other;
        res
    }
}

impl<Var, C: Coeff, T> Sub<T> for Series<Var, C>
where
    Series<Var, C>: SubAssign<T>,
{
    type Output = Series<Var, C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone + AddAssign>
    MulAssign<&'a Series<Var, C>> for Series<Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: &'a Series<Var, C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, Var, C: > MulAssign<SeriesSlice<'a, Var, C>> for Series<Var, C>
where
    Var: PartialEq + fmt::Debug,
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C> + Coeff + Clone + AddAssign,
{
    fn mul_assign(&mut self, other: SeriesSlice<'a, Var, C>) {
        assert_eq!(self.var, *other.var);
        self.min_pow += other.min_pow();
        let num_coeffs = min(self.coeffs.len(), other.coeffs.len());
        self.coeffs.truncate(num_coeffs);
        // compute Cauchy product
        for k in (1..self.coeffs.len()).rev() {
            let (c_k, c) = self.coeffs[..=k].split_last_mut().unwrap();
            *c_k *= &other.coeffs[0];
            for i in 1..=k {
                *c_k += &c[k - i] * &other.coeffs[i]
            }
        }
        self.coeffs.first_mut().map(|c| *c *= &other.coeffs[0]);
    }
}

impl<Var, C: Coeff> MulAssign for Series<Var, C>
where
    for<'a> Series<Var, C>: MulAssign<&'a Series<Var, C>>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: Series<Var, C>) {
        *self *= &other
    }
}

// TODO: somehow make multiplication symmetric?
impl<'a, Var, C: Coeff, T> Mul<T> for &'a Series<Var, C>
where
    Series<Var, C>: Clone + MulAssign<T>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res *= other;
        res
    }
}

impl<Var, C: Coeff, T> Mul<T> for Series<Var, C>
where
    Series<Var, C>: MulAssign<T>,
{
    type Output = Series<Var, C>;

    fn mul(mut self, other: T) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> DivAssign<&'a Series<Var, C>>
    for Series<Var, C>
where
    Series<Var, C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c Series<Var, C>: MulInverse<Output = Series<Var, C>>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s /= &s.clone();
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: &'a Series<Var, C>) {
        self.div_assign(other.as_slice(..));
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> DivAssign for Series<Var, C>
where
    Series<Var, C>: MulAssign + MulInverse<Output = Series<Var, C>>,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s /= s.clone();
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: Series<Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign> DivAssign<SeriesSlice<'a, Var, C>>
    for Series<Var, C>
where
    Series<Var, C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c Series<Var, C>: MulInverse<Output = Series<Var, C>>,
{
    fn div_assign(&mut self, other: SeriesSlice<'a, Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, Var, C: Coeff, T> Div<T> for &'a Series<Var, C>
where
    Series<Var, C>: Clone + DivAssign<T>,
{
    type Output = Series<Var, C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res /= other;
        res
    }
}

impl<Var, C: Coeff, T> Div<T> for Series<Var, C>
where
    Series<Var, C>: DivAssign<T>,
{
    type Output = Series<Var, C>;

    fn div(mut self, other: T) -> Self::Output {
        self /= other;
        self
    }
}

impl<Var, C: Coeff> Exp for Series<Var, C>
where
    for<'a> &'a C: Mul<Output = C>,
    for<'a> C: MulAssign<&'a C>,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
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
        Series::new(self.var, 0, coeff)
    }
}

impl<'a, Var, C: Coeff> Exp for &'a Series<Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
    for<'b> C: MulAssign<&'b C>,
    Var: Clone,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
{
    type Output = Series<Var, C>;

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

impl<Var, C: Coeff> Ln for Series<Var, C>
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
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = self.coeffs[0].clone();
        self.coeffs[0] = C::from(1);
        for i in 1..self.coeffs.len() {
            self.coeffs[i] /= &c_k0;
        }
        let a = self.coeffs;
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
        Series::new(self.var, 0, b)
    }
}

impl<'a, Var, C: Coeff> Ln for &'a Series<Var, C>
where
    for<'b> C: Div<&'b C, Output = C>,
    for<'b> &'b C: Mul<Output = C> + Ln<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
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
        self.as_slice(..).ln()
    }
}

impl<Var, C: Coeff, T> Pow<T> for Series<Var, C>
where
    Series<Var, C>: Ln,
    <Series<Var, C> as Ln>::Output: Mul<T>,
    <<Series<Var, C> as Ln>::Output as std::ops::Mul<T>>::Output: Exp,
{
    type Output = <<<Series<Var, C> as Ln>::Output as std::ops::Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

impl<'a, Var, C: Coeff, T> Pow<T> for &'a Series<Var, C>
where
    for<'b> SeriesSlice<'b, Var, C>: Ln<Output = Series<Var, C>>,
    Series<Var, C>: Mul<T>,
    <Series<Var, C> as Mul<T>>::Output: Exp,
{
    type Output = <<Series<Var, C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        self.as_slice(..).pow(exponent)
    }
}

impl<Var, C: Coeff + Clone> AddAssign<C> for Series<Var, C>
where
    C: AddAssign<C>,
{
    fn add_assign(&mut self, rhs: C) {
        if self.cutoff_pow() <= 0 || rhs == self.zero {
            return;
        }
        if self.min_pow() <= 0 {
            let idx = (-self.min_pow()) as usize;
            self.coeffs[idx] += rhs;
            if self.min_pow() == 0 {
                self.trim()
            }
        } else {
            let mut new_coeffs = vec![rhs];
            new_coeffs.resize(self.min_pow() as usize, self.zero.clone());
            new_coeffs.extend(self.coeffs.drain(..));
            self.coeffs = new_coeffs;
            self.min_pow = 0;
        }
    }
}

impl<Var, C: Coeff + Clone> SubAssign<C> for Series<Var, C>
where
    C: Neg<Output = C> + SubAssign<C>,
{
    fn sub_assign(&mut self, rhs: C) {
        if self.cutoff_pow() <= 0 || rhs == self.zero {
            return;
        }
        if self.min_pow() <= 0 {
            let idx = (-self.min_pow()) as usize;
            self.coeffs[idx] -= rhs;
            if self.min_pow() == 0 {
                self.trim()
            }
        } else {
            let mut new_coeffs = vec![-rhs];
            new_coeffs.resize(self.min_pow() as usize, self.zero.clone());
            new_coeffs.extend(self.coeffs.drain(..));
            self.coeffs = new_coeffs;
            self.min_pow = 0;
        }
    }
}

impl<'a, Var, C: Coeff> MulAssign<&'a C> for Series<Var, C>
where
    C: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff *= rhs
        }
    }
}

impl<Var, C: Coeff> MulAssign<C> for Series<Var, C>
where
    for<'a> Series<Var, C>: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: C) {
        *self *= &rhs
    }
}

impl<'a, Var, C: Coeff> DivAssign<&'a C> for Series<Var, C>
where
    C: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff /= rhs
        }
    }
}

impl<Var, C: Coeff> DivAssign<C> for Series<Var, C>
where
    for<'a> Series<Var, C>: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: C) {
        *self /= &rhs
    }
}

/// Data parts of a series
///
/// # Example
///
/// ```rust
/// // destructure a series
/// let s = series::Series::new("x", -1, vec![1,2,3]);
/// let series::SeriesParts{var, min_pow, coeffs} = s.into();
/// assert_eq!(var, "x");
/// assert_eq!(min_pow, -1);
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct SeriesParts<Var, C> {
    pub var: Var,
    pub min_pow: isize,
    pub coeffs: Vec<C>,
}

impl<Var, C: Coeff> From<Series<Var, C>> for SeriesParts<Var, C> {
    fn from(s: Series<Var, C>) -> Self {
        SeriesParts {
            var: s.var,
            min_pow: s.min_pow,
            coeffs: s.coeffs,
        }
    }
}

impl<Var, C: Coeff> From<SeriesParts<Var, C>> for Series<Var, C> {
    fn from(parts: SeriesParts<Var, C>) -> Self {
        Series::new(parts.var, parts.min_pow, parts.coeffs)
    }
}

impl<Var, C: Coeff> LnVarFree for Series<Var, C>
where
    for<'a> C: DivAssign<&'a C>,
    for<'a> &'a C: Mul<Output = C>,
    C: Clone + SubAssign + Ln<Output = C> + Div<Output = C> + Mul<Output = C>,
{
    type Output = Self;

    fn ln_var_free(mut self) -> Self {
        debug_assert!(self.min_pow == 0);
        assert!(!self.coeffs.is_empty());
        let c_k0 = self.coeffs[0].clone();
        self.coeffs[0] = C::from(1);
        for i in 1..self.coeffs.len() {
            self.coeffs[i] /= &c_k0;
        }
        let a = self.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = c_k0.ln();
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone());
            for i in 1..n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                let tmp = num_factor * (&a[n - i] * &b[i]);
                b[n] -= tmp;
            }
        }
        Series::new(self.var, 0, b)
    }
}

impl<Var, C: Coeff> Series<Var, C> {
    /// Calculate series to some integer power
    ///
    /// In contrast to the more general `pow` method this does _not_ require
    /// the variable type to be convertible to the coefficient type and gives
    /// an overall nicer output
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec![1.,3.,7.]);
    /// let s_to_minus_5 = series::Series::new("x", 5, vec![1.,-15.,100.]);
    /// assert_eq!(s.powi(-5), s_to_minus_5);
    /// ```
    pub fn powi(mut self, exp: i32) -> Self
    where
        for<'a> C: DivAssign<&'a C>,
        for<'a> &'a C: Mul<Output = C>,
        C: Clone
            + SubAssign
            + Ln<Output = C>
            + Div<Output = C>
            + Mul<Output = C>,
        Series<Var, C>: Mul<C, Output = Self>
            + Exp<Output = Self>
            + MulInverse<Output = Self>,
    {
        let new_min_pow = self.min_pow * (exp as isize);
        self.min_pow = 0;
        let pow = (self.ln_var_free() * C::from(exp)).exp();

        Series::new(pow.var, new_min_pow, pow.coeffs)
    }
}

impl<Var, C: Coeff> AddAssignHelper<Var, C> for Series<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn truncate_cutoff_pow<'a>(&mut self, other: SeriesSlice<'a, Var, C>) {
        if other.cutoff_pow() < self.cutoff_pow() {
            let to_remove = min(
                (self.cutoff_pow() - other.cutoff_pow()) as usize,
                self.coeffs.len(),
            );
            let new_size = self.coeffs.len() - to_remove;
            self.coeffs.truncate(new_size);
            debug_assert!(
                self.coeffs.is_empty()
                    || other.cutoff_pow() == self.cutoff_pow()
            );
        }
    }

    fn add_overlap<'a>(&mut self, other: SeriesSlice<'a, Var, C>) {
        let offset = self.min_pow();
        for (i, c) in self.coeffs.iter_mut().enumerate() {
            let power = offset + i as isize;
            *c += other.coeff(power).unwrap();
        }
    }

    fn num_leading<'a>(&mut self, other: SeriesSlice<'a, Var, C>) -> usize {
        min(
            (self.min_pow() - other.min_pow()) as usize,
            other.coeffs.len(),
        )
    }
}

impl<Var, C: Coeff> ExpCoeff for Series<Var, C>
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

impl<'a, Var, C: Coeff> ExpCoeff for SeriesSlice<'a, Var, C>
where
    for<'c> &'c C: Mul<Output = C>,
    for<'c> C: MulAssign<&'c C>,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
{
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        assert!(self.min_pow() >= 0);
        let mut b = Vec::with_capacity(min(self.coeffs.len(), 1));
        b.push(C::from(1));
        debug_assert!(self.cutoff_pow() >= 0);
        for n in 1..self.cutoff_pow() as usize {
            let mut b_n = C::from(0);
            for i in 1..=n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                let a_i = self.coeff(i as isize).unwrap();
                b_n += num_factor * (a_i * &b[n - i]);
            }
            b.push(b_n);
        }
        if self.min_pow() == 0 {
            let exp_a_0 = self.coeff(0).unwrap().clone().exp();
            for b_n in &mut b {
                *b_n *= &exp_a_0;
            }
        }
        b
    }
}
