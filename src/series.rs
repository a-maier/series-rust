use crate::ops::{Exp, Ln, Pow};
use crate::util::trim_start;
use crate::{traits::*, SeriesInParts};
use crate::{Coeff, IntoIter, Iter};
use crate::{SeriesIn, SeriesSlice};

use std::cmp::min;
use std::convert::From;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};

/// Laurent series in a single anonymous variable up to some power
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct Series<C: Coeff> {
    pub(crate) min_pow: isize,
    pub(crate) coeffs: Vec<C>,
}

impl<C: Coeff> Series<C> {
    /// Create a new series
    ///
    /// # Example
    ///
    /// This creates a series starting with the power -1 with
    /// coefficients 1, 2, 3. In other words, the series x^-1 + 2 +
    /// 3*x + O(x^2).
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// ```
    pub fn new(min_pow: isize, coeffs: Vec<C>) -> Self {
        let mut res = Series {
            min_pow,
            coeffs,
        };
        res.trim();
        res
    }

    /// Create a new series with a given cutoff power
    ///
    /// # Example
    ///
    /// This creates a series starting at power -1 with coefficients
    /// 1, 2, 3 and vanishing coefficients up to power 5. In other
    /// words, the series x^-1 + 2 + 3*x + O(x^5).
    /// ```rust
    /// let s = series::Series::with_cutoff(-1..5, vec!(1,2,3));
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the cutoff power is lower than the starting power
    ///
    pub fn with_cutoff(powers: Range<isize>, mut coeffs: Vec<C>) -> Self {
        let min_pow = powers.start;
        let cutoff_pow = powers.end;
        assert!(cutoff_pow >= min_pow);
        let len = (cutoff_pow - min_pow) as usize;
        // can't use resize here, because C is not Clone
        if len < coeffs.len() {
            coeffs.truncate(len)
        } else {
            let num_missing = len - coeffs.len();
            coeffs.reserve(num_missing);
            for _ in 0..num_missing {
                coeffs.push(C::zero());
            }
        }
        Series::new(min_pow, coeffs)
    }

    /// Turn into a series with a named expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3)).in_var("x");
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn in_var<Var>(self, var: Var) -> SeriesIn<Var, C> {
        let Self {
            min_pow,
            coeffs,
        } = self;
        SeriesIn::new(var, min_pow, coeffs)
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// assert_eq!(s.min_pow(), -1);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    /// Get the power of the expansion variable where the series is
    /// truncated.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3));
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
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// let mut iter = s.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        self.as_slice(..).iter()
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
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// assert_eq!(s.try_coeff(-5), None);
    /// assert_eq!(s.try_coeff(-2), None);
    /// assert_eq!(s.try_coeff(-1), Some(&1));
    /// assert_eq!(s.try_coeff(0), Some(&2));
    /// assert_eq!(s.try_coeff(1), Some(&3));
    /// assert_eq!(s.try_coeff(2), None);
    /// assert_eq!(s.try_coeff(5), None);
    /// ```
    pub fn try_coeff(&self, pow: isize) -> Option<&C> {
        self.as_slice(..).try_coeff(pow)
    }
}

impl<C: 'static + Coeff + Send + Sync> Series<C> {
    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// assert_eq!(s.coeff(-5), Some(&0));
    /// assert_eq!(s.coeff(-2), Some(&0));
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&2));
    /// assert_eq!(s.coeff(1), Some(&3));
    /// assert_eq!(s.coeff(2), None);
    /// assert_eq!(s.coeff(5), None);
    /// ```
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        self.as_slice(..).coeff(pow)
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, Range<isize>> for Series<C> {
    type Output = SeriesSlice<'a, C>;

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
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
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
        SeriesSlice::new(r.start, &self.coeffs[start..end])
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeInclusive<isize>> for Series<C> {
    type Output = SeriesSlice<'a, C>;

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
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
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
        SeriesSlice::new(start, &self.coeffs[ustart..=end])
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>> for Series<C> {
    type Output = SeriesSlice<'a, C>;

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
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
    /// let t = s.as_slice(..=1);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        SeriesSlice::new(self.min_pow, &self.coeffs[..=end])
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>> for Series<C> {
    type Output = SeriesSlice<'a, C>;

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
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
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
        SeriesSlice::new(r.start, &self.coeffs[start..])
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>> for Series<C> {
    type Output = SeriesSlice<'a, C>;

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
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
    /// let t = s.as_slice(..2);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        SeriesSlice::new(self.min_pow, &self.coeffs[..end])
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeFull> for Series<C> {
    type Output = SeriesSlice<'a, C>;

    /// A slice containing the complete series.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    /// let s = series::Series::new(-1, vec!(1,2,3,4));
    /// let t = s.as_slice(..);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert_eq!(t[-1], s[-1]);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        SeriesSlice::new(self.min_pow, &self.coeffs[r])
    }
}

impl<C: Coeff> Index<isize> for Series<C> {
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
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// assert_eq!(s[-1], 1);
    /// assert_eq!(s[0], 2);
    /// assert_eq!(s[1], 3);
    /// assert!(std::panic::catch_unwind(|| s[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| s[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow) as usize]
    }
}

impl<C: Coeff> std::iter::IntoIterator for Series<C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the series coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec!(1,2,3));
    /// let mut iter = s.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        (self.min_pow..).zip(self.coeffs)
    }
}

impl<'a, C: Coeff + SubAssign> MulInverse for &'a Series<C>
where
    C: Coeff + SubAssign,
    for<'c> &'c C: Div<Output = C> + Mul<Output = C>,
{
    type Output = Series<C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new(-1, vec!(1.,2.,3.));
    /// let s_inv = (&s).mul_inverse();
    /// let one = series::Series::new(0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        self.as_slice(..).mul_inverse()
    }
}

impl<C: Coeff + SubAssign> MulInverse for Series<C>
where
    for<'a> &'a Series<C>: MulInverse<Output = Series<C>>,
{
    type Output = Series<C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new(-1, vec!(1.,2.,3.));
    /// let s_inv = s.clone().mul_inverse();
    /// let one = series::Series::new(0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self).mul_inverse()
    }
}

impl<C: Coeff> Series<C> {
    fn trim(&mut self) {
        self.min_pow += trim_start(&mut self.coeffs, &C::zero()) as isize;
    }
}

impl<C: Coeff + Neg<Output = C>> Neg for Series<C> {
    type Output = Series<C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new(-3, vec!(-1.,0.,3.));
    /// assert_eq!(-s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        Series::new(self.min_pow, neg_coeff)
    }
}

impl<'a, C: Coeff> Neg for &'a Series<C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Series<C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new(-3, vec!(-1.,0.,3.));
    /// assert_eq!(-&s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, C: Coeff + Clone> AddAssign<&'a Series<C>> for Series<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// let t = Series::new(-1, vec!(3., 4., 5.));
    /// let res = Series::new(-3, vec!(1.,0.,0.));
    /// s += &t;
    /// assert_eq!(res, s);
    /// ```
    fn add_assign(&mut self, other: &'a Series<C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, C: Coeff + Clone> AddAssign<SeriesSlice<'a, C>> for Series<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: SeriesSlice<'a, C>) {
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

impl<C: Coeff> AddAssign<Series<C>> for Series<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// let t = Series::new(-1, vec!(3., 4., 5.));
    /// let res = Series::new(-3, vec!(1.,0.,0.));
    /// s += t;
    /// assert_eq!(res, s);
    /// ```
    fn add_assign(&mut self, mut other: Series<C>) {
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

impl<'a, C: Coeff + Clone, Rhs> Add<Rhs> for &'a Series<C>
where
    Series<C>: AddAssign<Rhs>,
{
    type Output = Series<C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<C: Coeff + Clone, Rhs> Add<Rhs> for Series<C>
where
    Series<C>: AddAssign<Rhs>,
{
    type Output = Series<C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, C: Coeff> SubAssign<&'a Series<C>> for Series<C>
where
    for<'c> &'c Series<C>: Neg<Output = Series<C>>,
    Series<C>: AddAssign<Series<C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// let res = Series::new(0, vec!());
    /// s -= &s.clone();
    /// assert_eq!(res, s);
    /// ```
    fn sub_assign(&mut self, other: &'a Series<C>) {
        *self += -other;
    }
}

impl<'a, C: Coeff> SubAssign<SeriesSlice<'a, C>> for Series<C>
where
    for<'c> SeriesSlice<'c, C>: Neg<Output = Series<C>>,
    Series<C>: AddAssign<Series<C>>,
{
    fn sub_assign(&mut self, other: SeriesSlice<'a, C>) {
        *self += -other;
    }
}

impl<C: Coeff> SubAssign<Series<C>> for Series<C>
where
    Series<C>: AddAssign + Neg<Output = Series<C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// let res = Series::new(0, vec!());
    /// s -= s.clone();
    /// assert_eq!(res, s);
    /// ```
    fn sub_assign(&mut self, other: Series<C>) {
        *self += -other;
    }
}

// TODO: somehow make addition symmetric?
impl<'a, C: Coeff, T> Sub<T> for &'a Series<C>
where
    Series<C>: Clone + SubAssign<T>,
{
    type Output = Series<C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res -= other;
        res
    }
}

impl<C: Coeff, T> Sub<T> for Series<C>
where
    Series<C>: SubAssign<T>,
{
    type Output = Series<C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, C: Coeff + Clone + AddAssign> MulAssign<&'a Series<C>> for Series<C>
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
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new(-6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    fn mul_assign(&mut self, other: &'a Series<C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, C> MulAssign<SeriesSlice<'a, C>> for Series<C>
where
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C> + Coeff + Clone + AddAssign,
{
    fn mul_assign(&mut self, other: SeriesSlice<'a, C>) {
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
        if let Some(c0) = self.coeffs.first_mut() {
            *c0 *= &other.coeffs[0]
        }
    }
}

impl<C: Coeff> MulAssign for Series<C>
where
    for<'a> Series<C>: MulAssign<&'a Series<C>>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new(-6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    fn mul_assign(&mut self, other: Series<C>) {
        *self *= &other
    }
}

// TODO: somehow make multiplication symmetric?
impl<C: Coeff> Mul for Series<C>
where
    Series<C>: MulAssign,
{
    type Output = Series<C>;

    fn mul(mut self, other: Series<C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a Series<C>> for Series<C>
where
    Series<C>: MulAssign<SeriesSlice<'a, C>>,
{
    type Output = Series<C>;

    fn mul(self, other: &'a Series<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, C: Coeff> Mul<SeriesSlice<'a, C>> for Series<C>
where
    Series<C>: MulAssign<SeriesSlice<'a, C>>,
{
    type Output = Series<C>;

    fn mul(mut self, other: SeriesSlice<'a, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff> Mul<C> for Series<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Series<C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a C> for Series<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Series<C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff, T> Mul<T> for &'a Series<C>
where
    SeriesSlice<'a, C>: Mul<T, Output = Series<C>>,
{
    type Output = Series<C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, C: Coeff + SubAssign> DivAssign<&'a Series<C>> for Series<C>
where
    Series<C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c Series<C>: MulInverse<Output = Series<C>>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// s /= &s.clone();
    /// let res = Series::new(0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    fn div_assign(&mut self, other: &'a Series<C>) {
        self.div_assign(other.as_slice(..));
    }
}

impl<C: Coeff + SubAssign> DivAssign for Series<C>
where
    Series<C>: MulAssign + MulInverse<Output = Series<C>>,
    for<'a> &'a C: Div<Output = C> + Mul<Output = C>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new(-3, vec!(1.,0.,-3.));
    /// s /= s.clone();
    /// let res = Series::new(0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: Series<C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, C: Coeff + SubAssign> DivAssign<SeriesSlice<'a, C>> for Series<C>
where
    Series<C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c Series<C>: MulInverse<Output = Series<C>>,
{
    fn div_assign(&mut self, other: SeriesSlice<'a, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, C: Coeff, T> Div<T> for &'a Series<C>
where
    Series<C>: Clone + DivAssign<T>,
{
    type Output = Series<C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res /= other;
        res
    }
}

impl<C: Coeff, T> Div<T> for Series<C>
where
    Series<C>: DivAssign<T>,
{
    type Output = Series<C>;

    fn div(mut self, other: T) -> Self::Output {
        self /= other;
        self
    }
}

impl<C: Coeff> Exp for Series<C>
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
        Series::new(0, coeff)
    }
}

impl<'a, C: Coeff> Exp for &'a Series<C>
where
    for<'b> &'b C: Mul<Output = C>,
    for<'b> C: MulAssign<&'b C>,
    C: Clone
        + Div<Output = C>
        + Mul<Output = C>
        + AddAssign
        + Exp<Output = C>
        + From<i32>,
{
    type Output = Series<C>;

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

impl<C: Coeff> Ln for Series<C>
where
    for<'a> C: DivAssign<&'a C>,
    for<'a> &'a C: Mul<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
        + Ln<Output = C>
        + From<i32>,
{
    type Output = Self;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has only vanishing coefficients or does
    /// not start with power 0. Adjoin a variable with `in_var` to
    /// compute the logarithm of a series with a non-vanishing leading
    /// power.
    fn ln(self) -> Self {
        assert_eq!(self.min_pow(), 0);
        assert!(!self.coeffs.is_empty());
        self.ln_var_free()
    }
}

impl<'a, C: Coeff> Ln for &'a Series<C>
where
    for<'b> C: Div<&'b C, Output = C>,
    for<'b> &'b C: Mul<Output = C> + Ln<Output = C>,
    C: Clone
        + SubAssign
        + Add<Output = C>
        + Mul<Output = C>
        + Div<Output = C>
        + From<i32>,
{
    type Output = Series<C>;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has only vanishing coefficients or does
    /// not start with power 0. Adjoin a variable with `in_var` to
    /// compute the logarithm of a series with a non-vanishing leading
    /// power.
    fn ln(self) -> Self::Output {
        self.as_slice(..).ln()
    }
}

impl<C: Coeff, T> Pow<T> for Series<C>
where
    Series<C>: Ln,
    <Series<C> as Ln>::Output: Mul<T>,
    <<Series<C> as Ln>::Output as std::ops::Mul<T>>::Output: Exp,
{
    type Output = <<<Series<C> as Ln>::Output as std::ops::Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

impl<'a, C: Coeff, T> Pow<T> for &'a Series<C>
where
    for<'b> SeriesSlice<'b, C>: Ln<Output = Series<C>>,
    Series<C>: Mul<T>,
    <Series<C> as Mul<T>>::Output: Exp,
{
    type Output = <<Series<C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        self.as_slice(..).pow(exponent)
    }
}

impl<C: Coeff + Clone> AddAssign<C> for Series<C>
where
    C: AddAssign<C>,
{
    fn add_assign(&mut self, rhs: C) {
        if self.cutoff_pow() <= 0 || rhs.is_zero() {
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
            new_coeffs.resize(self.min_pow() as usize, C::zero());
            new_coeffs.append(&mut self.coeffs);
            self.coeffs = new_coeffs;
            self.min_pow = 0;
        }
    }
}

impl<C: Coeff + Clone> SubAssign<C> for Series<C>
where
    C: Neg<Output = C> + SubAssign<C>,
{
    fn sub_assign(&mut self, rhs: C) {
        if self.cutoff_pow() <= 0 || rhs.is_zero() {
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
            new_coeffs.resize(self.min_pow() as usize, C::zero());
            new_coeffs.append(&mut self.coeffs);
            self.coeffs = new_coeffs;
            self.min_pow = 0;
        }
    }
}

impl<'a, C: Coeff> MulAssign<&'a C> for Series<C>
where
    C: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff *= rhs
        }
    }
}

impl<C: Coeff> MulAssign<C> for Series<C>
where
    for<'a> Series<C>: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: C) {
        *self *= &rhs
    }
}

impl<'a, C: Coeff> DivAssign<&'a C> for Series<C>
where
    C: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff /= rhs
        }
    }
}

impl<C: Coeff> DivAssign<C> for Series<C>
where
    for<'a> Series<C>: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: C) {
        *self /= &rhs
    }
}

impl<'a, C: Coeff + Clone> From<SeriesSlice<'a, C>> for Series<C> {
    fn from(s: SeriesSlice<'a, C>) -> Self {
        Series::new(s.min_pow, s.coeffs.to_vec())
    }
}

impl<C: Coeff, Var> From<SeriesIn<Var, C>> for Series<C> {
    fn from(source: SeriesIn<Var, C>) -> Self {
        let SeriesInParts {
            var: _,
            min_pow,
            coeffs,
        } = source.into();
        Self {
            min_pow,
            coeffs,
        }
    }
}

/// Data parts of a series
///
/// # Example
///
/// ```rust
/// // destructure a series
/// let s = series::Series::new(-1, vec![1,2,3]);
/// let series::SeriesParts{min_pow, coeffs} = s.into();
/// assert_eq!(min_pow, -1);
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct SeriesParts<C> {
    pub min_pow: isize,
    pub coeffs: Vec<C>,
}

impl<C: Coeff> From<Series<C>> for SeriesParts<C> {
    fn from(s: Series<C>) -> Self {
        SeriesParts {
            min_pow: s.min_pow,
            coeffs: s.coeffs,
        }
    }
}

impl<C: Coeff> From<SeriesParts<C>> for Series<C> {
    fn from(parts: SeriesParts<C>) -> Self {
        Series::new(parts.min_pow, parts.coeffs)
    }
}

impl<C: Coeff> LnVarFree for Series<C>
where
    for<'a> C: DivAssign<&'a C>,
    for<'a> &'a C: Mul<Output = C>,
    C: Clone
        + SubAssign
        + Ln<Output = C>
        + Div<Output = C>
        + Mul<Output = C>
        + From<i32>,
{
    type Output = Self;

    fn ln_var_free(mut self) -> Self {
        debug_assert!(self.min_pow == 0);
        assert!(!self.coeffs.is_empty());
        let c_k0 = self.coeffs[0].clone();
        self.coeffs[0] = C::one();
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
        Series::new(0, b)
    }
}

impl<C: Coeff> Series<C> {
    /// Calculate series to some integer power
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new(-1, vec![1.,3.,7.]);
    /// let s_to_minus_5 = series::Series::new(5, vec![1.,-15.,100.]);
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
            + Mul<Output = C>
            + From<i32>,
        Series<C>: Mul<C, Output = Self>
            + Exp<Output = Self>
            + MulInverse<Output = Self>,
    {
        if self.coeffs.is_empty() {
            self.min_pow *= exp as isize;
            return self;
        }
        let new_min_pow = self.min_pow * (exp as isize);
        self.min_pow = 0;
        let pow = (self.ln_var_free() * C::from(exp)).exp();

        Series::new(new_min_pow, pow.coeffs)
    }
}

impl<C: Coeff> Series<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn truncate_cutoff_pow(&mut self, other: SeriesSlice<'_, C>) {
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

    fn add_overlap(&mut self, other: SeriesSlice<'_, C>) {
        let offset = self.min_pow();
        for (i, c) in self.coeffs.iter_mut().enumerate() {
            let power = offset + i as isize;
            if let Some(coeff) = other.try_coeff(power) {
                *c += coeff
            }
        }
    }

    fn num_leading(&mut self, other: SeriesSlice<'_, C>) -> usize {
        min(
            (self.min_pow() - other.min_pow()) as usize,
            other.coeffs.len(),
        )
    }
}

impl<C: Coeff> ExpCoeff for Series<C>
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
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        self.as_slice(..).exp_coeff()
    }
}

impl<'a, C: Coeff> ExpCoeff for SeriesSlice<'a, C>
where
    for<'c> &'c C: Mul<Output = C>,
    for<'c> C: MulAssign<&'c C>,
    C: Clone
        + Div<Output = C>
        + Mul<Output = C>
        + AddAssign
        + Exp<Output = C>
        + From<i32>,
{
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        assert!(self.min_pow() >= 0);
        let mut b = Vec::with_capacity(min(self.coeffs.len(), 1));
        b.push(C::one());
        debug_assert!(self.cutoff_pow() >= 0);
        for n in 1..self.cutoff_pow() as usize {
            let mut b_n = C::zero();
            for i in 1..=n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                if let Some(a_i) = self.try_coeff(i as isize) {
                    b_n += num_factor * (a_i * &b[n - i]);
                }
            }
            b.push(b_n);
        }
        if self.min_pow() == 0 {
            let exp_a_0 = self.try_coeff(0).cloned()
                .unwrap_or_else(|| C::zero()).exp();
            for b_n in &mut b {
                *b_n *= &exp_a_0;
            }
        }
        b
    }
}
