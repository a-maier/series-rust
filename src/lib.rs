/// Laurent series with support for usual arithmetic operations and some
/// common functions
#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
use std::cmp::min;
use std::convert::From;
use std::fmt;
use std::iter::Zip;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub,
    SubAssign, Index, IndexMut,
    Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive,
};

pub mod ops;
use self::ops::{Exp, Ln, Pow};

/// Minimum requirements on series coefficients
pub trait Coeff: From<i32> + PartialEq {}
impl<T: From<i32> + PartialEq> Coeff for T {}

/// Immutable `Series` iterator.
///
/// This `struct` is created by the `iter` method on `Series`
pub type Iter<'a, C> = Zip<RangeFrom<isize>, std::slice::Iter<'a, C>>;
/// Mutable `Series` iterator.
///
/// This `struct` is created by the `iter_mut` method on `Series`
pub type IterMut<'a, C> = Zip<RangeFrom<isize>, std::slice::IterMut<'a, C>>;
/// An iterator that moves out of a vector.
///
/// This `struct` is created by the `into_iter` method on `Series`
pub type IntoIter<C> = Zip<RangeFrom<isize>, std::vec::IntoIter<C>>;

/// Laurent series in a single variable up to some power
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct Series<Var, C: Coeff> {
    var: Var,
    min_pow: isize,
    coeffs: Vec<C>,
    zero: C, // TODO: evil hack, learn how to do this better
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

    /// Iterator over the series coefficients.
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
    /// for (pow, coeff) in s.iter_mut() {
    ///     *coeff += 1
    /// }
    /// let inc = series::Series::new("x", -1, vec!(2,3,4));
    /// assert_eq!(s, inc);
    /// ```
    pub fn iter_mut(&mut self) -> IterMut<C> {
        (self.min_pow..).zip(self.coeffs.iter_mut())
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

pub trait AsSlice<'a, Var, C: Coeff, T> {
    // TODO:
    // we would like to use an associated Output type
    // as return type, but that makes it almost impossible
    // to have the lifetime parameter in the return type
    fn as_slice(&'a self, index: T) -> SeriesSlice<'a, Var, C>;
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, Range<isize>> for Series<Var, C> {
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
    fn as_slice(&'a self, index: Range<isize>) -> SeriesSlice<'a, Var, C> {
        let mut index = index;
        let slice_min_pow = index.start;
        index.start -= self.min_pow;
        index.end -= self.min_pow;
        let index = (index.start as usize)..(index.end as usize);
        SeriesSlice{
            var: &self.var,
            min_pow: slice_min_pow,
            coeffs: &self.coeffs[index],
            zero: &self.zero,
        }
    }
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, RangeFrom<isize>> for Series<Var, C> {
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
    fn as_slice(&'a self, index: RangeFrom<isize>) -> SeriesSlice<'a, Var, C> {
        let mut index = index;
        let slice_min_pow = index.start;
        index.start -= self.min_pow;
        let index = (index.start as usize)..;
        SeriesSlice{
            var: &self.var,
            min_pow: slice_min_pow,
            coeffs: &self.coeffs[index],
            zero: &self.zero,
        }
    }
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, RangeFull> for Series<Var, C> {
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
    fn as_slice(&'a self, _index: RangeFull) -> SeriesSlice<'a, Var, C> {
        SeriesSlice{
            var: &self.var,
            min_pow: self.min_pow,
            coeffs: self.coeffs.as_slice(),
            zero: &self.zero,
        }
    }
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, RangeInclusive<isize>> for Series<Var, C> {
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
    fn as_slice(&'a self, index: RangeInclusive<isize>) -> SeriesSlice<'a, Var, C> {
        let (mut start, mut end) = index.into_inner();
        let slice_min_pow = start;
        start -= self.min_pow;
        end -= self.min_pow;
        let index = (start as usize)..=(end as usize);
        SeriesSlice{
            var: &self.var,
            min_pow: slice_min_pow,
            coeffs: &self.coeffs[index],
            zero: &self.zero,
        }
    }
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, RangeTo<isize>> for Series<Var, C> {
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
    fn as_slice(&'a self, index: RangeTo<isize>) -> SeriesSlice<'a, Var, C> {
        let mut index = index;
        index.end -= self.min_pow;
        let index = RangeTo{
            end: index.end as usize,
        };
        SeriesSlice{
            var: &self.var,
            min_pow: self.min_pow,
            coeffs: &self.coeffs[index],
            zero: &self.zero,
        }
    }
}

impl<'a, Var, C: Coeff> AsSlice<'a, Var, C, RangeToInclusive<isize>> for Series<Var, C> {
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
    fn as_slice(&'a self, index: RangeToInclusive<isize>) -> SeriesSlice<'a, Var, C> {
        let mut index = index;
        index.end -= self.min_pow;
        let index = ..=(index.end as usize);
        SeriesSlice{
            var: &self.var,
            min_pow: self.min_pow,
            coeffs: &self.coeffs[index],
            zero: &self.zero,
        }
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

/// Multiplicative inverse
pub trait MulInverse {
    type Output;

    fn mul_inverse(self) -> Self::Output;
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
        let idx = {
            let non_zero = self
                .coeffs
                .iter()
                .enumerate()
                .find(|&(_, c)| *c != C::from(0));
            match non_zero {
                Some((idx, _)) => Some(idx),
                _ => None,
            }
        };
        match idx {
            Some(idx) => {
                if idx > 0 {
                    self.min_pow += idx as isize;
                    self.coeffs.drain(0..idx);
                }
            }
            None => {
                self.min_pow += self.coeffs.len() as isize;
                self.coeffs = vec![];
            }
        }
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

// another dubious helper trait that only serves to prevent obscure
// compiler errors in rust 1.36.0
trait MulNaive<'a, Var, C: Coeff>{
    fn mul_naive(self, b: SeriesSlice<'a, Var, C>) -> Series<Var, C>;
}

impl <'a, 'b, Var, C: Coeff> MulNaive<'b, Var, C> for SeriesSlice<'a, Var, C>
where
    for<'c,'d> &'c C: Mul<&'d C, Output = C>,
    C: AddAssign + Mul<Output = C> + Clone,
    Var: Clone + PartialEq + fmt::Debug,
{
    fn mul_naive(self, b: SeriesSlice<'b, Var, C>) -> Series<Var, C> {
        let a = self;
        assert_eq!(a.var, b.var);
        let res_len = a.coeffs.len() + b.coeffs.len();
        let mut res_coeffs = Vec::with_capacity(res_len);
        for _ in 0..res_len {
            res_coeffs.push(C::from(0));
        }
        for (i, a) in a.coeffs.iter().enumerate() {
            for (j, b) in b.coeffs.iter().enumerate() {
                res_coeffs[i+j] += a.clone() * b.clone()
            }
        }
        Series::new(a.var.clone(), a.min_pow() + b.min_pow(), res_coeffs)
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
    <Series<Var, C> as ops::Ln>::Output: Mul<T>,
    <<Series<Var, C> as ops::Ln>::Output as std::ops::Mul<T>>::Output: Exp,
{
    type Output = <<<Series<Var, C> as ops::Ln>::Output as std::ops::Mul<T>>::Output as ops::Exp>::Output;

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

// logarithm for series starting with var^0
trait LnVarFree {
    type Output;
    fn ln_var_free(self) -> Self::Output;
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

/// View into a Laurent series
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct SeriesSlice<'a, Var, C: Coeff> {
    var: &'a Var,
    min_pow: isize,
    coeffs: &'a [C],
    zero: &'a C,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, Var, C: Coeff> std::marker::Copy for SeriesSlice<'a, Var, C>
{}

impl<'a, Var, C: Coeff> std::clone::Clone for SeriesSlice<'a, Var, C>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, Var, C: Coeff> SeriesSlice<'a, Var, C> {
    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    pub fn cutoff_pow(&self) -> isize {
        self.min_pow + (self.coeffs.len() as isize)
    }

    pub fn coeff(&self, pow: isize) -> Option<&C> {
        if pow < self.min_pow() {
            return Some(self.zero); // TODO this is a bad hack
        }
        if pow >= self.cutoff_pow() {
            return None;
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(&self.coeffs[idx])
    }

    pub fn iter(&self) -> Iter<C> {
        (self.min_pow..).zip(self.coeffs.iter())
    }

    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos + self.min_pow()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = SeriesSlice{
            var: self.var,
            min_pow: self.min_pow(),
            coeffs: lower,
            zero: self.zero
        };
        let upper = SeriesSlice{
            var: self.var,
            min_pow: pos,
            coeffs: upper,
            zero: self.zero
        };
        (lower, upper)
    }
}

impl<'a, Var, C: Coeff> Index<isize> for SeriesSlice<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index-self.min_pow) as usize]
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> SeriesSlice<'a, Var, C> {
    pub fn to_owned(&self) -> Series<Var, C> {
        Series::new(self.var.clone(), self.min_pow, self.coeffs.to_vec())
    }
}

impl<'a, Var, C> MulInverse for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    C: Coeff + SubAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
{
    type Output = Series<Var, C>;

    fn mul_inverse(self) -> Self::Output {
        let inv_min_pow = -self.min_pow;
        if self.coeffs.is_empty() {
            return Series::new(self.var.clone(), inv_min_pow, vec![]);
        }
        let a: Vec<_> =
            self.coeffs.iter().map(|c| c / &self.coeffs[0]).collect();
        let mut b = Vec::with_capacity(a.len());
        b.push(C::from(1));
        for n in 1..a.len() {
            let mut b_n = C::from(0);
            for i in 0..n {
                b_n -= &a[n - i] * &b[i];
            }
            b.push(b_n);
        }
        let inv_coeffs: Vec<_> =
            b.iter().map(|b| b / &self.coeffs[0]).collect();
        Series::new(self.var.clone(), inv_min_pow, inv_coeffs)
    }
}

impl<'a, Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for SeriesSlice<'a, Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if *coeff == C::from(0) {
                continue;
            }
            let cur_pow = self.min_pow() + i as isize;
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "({})", coeff)?;
            if cur_pow != 0 {
                write!(f, "*{}", self.var)?;
                if cur_pow != 1 {
                    write!(f, "^{}", cur_pow)?;
                }
            }
        }
        if !self.coeffs.is_empty() {
            write!(f, " + ")?;
        }
        write!(f, "O({}^{})", self.var, self.cutoff_pow())
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for SeriesSlice<'a, Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Series<Var, C>;

    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Series::new(self.var.clone(), self.min_pow, neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for SeriesSlice<'a, Var, C>
where
    Series<Var, C>: AddAssign<Rhs>,
{
    type Output = Series<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.to_owned();
        res += other;
        res
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

// Spurious trait, needed for rust 1.28
// direct implementation used to work in 1.24
trait AddAssignHelper<Var, C: Coeff> {
    fn truncate_cutoff_pow<'a>(&mut self, other: SeriesSlice<'a, Var, C>);
    fn add_overlap<'a>(&mut self, other: SeriesSlice<'a, Var, C>);
    fn num_leading<'a>(&mut self, other: SeriesSlice<'a, Var, C>) -> usize;
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

impl<'a, Var, C: Coeff, T> Sub<T> for SeriesSlice<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    Series<Var, C>: SubAssign<T>,
{
    type Output = Series<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.to_owned();
        res -= other;
        res
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

impl<'a, Var, C: Coeff, T> Mul<T> for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: MulAssign<T>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        let mut res = self.to_owned();
        res *= other;
        res
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

impl<'a, Var, C: Coeff, T> Div<T> for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: DivAssign<T>,
{
    type Output = Series<Var, C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = self.to_owned();
        res /= other;
        res
    }
}

impl<'a, Var, C: Coeff> Exp for SeriesSlice<'a, Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
    for<'b> C: MulAssign<&'b C>,
    Var: Clone,
    C: Clone + Div<Output = C> + Mul<Output = C> + AddAssign + Exp<Output = C>,
{
    type Output = Series<Var, C>;

    fn exp(self) -> Self::Output {
        Series::new(self.var.clone(), 0, self.exp_coeff())
    }
}

// TODO: understand why there is a compiler error when removing this trait
//       (with rust 1.28, ok with 1.24)
//       and just implement the method
trait ExpCoeff {
    type Output;
    fn exp_coeff(&self) -> Self::Output;
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

impl<'a, Var, C: Coeff> Ln for SeriesSlice<'a, Var, C>
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
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = &self.coeffs[0];
        // self.coeffs[0] = C::from(1);
        // for i in 1..self.coeffs.len() {
        //     self.coeffs[i] /= &c_k0;
        // }
        let a = &self.coeffs;
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

impl<'a, Var, C: Coeff, T> Pow<T> for SeriesSlice<'a, Var, C>
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tst_series() {
        let var = String::from("x");
        let min_pow = -10;
        let coeffs = vec![];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-11), Some(&0));
        assert_eq!(s.coeff(-10), None);

        let min_pow = -3;
        let coeffs = vec![1., 2., 3.];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-4), Some(&0.));
        assert_eq!(s.coeff(-3), Some(&1.));
        assert_eq!(s.coeff(-2), Some(&2.));
        assert_eq!(s.coeff(-1), Some(&3.));
        assert_eq!(s.coeff(0), None);

        let min_pow = -2;
        let coeffs = vec![0., 0., 3.];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow + 2);
        assert_eq!(s.coeff(-2), Some(&0.));
        assert_eq!(s.coeff(-1), Some(&0.));
        assert_eq!(s.coeff(0), Some(&3.));
        assert_eq!(s.coeff(1), None);

        let s = Series::new(var.clone(), -2, vec![0., 0., 1.]);
        let t = Series::new(var.clone(), 0, vec![1.]);
        assert_eq!(s, t);

        let s = Series::new(var.clone(), -3, vec![0., 0., 0.]);
        let t = Series::new(var.clone(), 0, vec![]);
        assert_eq!(s, t);
    }

    #[test]
    fn tst_series_with_cutoff() {
        let s = Series::with_cutoff("x", -10, 1, Vec::<i32>::new());
        let t = Series::new("x", 1, vec![]);
        assert_eq!(s, t);

        let s = Series::with_cutoff("x", 0, 5, vec![1, 2, 3]);
        let t = Series::new("x", 0, vec![1, 2, 3, 0, 0]);
        assert_eq!(s, t);

        let s = Series::with_cutoff("x", 0, 2, vec![1, 2, 3]);
        let t = Series::new("x", 0, vec![1, 2]);
        assert_eq!(s, t);
    }
    #[test]
    #[should_panic]
    fn tst_bad_cutoff() {
        let _ = Series::with_cutoff("x", 0, -2, vec![1, 2, 3]);
    }

    #[test]
    fn tst_display() {
        // let s = Series::new("x", -10, vec!());
        // assert_eq!(format!("{}", s), "O(x^-10)");
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-3 + (-3)*x^-1 + O(x^0)");
        let s = Series::new("x", -1, vec![1., 2., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-1 + (2) + (-3)*x + O(x^2)");
    }

    #[test]
    fn tst_neg() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![-1., 0., 3.]);
        assert_eq!(res, -&s);
        assert_eq!(res, -s);
    }

    #[test]
    fn tst_add() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![2., 0., -6.]);
        assert_eq!(res, &s + &s);
        assert_eq!(res, &s + s.clone());
        assert_eq!(res, s.clone() + &s);
        assert_eq!(res, s.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        assert_eq!(s, &s + &t);
        assert_eq!(s, &t + &s);
        assert_eq!(s, &s + t.clone());
        assert_eq!(s, &t + s.clone());
        assert_eq!(s, s.clone() + &t);
        assert_eq!(s, t.clone() + &s);
        assert_eq!(s, s.clone() + t.clone());
        assert_eq!(s, t.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -3, vec![-1., 0., 3.]);
        let res = Series::new("x", 0, vec![]);
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());
    }

    #[test]
    fn tst_add_assign() {
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![2., 0., -6.]);
        s += s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -1, vec![3., 4., 5.]);
        let t = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = s.clone();
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        s += &t;
        assert_eq!(s, res);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(s, res);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -3, vec![-1., 0., 3.]);
        let res = Series::new("x", 0, vec![]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_sub() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![]);
        assert_eq!(res, &s - &s);
        assert_eq!(res, &s - s.clone());
        assert_eq!(res, s.clone() - &s);
        assert_eq!(res, s.clone() - s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![-3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        assert_eq!(res, &s - &t);
        assert_eq!(res, &s - t.clone());
        assert_eq!(res, s.clone() - &t);
        assert_eq!(res, s - t);

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        assert_eq!(s, &s - &t);
        assert_eq!(s, &s - t.clone());
        assert_eq!(s, s.clone() - &t);
        assert_eq!(s, s.clone() - t);
    }

    #[test]
    fn tst_sub_assign() {
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![]);
        s -= s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![-3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = s.clone();
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_mul() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -6, vec![1., 0., -6.]);
        assert_eq!(res, &s * &s);
        assert_eq!(res, &s * s.clone());
        assert_eq!(res, s.clone() * &s);
        assert_eq!(res, s.clone() * s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5., 7.]);
        let res = Series::new("x", -4, vec![3., 4., -4.]);
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);

        let s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, s.clone() * &t);
        assert_eq!(res, t.clone() * &s);
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);
    }

    #[test]
    fn tst_mul_naive() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -6, vec![1., 0., -6.]);
        let square = s.as_slice().mul_naive(s.as_slice());
        assert_eq!(square.cutoff_pow(), 0);
        let square = square + Series::new("x", -3, vec![]);
        assert_eq!(square, res);

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5., 7.]);
        let res = Series::new("x", -4, vec![3., 4., -4.]);
        let prod = s.as_slice().mul_naive(t.as_slice());
        assert_eq!(prod.cutoff_pow(), 3);
        let prod = prod + Series::new("x", -1, vec![]);
        assert_eq!(prod, res);

        let s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        let prod = s.as_slice().mul_naive(t.as_slice());
        assert_eq!(prod.cutoff_pow(), 6);
        let prod = prod + Series::new("x", 3, vec![]);
        assert_eq!(prod, res);
    }

    #[test]
    fn tst_mul_assign() {
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s *= s.clone();
        let res = Series::new("x", -6, vec![1., 0., -6.]);
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5., 7.]);
        s *= &t;
        let res = Series::new("x", -4, vec![3., 4., -4.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s *= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        s *= &t;
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        s *= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_div() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, &s / &s);
        assert_eq!(res, &s / s.clone());
        assert_eq!(res, s.clone() / &s);
        assert_eq!(res, s.clone() / s.clone());

        // disabled for floating-point inaccuracy
        // let s = Series::new("x", -3, vec!(1.,0.,-3.));
        // let t = Series::new("x", -1, vec!(3., 4., 5., 7.));
        // let res = Series::new("x", -2, vec!(1./3.,-4./9.,-26./27.));
        // assert_eq!(res, &s / &t);
        // assert_eq!(res, &t / &s);
        // assert_eq!(res, s / t);

        let s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        let res = Series::new("x", -6, vec![1., 14., 43.]);
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);

        let s = Series::new("x", 1, vec![1., 7., -3.]);
        let t = Series::new("x", 5, vec![]);
        let res = Series::new("x", -4, vec![]);
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);
    }

    #[test]
    fn tst_div_assign() {
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s /= s.clone();
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        s /= &t;
        let res = Series::new("x", -6, vec![1., 14., 43.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        s /= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", 1, vec![1., 7., -3.]);
        let t = Series::new("x", 5, vec![]);
        s /= &t;
        let res = Series::new("x", -4, vec![]);
        assert_eq!(res, s);
        let mut s = Series::new("x", 1, vec![1., 7., -3.]);
        s /= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_var() {
        let _ = Series::new(String::from("x"), -3, vec![1., 0., -3.]);
        let _ = Series::new('j', -3, vec![1., 0., -3.]);
        let _ = Series::new(8, -3, vec![1., 0., -3.]);
    }

    #[derive(Debug, Clone, PartialEq)]
    struct Mystr<'a>(&'a str);

    impl<'a> From<Mystr<'a>> for f64 {
        fn from(_s: Mystr<'a>) -> f64 {
            panic!("can't turn str to f64")
        }
    }

    #[test]
    fn tst_ln() {
        let s = Series::new(Mystr("x"), 0, vec![1., 7., -3.]);
        let res = Series::new(Mystr("x"), 1, vec![7., -55. / 2.]);
        assert_eq!(res, (&s).ln());
        assert_eq!(res, s.ln());

        let s = Series::new(Mystr("x"), 0, vec![4., 7., -3.]);
        let res =
            Series::new(Mystr("x"), 0, vec![4_f64.ln(), 7. / 4., -73. / 32.]);
        assert_eq!(res, (&s).ln());
        assert_eq!(res, s.ln());
    }

    #[test]
    fn tst_exp() {
        let s = Series::new("x", 1, vec![7., -3.]);
        let res = Series::new("x", 0, vec![1., 7., 43. / 2.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 2, vec![0.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 0, vec![5., 11., -7.]);
        let e5 = 5_f64.exp();
        let res = Series::new("x", 0, vec![e5, e5 * 11., e5 * 107. / 2.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());
    }

    #[test]
    fn tst_pow() {
        let base = Series::new(Mystr("x"), 0, vec![1., 7., 0.]);
        let exp = Series::new(Mystr("x"), -1, vec![1., -5., 43.]);
        let e7 = 7_f64.exp();
        let res = Series::new(Mystr("x"), 0, vec![e7, -119. / 2. * e7]);
        assert_eq!(res, (&base).pow(&exp));
        assert_eq!(res, (&base).pow(exp.clone()));
        assert_eq!(res, base.clone().pow(&exp));
        assert_eq!(res, base.pow(exp));

        let base = Series::new(Mystr("x"), 0, vec![2., 7., 0.]);
        let exp = Series::new(Mystr("x"), 0, vec![3., -5., 11.]);
        // rescale result so we can use round and still get decent precision
        let rescale = Series::new(Mystr("x"), 0, vec![1e13, 0., 0., 0.]);
        let test = &rescale * &base.pow(exp);
        let ln2 = 2_f64.ln();
        let res = Series::new(
            Mystr("x"),
            0,
            vec![8., 84. - 40. * ln2, 154. + ln2 * (-332. + 100. * ln2)],
        );
        let res = rescale * res;
        assert_eq!(res.min_pow(), test.min_pow());
        assert_eq!(res.cutoff_pow(), test.cutoff_pow());
        for i in res.min_pow()..res.cutoff_pow() {
            assert_eq!(
                res.coeff(i).unwrap().round(),
                test.coeff(i).unwrap().round()
            );
        }
    }

    #[test]
    fn tst_scalar() {
        let s = Series::new("x", -3, vec![1., 0., -2.]);
        let res = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(res, &s / 2.);
        let mut s = s;
        s /= 2.;
        assert_eq!(res, s);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -3, vec![1., 0., -2.]);
        assert_eq!(res, &s * 2.);
        let mut s = s;
        s *= 2.;
        assert_eq!(res, s);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(s, &s + 2.);
        let s = Series::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -2, vec![1. / 2., 0., 1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Series::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", 0, vec![2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Series::new("x", 0, vec![-2., 0., -1.]);
        let res = Series::new("x", 2, vec![-1.]);
        assert_eq!(res, s + 2.);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(s, &s - 2.);
        let s = Series::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -2, vec![1. / 2., 0., -3.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Series::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", 0, vec![-2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Series::new("x", 0, vec![2., 0., -1.]);
        let res = Series::new("x", 2, vec![-1.]);
        assert_eq!(res, s - 2.);

        let base = Series::new(Mystr("x"), 0, vec![1., 7., 0.]);
        assert_eq!(base, (&base).pow(1.));
        assert_eq!(base, (&base).pow(&1.));
        let res = Series::new(Mystr("x"), 0, vec![1., 21., 147.]);
        assert_eq!(res, (&base).pow(3.));
        assert_eq!(res, base.pow(&3.));
    }

    #[test]
    #[should_panic]
    fn tst_bad_add() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s + t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_sub() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s - t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_mul() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s * t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_div() {
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s / t;
    }
}
