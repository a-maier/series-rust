use crate::ops::{Exp, Ln, Pow};
use crate::util::trim_start_zero;
use crate::{Coeff, IntoIter, Iter};
use crate::{Series, anon_series_slice::AnonSeriesSlice};
use crate::{SeriesParts, traits::*};

use std::cmp::{Ordering, min};
use std::convert::From;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};

/// Laurent series in a single anonymous variable up to some power
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct AnonSeries<C: Coeff> {
    pub(crate) min_pow: isize,
    pub(crate) coeffs: Vec<C>,
}

impl<C: Coeff> AnonSeries<C> {
    /// Create a new series
    ///
    /// # Example
    ///
    /// This creates a series starting with the power -1 with
    /// coefficients 1, 2, 3. In other words, the series x^-1 + 2 +
    /// 3*x + O(x^2).
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// ```
    pub fn new(min_pow: isize, coeffs: Vec<C>) -> Self {
        let mut res = AnonSeries { min_pow, coeffs };
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
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::with_cutoff(-1..5, vec![1, 2, 3]);
    /// ```
    pub fn with_cutoff(powers: Range<isize>, mut coeffs: Vec<C>) -> Self {
        let min_pow = powers.start;
        let cutoff_pow = powers.end;
        if cutoff_pow < min_pow {
            return AnonSeries::new(cutoff_pow, vec![]);
        }
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
        AnonSeries::new(min_pow, coeffs)
    }

    /// Turn into a series with a named expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]).in_var("x");
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn in_var<Var>(self, var: Var) -> Series<Var, C> {
        Series { series: self, var }
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// assert_eq!(s.cutoff_pow(), 2);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.as_slice(..).cutoff_pow()
    }

    /// Get the number of known coefficients in the series.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::with_cutoff(-1..5, vec!(1,2,3));
    /// assert_eq!(s.len(), 6);
    /// // This holds true for any series
    /// assert_eq!(s.len(), (s.cutoff_pow() - s.min_pow()) as usize);
    /// ```
    pub fn len(&self) -> usize {
        self.as_slice(..).len()
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// let mut iter = s.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<'_, C> {
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
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-1, vec![1,2,3]);
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
        if pow < self.min_pow() {
            return self.apply_at_new_front(pow, f);
        }
        let index = (pow - self.min_pow()) as usize;
        f(&mut self.coeffs[index]);
        if index == 0 && self.coeffs[0].is_zero() {
            self.trim()
        }
    }

    fn apply_at_new_front<F: FnOnce(&mut C)>(&mut self, pow: isize, f: F) {
        let mut c = C::zero();
        f(&mut c);
        if c.is_zero() {
            return;
        }
        let nnew = (self.min_pow() - pow) as usize;
        let new_len = self.coeffs.len() + nnew;
        self.coeffs.push(c);
        self.coeffs.resize_with(new_len, || C::zero());
        self.coeffs.rotate_right(nnew);
        self.min_pow = pow;
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
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// s.for_each(|_, c| *c *= *c);
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&4));
    /// assert_eq!(s.coeff(1), Some(&9));
    /// assert_eq!(s.coeff(2), Some(&16));
    /// ```
    pub fn for_each<F>(&mut self, mut f: F)
    where
        F: FnMut(isize, &mut C),
    {
        let min_pow = self.min_pow;
        for (n, c) in &mut self.coeffs.iter_mut().enumerate() {
            f(min_pow + n as isize, c)
        }
        self.trim();
    }

    /// Transform all coefficients
    ///
    /// `f(p, c)` is applied to each term, where `p` is the power of
    /// the variable and `c` a coefficient. `p` takes all values in
    /// the range `min_pow()..cutoff_pow()`.
    ///
    /// # Example
    ///
    /// Replace each coefficient by its square
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let s = s.map(|_, c| c * c);
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&4));
    /// assert_eq!(s.coeff(1), Some(&9));
    /// assert_eq!(s.coeff(2), Some(&16));
    /// ```
    pub fn map<D: Coeff, F>(self, mut f: F) -> AnonSeries<D>
    where
        F: FnMut(isize, C) -> D,
    {
        let min_pow = self.min_pow;
        let coeffs = self.into_iter().map(|(pow, c)| f(pow, c)).collect();
        AnonSeries::new(min_pow, coeffs)
    }
}

impl<C: 'static + Coeff + Send + Sync> AnonSeries<C> {
    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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

impl<'a, C: Coeff> AsSlice<Range<isize>> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

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
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(0..2);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(self, r: Range<isize>) -> Self::Output {
        let start = (r.start - self.min_pow()) as usize;
        let end = (r.end - self.min_pow()) as usize;
        AnonSeriesSlice::new(r.start, &self.coeffs[start..end])
    }
}

impl<'a, C: Coeff> AsSlice<RangeInclusive<isize>> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

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
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(0..=1);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(self, r: RangeInclusive<isize>) -> Self::Output {
        let (start, end) = r.into_inner();
        let ustart = (start - self.min_pow()) as usize;
        let end = (end - self.min_pow()) as usize;
        AnonSeriesSlice::new(start, &self.coeffs[ustart..=end])
    }
}

impl<'a, C: Coeff> AsSlice<RangeToInclusive<isize>> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(..=1);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(self, r: RangeToInclusive<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        AnonSeriesSlice::new(self.min_pow, &self.coeffs[..=end])
    }
}

impl<'a, C: Coeff> AsSlice<RangeFrom<isize>> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the lower bound is smaller than the leading power.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(0..);
    /// assert_eq!(t.min_pow(), 0);
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert!(std::panic::catch_unwind(|| t[-1]).is_err());
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(self, r: RangeFrom<isize>) -> Self::Output {
        let start = (r.start - self.min_pow()) as usize;
        AnonSeriesSlice::new(r.start, &self.coeffs[start..])
    }
}

impl<'a, C: Coeff> AsSlice<RangeTo<isize>> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

    /// A slice of the series truncated to the given range of powers.
    ///
    /// # Panics
    ///
    /// Panics if the upper bound is at least as big as the cut-off power.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(..2);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), 2);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert!(std::panic::catch_unwind(|| t[2]).is_err());
    /// ```
    fn as_slice(self, r: RangeTo<isize>) -> Self::Output {
        let end = (r.end - self.min_pow()) as usize;
        AnonSeriesSlice::new(self.min_pow, &self.coeffs[..end])
    }
}

impl<'a, C: Coeff> AsSlice<RangeFull> for &'a AnonSeries<C> {
    type Output = AnonSeriesSlice<'a, C>;

    /// A slice containing the complete series.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3, 4]);
    /// let t = s.as_slice(..);
    /// assert_eq!(t.min_pow(), s.min_pow());
    /// assert_eq!(t.cutoff_pow(), s.cutoff_pow());
    /// assert_eq!(t[-1], s[-1]);
    /// assert_eq!(t[0], s[0]);
    /// assert_eq!(t[1], s[1]);
    /// assert_eq!(t[2], s[2]);
    /// ```
    fn as_slice(self, r: RangeFull) -> Self::Output {
        AnonSeriesSlice::new(self.min_pow, &self.coeffs[r])
    }
}

impl<C: Coeff> Index<isize> for AnonSeries<C> {
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
    /// # use series::anon_series::AnonSeries;
    /// use series::AsSlice;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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

impl<C: Coeff> std::iter::IntoIterator for AnonSeries<C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the series coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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

impl<C: Coeff + SubAssign> MulInverse for &AnonSeries<C>
where
    C: Coeff + SubAssign,
    for<'c> &'c C: Div<Output = C> + Mul<Output = C>,
{
    type Output = AnonSeries<C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::MulInverse;
    /// let s = AnonSeries::new(-1, vec![1., 2., 3.]);
    /// let s_inv = (&s).mul_inverse();
    /// let one = AnonSeries::new(0, vec![1., 0., 0.]);
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        self.as_slice(..).mul_inverse()
    }
}

impl<C: Coeff + SubAssign> MulInverse for AnonSeries<C>
where
    for<'a> &'a AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    type Output = AnonSeries<C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// use series::MulInverse;
    /// let s = AnonSeries::new(-1, vec![1., 2., 3.]);
    /// let s_inv = s.clone().mul_inverse();
    /// let one = AnonSeries::new(0, vec![1., 0., 0.]);
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self).mul_inverse()
    }
}

impl<C: Coeff> AnonSeries<C> {
    fn trim(&mut self) {
        self.min_pow += trim_start_zero(&mut self.coeffs) as isize;
    }
}

impl<C: Coeff + Neg<Output = C>> Neg for AnonSeries<C> {
    type Output = AnonSeries<C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let minus_s = AnonSeries::new(-3, vec![-1., 0., 3.]);
    /// assert_eq!(-s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        AnonSeries::new(self.min_pow, neg_coeff)
    }
}

impl<C: Coeff> Neg for &AnonSeries<C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = AnonSeries<C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let minus_s = AnonSeries::new(-3, vec![-1., 0., 3.]);
    /// assert_eq!(-&s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, C: Coeff + Clone> AddAssign<&'a AnonSeries<C>> for AnonSeries<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let t = AnonSeries::new(-1, vec![3., 4., 5.]);
    /// let res = AnonSeries::new(-3, vec![1., 0., 0.]);
    /// s += &t;
    /// assert_eq!(res, s);
    /// ```
    fn add_assign(&mut self, other: &'a AnonSeries<C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, C: Coeff + Clone> AddAssign<AnonSeriesSlice<'a, C>> for AnonSeries<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: AnonSeriesSlice<'a, C>) {
        if self.cutoff_pow() > other.cutoff_pow() {
            self.truncate_cutoff_pow(other);
        }
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

impl<C: AddAssign + Coeff> AddAssign for AnonSeries<C> {
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let t = AnonSeries::new(-1, vec![3., 4., 5.]);
    /// let res = AnonSeries::new(-3, vec![1., 0., 0.]);
    /// s += t;
    /// assert_eq!(res, s);
    /// ```
    fn add_assign(&mut self, mut other: AnonSeries<C>) {
        match self.cutoff_pow().cmp(&other.cutoff_pow()) {
            Ordering::Less => other.truncate_cutoff_pow(self.as_slice(..)),
            Ordering::Equal => {}
            Ordering::Greater => self.truncate_cutoff_pow(other.as_slice(..)),
        }
        if self.min_pow() > other.min_pow() {
            std::mem::swap(self, &mut other);
        }
        let offset = other.min_pow() - self.min_pow();
        let lhs = self.coeffs.iter_mut().skip(offset as usize);
        let rhs = other.coeffs.into_iter();
        debug_assert_eq!(lhs.size_hint(), rhs.size_hint());
        for (lhs, rhs) in lhs.zip(rhs) {
            lhs.add_assign(rhs);
        }
        self.trim();
    }
}

impl<C: Coeff + Clone, Rhs> Add<Rhs> for &AnonSeries<C>
where
    AnonSeries<C>: AddAssign<Rhs>,
{
    type Output = AnonSeries<C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<C: Coeff + Clone, Rhs> Add<Rhs> for AnonSeries<C>
where
    AnonSeries<C>: AddAssign<Rhs>,
{
    type Output = AnonSeries<C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, C: Coeff> SubAssign<&'a AnonSeries<C>> for AnonSeries<C>
where
    for<'c> &'c AnonSeries<C>: Neg<Output = AnonSeries<C>>,
    AnonSeries<C>: AddAssign<AnonSeries<C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let res = AnonSeries::new(0, vec![]);
    /// s -= &s.clone();
    /// assert_eq!(res, s);
    /// ```
    fn sub_assign(&mut self, other: &'a AnonSeries<C>) {
        *self += -other;
    }
}

impl<'a, C: Coeff> SubAssign<AnonSeriesSlice<'a, C>> for AnonSeries<C>
where
    for<'c> AnonSeriesSlice<'c, C>: Neg<Output = AnonSeries<C>>,
    AnonSeries<C>: AddAssign<AnonSeries<C>>,
{
    fn sub_assign(&mut self, other: AnonSeriesSlice<'a, C>) {
        *self += -other;
    }
}

impl<C: Coeff> SubAssign<AnonSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: AddAssign + Neg<Output = AnonSeries<C>>,
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// let res = AnonSeries::new(0, vec![]);
    /// s -= s.clone();
    /// assert_eq!(res, s);
    /// ```
    fn sub_assign(&mut self, other: AnonSeries<C>) {
        *self += -other;
    }
}

// TODO: somehow make addition symmetric?
impl<C: Coeff, T> Sub<T> for &AnonSeries<C>
where
    AnonSeries<C>: Clone + SubAssign<T>,
{
    type Output = AnonSeries<C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res -= other;
        res
    }
}

impl<C: Coeff, T> Sub<T> for AnonSeries<C>
where
    AnonSeries<C>: SubAssign<T>,
{
    type Output = AnonSeries<C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, C: Coeff + Clone + AddAssign> MulAssign<&'a AnonSeries<C>>
    for AnonSeries<C>
where
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// s *= &s.clone();
    /// let res = AnonSeries::new(-6, vec![1., 0., -6.]);
    /// assert_eq!(res, s);
    /// ```
    fn mul_assign(&mut self, other: &'a AnonSeries<C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, C> MulAssign<AnonSeriesSlice<'a, C>> for AnonSeries<C>
where
    for<'b> &'b C: Mul<Output = C>,
    C: MulAssign<&'a C> + Coeff + Clone + AddAssign,
{
    fn mul_assign(&mut self, other: AnonSeriesSlice<'a, C>) {
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

impl<C: Coeff> MulAssign for AnonSeries<C>
where
    for<'a> AnonSeries<C>: MulAssign<&'a AnonSeries<C>>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// s *= &s.clone();
    /// let res = AnonSeries::new(-6, vec![1., 0., -6.]);
    /// assert_eq!(res, s);
    /// ```
    fn mul_assign(&mut self, other: AnonSeries<C>) {
        *self *= &other
    }
}

// TODO: somehow make multiplication symmetric?
impl<C: Coeff> Mul for AnonSeries<C>
where
    AnonSeries<C>: MulAssign,
{
    type Output = AnonSeries<C>;

    fn mul(mut self, other: AnonSeries<C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a AnonSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: MulAssign<AnonSeriesSlice<'a, C>>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: &'a AnonSeries<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, C: Coeff> Mul<AnonSeriesSlice<'a, C>> for AnonSeries<C>
where
    AnonSeries<C>: MulAssign<AnonSeriesSlice<'a, C>>,
{
    type Output = AnonSeries<C>;

    fn mul(mut self, other: AnonSeriesSlice<'a, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff> Mul<C> for AnonSeries<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = AnonSeries<C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a C> for AnonSeries<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = AnonSeries<C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff, T> Mul<T> for &'a AnonSeries<C>
where
    AnonSeriesSlice<'a, C>: Mul<T, Output = AnonSeries<C>>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, C: Coeff + SubAssign> DivAssign<&'a AnonSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// s /= &s.clone();
    /// let res = AnonSeries::new(0, vec![1., 0., 0.]);
    /// assert_eq!(res, s);
    /// ```
    fn div_assign(&mut self, other: &'a AnonSeries<C>) {
        self.div_assign(other.as_slice(..));
    }
}

impl<C: Coeff + SubAssign> DivAssign for AnonSeries<C>
where
    AnonSeries<C>: MulAssign + MulInverse<Output = AnonSeries<C>>,
    for<'a> &'a C: Div<Output = C> + Mul<Output = C>,
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let mut s = AnonSeries::new(-3, vec![1., 0., -3.]);
    /// s /= s.clone();
    /// let res = AnonSeries::new(0, vec![1., 0., 0.]);
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(&mut self, other: AnonSeries<C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, C: Coeff + SubAssign> DivAssign<AnonSeriesSlice<'a, C>>
    for AnonSeries<C>
where
    AnonSeries<C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    fn div_assign(&mut self, other: AnonSeriesSlice<'a, C>) {
        *self *= other.mul_inverse();
    }
}

impl<C: Coeff, T> Div<T> for &AnonSeries<C>
where
    AnonSeries<C>: Clone + DivAssign<T>,
{
    type Output = AnonSeries<C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res /= other;
        res
    }
}

impl<C: Coeff, T> Div<T> for AnonSeries<C>
where
    AnonSeries<C>: DivAssign<T>,
{
    type Output = AnonSeries<C>;

    fn div(mut self, other: T) -> Self::Output {
        self /= other;
        self
    }
}

impl<C: Coeff> Exp for AnonSeries<C>
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
        AnonSeries::new(0, coeff)
    }
}

impl<C: Coeff> Exp for &AnonSeries<C>
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
    type Output = AnonSeries<C>;

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

impl<C: Coeff> Ln for AnonSeries<C>
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

impl<C: Coeff> Ln for &AnonSeries<C>
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
    type Output = AnonSeries<C>;

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

impl<C: Coeff, T> Pow<T> for AnonSeries<C>
where
    AnonSeries<C>: Ln,
    <AnonSeries<C> as Ln>::Output: Mul<T>,
    <<AnonSeries<C> as Ln>::Output as std::ops::Mul<T>>::Output: Exp,
{
    type Output = <<<AnonSeries<C> as Ln>::Output as std::ops::Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

impl<C: Coeff, T> Pow<T> for &AnonSeries<C>
where
    for<'b> AnonSeriesSlice<'b, C>: Ln<Output = AnonSeries<C>>,
    AnonSeries<C>: Mul<T>,
    <AnonSeries<C> as Mul<T>>::Output: Exp,
{
    type Output = <<AnonSeries<C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        self.as_slice(..).pow(exponent)
    }
}

macro_rules! impl_add_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, C: Coeff> AddAssign<$rhs> for AnonSeries<C>
            where
                C: AddAssign<$rhs>,
            {
                fn add_assign(&mut self, rhs: $rhs) {
                    if self.cutoff_pow() <= 0 || rhs.is_zero() {
                        return;
                    }
                    if self.min_pow() > 0 {
                        self.coeffs.splice(
                            0..0,
                            std::iter::repeat_with(|| C::zero())
                                .take(self.min_pow() as usize)
                        );
                        self.min_pow = 0;
                    }
                    let idx = (-self.min_pow()) as usize;
                    self.coeffs[idx].add_assign(rhs);
                    if self.min_pow() == 0 {
                        self.trim()
                    }
                }
            }
        )*
    };
}

impl_add_assign_const!(C, &'a C);

macro_rules! impl_sub_assign_const {
    ($($rhs:ty), *) => {
        $(
            impl<'a, C: Coeff> SubAssign<$rhs> for AnonSeries<C>
            where
                C: SubAssign<$rhs>,
            {
                fn sub_assign(&mut self, rhs: $rhs) {
                    if self.cutoff_pow() <= 0 || rhs.is_zero() {
                        return;
                    }
                    if self.min_pow() > 0 {
                        self.coeffs.splice(
                            0..0,
                            std::iter::repeat_with(|| C::zero())
                                .take(self.min_pow() as usize)
                        );
                        self.min_pow = 0;
                    }
                    let idx = (-self.min_pow()) as usize;
                    self.coeffs[idx].sub_assign(rhs);
                    if self.min_pow() == 0 {
                        self.trim()
                    }
                }
            }
        )*
    };
}

impl_sub_assign_const!(C, &'a C);

impl<'a, C: Coeff> MulAssign<&'a C> for AnonSeries<C>
where
    C: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff *= rhs
        }
    }
}

impl<C: Coeff> MulAssign<C> for AnonSeries<C>
where
    for<'a> AnonSeries<C>: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: C) {
        *self *= &rhs
    }
}

impl<'a, C: Coeff> DivAssign<&'a C> for AnonSeries<C>
where
    C: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a C) {
        for coeff in &mut self.coeffs {
            *coeff /= rhs
        }
    }
}

impl<C: Coeff> DivAssign<C> for AnonSeries<C>
where
    for<'a> AnonSeries<C>: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: C) {
        *self /= &rhs
    }
}

impl<'a, C: Coeff + Clone> From<AnonSeriesSlice<'a, C>> for AnonSeries<C> {
    fn from(s: AnonSeriesSlice<'a, C>) -> Self {
        AnonSeries::new(s.min_pow, s.coeffs.to_vec())
    }
}

impl<C: Coeff, Var> From<Series<Var, C>> for AnonSeries<C> {
    fn from(source: Series<Var, C>) -> Self {
        let SeriesParts {
            var: _,
            min_pow,
            coeffs,
        } = source.into();
        Self { min_pow, coeffs }
    }
}

/// Data parts of a series
///
/// # Example
///
/// ```rust
/// // destructure a series
/// # use series::anon_series::{AnonSeries, AnonSeriesParts};
/// let s = AnonSeries::new(-1, vec![1,2,3]);
/// let AnonSeriesParts{min_pow, coeffs} = s.into();
/// assert_eq!(min_pow, -1);
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct AnonSeriesParts<C> {
    pub min_pow: isize,
    pub coeffs: Vec<C>,
}

impl<C: Coeff> From<AnonSeries<C>> for AnonSeriesParts<C> {
    fn from(s: AnonSeries<C>) -> Self {
        AnonSeriesParts {
            min_pow: s.min_pow,
            coeffs: s.coeffs,
        }
    }
}

impl<C: Coeff> From<AnonSeriesParts<C>> for AnonSeries<C> {
    fn from(parts: AnonSeriesParts<C>) -> Self {
        AnonSeries::new(parts.min_pow, parts.coeffs)
    }
}

impl<C: Coeff> LnVarFree for AnonSeries<C>
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
        AnonSeries::new(0, b)
    }
}

impl<C: Coeff> AnonSeries<C> {
    /// Calculate series to some integer power
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::anon_series::AnonSeries;
    /// let s = AnonSeries::new(-1, vec![1.,3.,7.]);
    /// let s_to_minus_5 = AnonSeries::new(5, vec![1.,-15.,100.]);
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
        AnonSeries<C>: Mul<C, Output = Self>
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

        AnonSeries::new(new_min_pow, pow.coeffs)
    }
}

impl<C: Coeff> AnonSeries<C> {
    fn truncate_cutoff_pow(&mut self, other: AnonSeriesSlice<'_, C>) {
        debug_assert!(other.cutoff_pow() < self.cutoff_pow());
        let to_remove = (self.cutoff_pow() - other.cutoff_pow()) as usize;
        if to_remove <= self.coeffs.len() {
            self.coeffs.truncate(self.coeffs.len() - to_remove);
        } else {
            self.coeffs.clear();
            self.min_pow = other.cutoff_pow()
        }
        debug_assert!(other.cutoff_pow() == self.cutoff_pow());
    }
}

impl<C: Coeff> AnonSeries<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_overlap(&mut self, other: AnonSeriesSlice<'_, C>) {
        let offset = self.min_pow();
        for (i, c) in self.coeffs.iter_mut().enumerate() {
            let power = offset + i as isize;
            if let Some(coeff) = other.try_coeff(power) {
                *c += coeff
            }
        }
    }

    fn num_leading(&mut self, other: AnonSeriesSlice<'_, C>) -> usize {
        min(
            (self.min_pow() - other.min_pow()) as usize,
            other.coeffs.len(),
        )
    }
}

impl<C: Coeff> ExpCoeff for AnonSeries<C>
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

impl<C: Coeff> ExpCoeff for AnonSeriesSlice<'_, C>
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
            let exp_a_0 = self
                .try_coeff(0)
                .cloned()
                .unwrap_or_else(|| C::zero())
                .exp();
            for b_n in &mut b {
                *b_n *= &exp_a_0;
            }
        }
        b
    }
}
