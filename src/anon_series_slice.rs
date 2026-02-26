use crate::ops::{Exp, Ln, Pow};
use crate::traits::{AsSlice, ExpCoeff, MulInverse};
use crate::util::trim_slice_start_zero;
use crate::zero_ref::zero_ref;
use crate::{Coeff, Iter, SeriesSlice, anon_series::AnonSeries};

use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// View into a Laurent series with an anonymous variable
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct AnonSeriesSlice<'a, C: Coeff> {
    pub(crate) min_pow: isize,
    pub(crate) coeffs: &'a [C],
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<C: Coeff> std::marker::Copy for AnonSeriesSlice<'_, C> {}

impl<C: Coeff> std::clone::Clone for AnonSeriesSlice<'_, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C: Coeff> AnonSeriesSlice<'a, C> {
    pub(super) fn new(min_pow: isize, coeffs: &'a [C]) -> Self {
        let mut res = AnonSeriesSlice { min_pow, coeffs };
        res.trim();
        res
    }

    fn trim(&mut self) {
        let removed = trim_slice_start_zero(&mut self.coeffs);
        self.min_pow += removed as isize;
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// assert_eq!(s.as_slice(..).min_pow(), -1);
    /// assert_eq!(s.as_slice(0..).min_pow(), 0);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    /// Get the power of the expansion variable where the slice is
    /// truncated.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// assert_eq!(s.as_slice(..).cutoff_pow(), 2);
    /// assert_eq!(s.as_slice(..1).cutoff_pow(), 1);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.min_pow + (self.coeffs.len() as isize)
    }

    /// Get the number of known coefficients in the series.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    /// let s = AnonSeries::with_cutoff(-1..5, vec![1, 2, 3]);
    /// assert_eq!(s.as_slice(..).len(), 6);
    /// // This holds true for any series
    /// assert_eq!(s.as_slice(..).len(), (s.cutoff_pow() - s.min_pow()) as usize);
    /// ```
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// let slice = s.as_slice(..);
    /// let mut iter = slice.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<'a, C> {
        (self.min_pow..).zip(self.coeffs.iter())
    }

    /// Split a series slice into two at the given power
    /// of the expansion variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// let (lower, upper) = s.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), -1);
    /// assert_eq!(upper.min_pow(), 0);
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos - self.min_pow()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = AnonSeriesSlice {
            min_pow: self.min_pow(),
            coeffs: lower,
        };
        let upper = AnonSeriesSlice {
            min_pow: pos,
            coeffs: upper,
        };
        (lower, upper)
    }

    /// Turn into a slice with a named expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
    /// let s = s.as_slice(..).in_var(&"x");
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn in_var<Var>(self, var: &'a Var) -> SeriesSlice<'a, Var, C> {
        SeriesSlice { var, series: self }
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
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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
        if pow >= self.min_pow() && pow < self.cutoff_pow() {
            Some(self.coeff_in_range(pow))
        } else {
            None
        }
    }

    fn coeff_in_range(&self, pow: isize) -> &'a C {
        debug_assert!(pow >= self.min_pow());
        debug_assert!(pow < self.cutoff_pow());
        let idx = (pow - self.min_pow()) as usize;
        &self.coeffs[idx]
    }
}

impl<'a, C: 'static + Coeff + Send + Sync> AnonSeriesSlice<'a, C> {
    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{anon_series::AnonSeries, AsSlice};
    ///
    /// let s = AnonSeries::new(-1, vec![1, 2, 3]);
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
        if pow < self.min_pow() {
            return Some(zero_ref());
        }
        if pow >= self.cutoff_pow() {
            return None;
        }
        Some(self.coeff_in_range(pow))
    }
}

impl<C: Coeff> Index<isize> for AnonSeriesSlice<'_, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow) as usize]
    }
}

impl<C> MulInverse for AnonSeriesSlice<'_, C>
where
    C: Coeff + SubAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
{
    type Output = AnonSeries<C>;

    fn mul_inverse(self) -> Self::Output {
        let inv_min_pow = -self.min_pow;
        if self.coeffs.is_empty() {
            return AnonSeries::new(inv_min_pow, vec![]);
        }
        let a: Vec<_> =
            self.coeffs.iter().map(|c| c / &self.coeffs[0]).collect();
        let mut b = Vec::with_capacity(a.len());
        b.push(C::one());
        for n in 1..a.len() {
            let mut b_n = C::zero();
            for i in 0..n {
                b_n -= &a[n - i] * &b[i];
            }
            b.push(b_n);
        }
        let inv_coeffs: Vec<_> =
            b.iter().map(|b| b / &self.coeffs[0]).collect();
        AnonSeries::new(inv_min_pow, inv_coeffs)
    }
}

impl<C: Coeff> Neg for AnonSeriesSlice<'_, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = AnonSeries<C>;

    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        AnonSeries::new(self.min_pow, neg_coeff)
    }
}

impl<C: Coeff + Clone, Rhs> Add<Rhs> for AnonSeriesSlice<'_, C>
where
    AnonSeries<C>: AddAssign<Rhs>,
{
    type Output = AnonSeries<C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = AnonSeries::from(self);
        res += other;
        res
    }
}

impl<C: Coeff, T> Sub<T> for AnonSeriesSlice<'_, C>
where
    C: Clone,
    AnonSeries<C>: SubAssign<T>,
{
    type Output = AnonSeries<C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = AnonSeries::from(self);
        res -= other;
        res
    }
}

impl<'a, C: Coeff> Mul for AnonSeriesSlice<'a, C>
where
    C: Clone,
    AnonSeries<C>: Mul<AnonSeriesSlice<'a, C>, Output = AnonSeries<C>>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: AnonSeriesSlice<'a, C>) -> Self::Output {
        AnonSeries::from(self) * other
    }
}

impl<C: Coeff> Mul<AnonSeries<C>> for AnonSeriesSlice<'_, C>
where
    C: Clone,
    AnonSeries<C>: MulAssign<AnonSeries<C>>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: AnonSeries<C>) -> Self::Output {
        AnonSeries::from(self) * other
    }
}

impl<'b, C: Coeff> Mul<&'b AnonSeries<C>> for AnonSeriesSlice<'_, C>
where
    C: Clone,
    for<'c> AnonSeries<C>: Mul<AnonSeriesSlice<'c, C>, Output = AnonSeries<C>>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: &'b AnonSeries<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<C: Coeff> Mul<C> for AnonSeriesSlice<'_, C>
where
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &other).collect();
        AnonSeries::new(self.min_pow(), coeffs)
    }
}

impl<'b, C: Coeff> Mul<&'b C> for AnonSeriesSlice<'_, C>
where
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = AnonSeries<C>;

    fn mul(self, other: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * other).collect();
        AnonSeries::new(self.min_pow(), coeffs)
    }
}

impl<C: Coeff, T> Div<T> for AnonSeriesSlice<'_, C>
where
    C: Clone,
    AnonSeries<C>: DivAssign<T>,
{
    type Output = AnonSeries<C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = AnonSeries::from(self);
        res /= other;
        res
    }
}

impl<C: Coeff> Exp for AnonSeriesSlice<'_, C>
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

    fn exp(self) -> Self::Output {
        AnonSeries::new(0, self.exp_coeff())
    }
}

impl<C: Coeff> Ln for AnonSeriesSlice<'_, C>
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
        assert_eq!(self.min_pow(), 0);
        assert!(!self.coeffs.is_empty());
        let c_k0 = &self.coeffs[0];
        // self.coeffs[0] = C::one();
        // for i in 1..self.coeffs.len() {
        //     self.coeffs[i] /= &c_k0;
        // }
        let a = &self.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = c_k0.ln();
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone() / c_k0);
            for i in 1..n {
                let num_factor = C::from(i as i32) / C::from(n as i32);
                let tmp = num_factor * (&a[n - i] * &b[i]) / c_k0;
                b[n] -= tmp;
            }
        }
        AnonSeries::new(0, b)
    }
}

impl<C: Coeff, T> Pow<T> for AnonSeriesSlice<'_, C>
where
    for<'b> AnonSeriesSlice<'b, C>: Ln<Output = AnonSeries<C>>,
    AnonSeries<C>: Mul<T>,
    <AnonSeries<C> as Mul<T>>::Output: Exp,
{
    type Output = <<AnonSeries<C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}

impl<'a, C: Coeff, Var> From<SeriesSlice<'a, Var, C>>
    for AnonSeriesSlice<'a, C>
{
    fn from(source: SeriesSlice<'a, Var, C>) -> Self {
        source.series
    }
}
