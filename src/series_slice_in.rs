use crate::ops::{Exp, Ln, Pow};
use crate::traits::{AsSlice, ExpCoeff, MulInverse};
use crate::util::trim_slice_start;
use crate::{Coeff, Iter, PolynomialSliceIn, SeriesIn};

use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// View into a Laurent series
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct SeriesSliceIn<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) min_pow: isize,
    pub(crate) coeffs: &'a [C],
    pub(crate) zero: &'a C,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, Var, C: Coeff> std::marker::Copy for SeriesSliceIn<'a, Var, C> {}

impl<'a, Var, C: Coeff> std::clone::Clone for SeriesSliceIn<'a, Var, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, Var, C: Coeff> SeriesSliceIn<'a, Var, C> {
    pub(super) fn new(
        var: &'a Var,
        min_pow: isize,
        coeffs: &'a [C],
        zero: &'a C,
    ) -> Self {
        let mut res = SeriesSliceIn {
            var,
            min_pow,
            coeffs,
            zero,
        };
        res.trim();
        res
    }

    fn trim(&mut self) {
        let (coeffs, removed) = trim_slice_start(self.coeffs, self.zero);
        self.coeffs = coeffs;
        self.min_pow += removed as isize;
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec![1,2,3]);
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
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.as_slice(..).cutoff_pow(), 2);
    /// assert_eq!(s.as_slice(..1).cutoff_pow(), 1);
    /// ```
    pub fn cutoff_pow(&self) -> isize {
        self.min_pow + (self.coeffs.len() as isize)
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
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
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
            return Some(self.zero); // TODO this is a bad hack
        }
        if pow >= self.cutoff_pow() {
            return None;
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(&self.coeffs[idx])
    }

    /// Iterator over the series powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..);
    /// let mut iter = slice.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        (self.min_pow..).zip(self.coeffs.iter())
    }

    /// Split a series slice into two at the given power
    /// of the expansion variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// let (lower, upper) = s.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), -1);
    /// assert_eq!(upper.min_pow(), 0);
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos - self.min_pow()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = SeriesSliceIn {
            var: self.var,
            min_pow: self.min_pow(),
            coeffs: lower,
            zero: self.zero,
        };
        let upper = SeriesSliceIn {
            var: self.var,
            min_pow: pos,
            coeffs: upper,
            zero: self.zero,
        };
        (lower, upper)
    }

    /// View as polynomial slice
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// let slice = s.as_slice(..).as_poly();
    /// let p = series::PolynomialIn::from(s.clone());
    /// assert_eq!(slice, p.as_slice(..));
    /// ```
    pub fn as_poly(&self) -> PolynomialSliceIn<'a, Var, C> {
        PolynomialSliceIn::new(self.var, self.min_pow, self.coeffs, self.zero)
    }

    /// Get the expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::SeriesIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.as_slice(..).var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        self.var
    }
}

impl<'a, Var, C: Coeff> Index<isize> for SeriesSliceIn<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow) as usize]
    }
}

impl<'a, Var, C> MulInverse for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    C: Coeff + SubAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    fn mul_inverse(self) -> Self::Output {
        let inv_min_pow = -self.min_pow;
        if self.coeffs.is_empty() {
            return SeriesIn::new(self.var.clone(), inv_min_pow, vec![]);
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
        SeriesIn::new(self.var.clone(), inv_min_pow, inv_coeffs)
    }
}

impl<'a, Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for SeriesSliceIn<'a, Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if !self.coeffs.is_empty() {
            self.as_poly().fmt(f)?;
            write!(f, " + ")?;
        }
        write!(f, "O({}^{})", self.var, self.cutoff_pow())
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for SeriesSliceIn<'a, Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        SeriesIn::new(self.var.clone(), self.min_pow, neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs>
    for SeriesSliceIn<'a, Var, C>
where
    SeriesIn<Var, C>: AddAssign<Rhs>,
{
    type Output = SeriesIn<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = SeriesIn::from(self);
        res += other;
        res
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for SeriesSliceIn<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    SeriesIn<Var, C>: SubAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = SeriesIn::from(self);
        res -= other;
        res
    }
}

impl<'a, Var, C: Coeff> Mul for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    SeriesIn<Var, C>: Mul<SeriesSliceIn<'a, Var, C>, Output = SeriesIn<Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: SeriesSliceIn<'a, Var, C>) -> Self::Output {
        SeriesIn::from(self) * other
    }
}

impl<'a, Var, C: Coeff> Mul<SeriesIn<Var, C>> for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    SeriesIn<Var, C>: MulAssign<SeriesIn<Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: SeriesIn<Var, C>) -> Self::Output {
        SeriesIn::from(self) * other
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b SeriesIn<Var, C>>
    for SeriesSliceIn<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    for<'c> SeriesIn<Var, C>:
        Mul<SeriesSliceIn<'c, Var, C>, Output = SeriesIn<Var, C>>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: &'b SeriesIn<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var, C: Coeff> Mul<C> for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &other).collect();
        SeriesIn::new(self.var.clone(), self.min_pow(), coeffs)
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b C> for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = SeriesIn<Var, C>;

    fn mul(self, other: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * other).collect();
        SeriesIn::new(self.var.clone(), self.min_pow(), coeffs)
    }
}

impl<'a, Var, C: Coeff, T> Div<T> for SeriesSliceIn<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    SeriesIn<Var, C>: DivAssign<T>,
{
    type Output = SeriesIn<Var, C>;

    fn div(self, other: T) -> Self::Output {
        let mut res = SeriesIn::from(self);
        res /= other;
        res
    }
}

impl<'a, Var, C: Coeff> Exp for SeriesSliceIn<'a, Var, C>
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

    fn exp(self) -> Self::Output {
        SeriesIn::new(self.var.clone(), 0, self.exp_coeff())
    }
}

impl<'a, Var, C: Coeff> Ln for SeriesSliceIn<'a, Var, C>
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
    type Output = SeriesIn<Var, C>;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(self) -> Self::Output {
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = &self.coeffs[0];
        // self.coeffs[0] = C::one();
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
        SeriesIn::new(self.var.clone(), 0, b)
    }
}

impl<'a, Var, C: Coeff, T> Pow<T> for SeriesSliceIn<'a, Var, C>
where
    for<'b> SeriesSliceIn<'b, Var, C>: Ln<Output = SeriesIn<Var, C>>,
    SeriesIn<Var, C>: Mul<T>,
    <SeriesIn<Var, C> as Mul<T>>::Output: Exp,
{
    type Output = <<SeriesIn<Var, C> as Mul<T>>::Output as Exp>::Output;

    fn pow(self, exponent: T) -> Self::Output {
        (self.ln() * exponent).exp()
    }
}
