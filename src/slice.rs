use crate::{Series, Coeff, Iter, PolynomialSlice};
use crate::traits::{MulInverse, ExpCoeff, AsSlice};
use crate::ops::{Exp, Ln, Pow};
use crate::util::trim_slice_start;

use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub,
    SubAssign, Index
};

/// View into a Laurent series
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct SeriesSlice<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) min_pow: isize,
    pub(crate) coeffs: &'a [C],
    pub(crate) zero: &'a C,
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
    pub(super) fn new(
        var: &'a Var,
        min_pow: isize,
        coeffs: &'a [C],
        zero: &'a C,
    ) -> Self {
        let mut res = SeriesSlice{var, min_pow, coeffs, zero};
        res.trim();
        res
    }

    fn trim(&mut self) {
        let (coeffs, removed) = trim_slice_start(self.coeffs, self.zero);
        self.coeffs = coeffs;
        self.min_pow += removed as isize;
    }

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

    pub fn as_poly(&self) -> PolynomialSlice<'a, Var, C> {
        PolynomialSlice::new(self.var, self.min_pow, self.coeffs, self.zero)
    }
}

impl<'a, Var, C: Coeff> Index<isize> for SeriesSlice<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index-self.min_pow) as usize]
    }
}

// TODO: impl std::borrow::ToOwned
// but this conflicts with impl for T: Clone in std
impl<'a, Var: Clone, C: Coeff + Clone> SeriesSlice<'a, Var, C> {
    pub fn to_owned(&self) ->  Series<Var, C> {
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
        if !self.coeffs.is_empty() {
            self.as_poly().fmt(f)?;
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

impl<'a, Var, C: Coeff> Mul for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    C: Clone,
    Series<Var, C>: Mul<SeriesSlice<'a, Var, C>, Output=Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: SeriesSlice<'a, Var, C>) -> Self::Output {
        self.to_owned() * other
    }
}

impl<'a, Var, C: Coeff> Mul<Series<Var, C>> for SeriesSlice<'a, Var, C>
where
    Series<Var, C>: Mul<SeriesSlice<'a, Var, C>, Output=Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: Series<Var, C>) -> Self::Output {
        other * self
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b Series<Var, C>> for SeriesSlice<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    for<'c> Series<Var, C>: Mul<SeriesSlice<'c, Var, C>, Output=Series<Var, C>>,
{
    type Output = Series<Var, C>;

    fn mul(self, other: &'b Series<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var, C: Coeff> Mul<C> for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    for<'c>  &'c C: Mul<Output=C>
{
    type Output = Series<Var, C>;

    fn mul(self, other: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &other).collect();
        Series::new(
            self.var.clone(),
            self.min_pow(),
            coeffs
        )
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b C> for SeriesSlice<'a, Var, C>
where
    Var: Clone,
    for<'c>  &'c C: Mul<Output=C>
{
    type Output = Series<Var, C>;

    fn mul(self, other: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * other).collect();
        Series::new(
            self.var.clone(),
            self.min_pow(),
            coeffs
        )
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
