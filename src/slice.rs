pub use crate::{Series, Coeff, Iter, MulInverse};
use crate::{AddAssignHelper, ExpCoeff};
use crate::ops::{Exp, Ln, Pow};

use std::cmp::min;
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub,
    SubAssign, Index
};

/// View into a Laurent series
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct SeriesSlice<'a, Var, C: Coeff> {
    pub(super) var: &'a Var,
    pub(super) min_pow: isize,
    pub(super) coeffs: &'a [C],
    pub(super) zero: &'a C,
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