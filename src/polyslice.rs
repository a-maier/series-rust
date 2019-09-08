use crate::{Coeff, Iter, Polynomial};
use crate::util::{trim_slice_start, trim_slice_end};
use crate::poly::MulHelper;
use crate::traits::AsSlice;

use std::fmt;
use std::ops::{
    Add, AddAssign, Div, Mul, Neg,
    Sub, SubAssign, Index
};

/// View into a Laurent polynomial
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct PolynomialSlice<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: &'a [C],
    pub(crate) zero: &'a C,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, Var, C: Coeff> std::marker::Copy for PolynomialSlice<'a, Var, C>
{}

impl<'a, Var, C: Coeff> std::clone::Clone for PolynomialSlice<'a, Var, C>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, Var, C: Coeff> PolynomialSlice<'a, Var, C> {
    pub(super) fn new(
        var: &'a Var,
        min_pow: isize,
        coeffs: &'a [C],
        zero: &'a C,
    ) -> Self {
        let mut res = PolynomialSlice{
            var,
            min_pow: Some(min_pow),
            coeffs,
            zero
        };
        res.trim();
        res
    }

    fn trim(&mut self) {
        let (coeffs, _removed) = trim_slice_end(self.coeffs, self.zero);
        self.coeffs = coeffs;
        if self.coeffs.is_empty() {
            self.min_pow = None;
        } else {
            let (coeffs, removed) = trim_slice_start(self.coeffs, self.zero);
            self.coeffs = coeffs;
            self.min_pow.as_mut().map(|m| *m += removed as isize);
        }
    }

    pub fn min_pow(&self) -> Option<isize> {
        self.min_pow
    }

    pub fn max_pow(&self) -> Option<isize> {
        self.min_pow.map(|c| c + self.coeffs.len() as isize)
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn coeff(&self, pow: isize) -> &C {
        if let Some(min_pow) = self.min_pow() {
            let idx = pow - min_pow;
            if idx < 0 || idx > self.len() as isize {
                &self.zero
            }
            else {
                &self.coeffs[idx as usize]
            }
        } else {
            &self.zero
        }
    }

    pub fn iter(&self) -> Iter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.iter())
    }

    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos - self.min_pow().unwrap()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = PolynomialSlice{
            var: self.var,
            min_pow: self.min_pow(),
            coeffs: lower,
            zero: self.zero
        };
        let upper = PolynomialSlice{
            var: self.var,
            min_pow: Some(pos),
            coeffs: upper,
            zero: self.zero
        };
        (lower, upper)
    }
}

impl<'a, Var, C: Coeff> Index<isize> for PolynomialSlice<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index-self.min_pow.unwrap()) as usize]
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> PolynomialSlice<'a, Var, C> {
    fn to_owned(&self) -> Polynomial<Var, C> {
        Polynomial::new(
            self.var.clone(),
            self.min_pow.unwrap_or(0),
            self.coeffs.to_vec()
        )
    }
}

impl<'a, Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for PolynomialSlice<'a, Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(min_pow) = self.min_pow() {
            for (i, coeff) in self.coeffs.iter().enumerate() {
                if *coeff == C::from(0) {
                    continue;
                }
                let cur_pow = min_pow + i as isize;
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
        }
        Ok(())
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for PolynomialSlice<'a, Var, C>
where
    Polynomial<Var, C>: AddAssign<Rhs>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.to_owned();
        res += other;
        res
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for PolynomialSlice<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    Polynomial<Var, C>: SubAssign<T>,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.to_owned();
        res -= other;
        res
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<PolynomialSlice<'b, Var, C>> for PolynomialSlice<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for <'c> C: AddAssign,
    for <'c> Polynomial<Var, C>: AddAssign<&'c Polynomial<Var, C>> + SubAssign<&'c Polynomial<Var, C>>,
    Polynomial<Var, C>: AddAssign<Polynomial<Var, C>> + SubAssign<Polynomial<Var, C>>,
    for<'c> PolynomialSlice<'c, Var, C>: Add<Output = Polynomial<Var, C>>,
    for <'c> &'c C: Mul<Output=C>
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: PolynomialSlice<'b, Var, C>) -> Self::Output {
        assert_eq!(self.var, other.var);
        let mut res = Polynomial::new(self.var.clone(), 0, vec![]);
        res.add_prod(self, other);
        res
    }
}

impl<'a, Var, C: Coeff> Mul<Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    for<'b> PolynomialSlice<'a, Var, C>: Mul<PolynomialSlice<'b, Var, C>, Output=Polynomial<Var, C>>
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: Polynomial<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    PolynomialSlice<'a, Var, C>: Mul<PolynomialSlice<'b, Var, C>, Output=Polynomial<Var, C>>
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: &'b Polynomial<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var: Clone, C: Coeff> Mul<C> for PolynomialSlice<'a, Var, C>
where
    for<'b> &'b C: Mul<Output=C>
{
    type Output = Polynomial<Var, C>;

    fn mul(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Mul<&'b C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Mul<Output=C>
{
    type Output = Polynomial<Var, C>;

    fn mul(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, Var: Clone, C: Coeff> Div<C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Div<Output=C>
{
    type Output = Polynomial<Var, C>;

    fn div(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / &scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Div<&'b C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Div<Output=C>
{
    type Output = Polynomial<Var, C>;

    fn div(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}
