use crate::traits::{AsSlice, KaratsubaMul};
use crate::{Polynomial, PolynomialSlice};
use crate::{Coeff, Iter, PolynomialIn};

use std::fmt;
use std::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

/// View into a Laurent polynomial
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct PolynomialSliceIn<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) poly: PolynomialSlice<'a, C>,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, Var, C: Coeff> std::marker::Copy for PolynomialSliceIn<'a, Var, C> {}

impl<'a, Var, C: Coeff> std::clone::Clone for PolynomialSliceIn<'a, Var, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, Var, C: Coeff> PolynomialSliceIn<'a, Var, C> {
    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).min_pow(), Some(-1));
    /// assert_eq!(p.as_slice(0..).min_pow(), Some(0));
    /// let p = series::PolynomialIn::new("x", -1, vec![0]);
    /// assert_eq!(p.as_slice(..).min_pow(), None);
    /// ```
    pub fn min_pow(&self) -> Option<isize> {
        self.poly.min_pow()
    }

    /// Get the highest power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).max_pow(), Some(1));
    /// assert_eq!(p.as_slice(..1).max_pow(), Some(0));
    /// let p = series::PolynomialIn::new("x", -1, vec![0]);
    /// assert_eq!(p.max_pow(), None);
    /// ```
    pub fn max_pow(&self) -> Option<isize> {
        self.poly.max_pow()
    }

    /// Get the difference between the highest and the lowest power of
    /// the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).len(), 3);
    /// assert_eq!(p.as_slice(0..2).len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.poly.len()
    }

    /// Check if the polynomial is zero
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert!(!p.as_slice(..).is_empty());
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(0));
    /// assert!(p.as_slice(..).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.poly.is_empty()
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// let slice = p.as_slice(..);
    /// let mut iter = slice.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        self.poly.iter()
    }

    /// Split a polynomial slice into two at the given power
    /// of the polynomial variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// let (lower, upper) = p.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), Some(-1));
    /// assert_eq!(upper.min_pow(), Some(0));
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let (lower, upper) = self.poly.split_at(pos);
        let lower = lower.in_var(self.var());
        let upper = upper.in_var(self.var());
        (lower, upper)
    }

    /// Get the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let s = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.as_slice(..).var(), &"x");
    /// ```
    pub fn var(&self) -> &'a Var {
        self.var
    }

    pub fn try_coeff(&self, pow: isize) -> Option<&'a C> {
        self.poly.try_coeff(pow)
    }
}

impl<'a, Var, C: Coeff> Index<isize> for PolynomialSliceIn<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.poly[index]
    }
}

impl<'a, Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for PolynomialSliceIn<'a, Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(min_pow) = self.min_pow() {
            for (i, coeff) in self.poly.coeffs.iter().enumerate() {
                if *coeff == C::zero() {
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

impl<'a, Var: Clone, C: Coeff> Neg for PolynomialSliceIn<'a, Var, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn neg(self) -> Self::Output {
        self.poly.neg().in_var(self.var.clone())
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs>
    for PolynomialSliceIn<'a, Var, C>
where
    PolynomialIn<Var, C>: AddAssign<Rhs>,
{
    type Output = PolynomialIn<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = PolynomialIn::from(self);
        res += other;
        res
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for PolynomialSliceIn<'a, Var, C>
where
    C: Clone,
    Var: Clone,
    PolynomialIn<Var, C>: SubAssign<T>,
{
    type Output = PolynomialIn<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = PolynomialIn::from(self);
        res -= other;
        res
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<PolynomialSliceIn<'b, Var, C>>
    for PolynomialSliceIn<'a, Var, C>
where
    Var: Clone + fmt::Debug + PartialEq,
    PolynomialSlice<'a, C>: Mul<PolynomialSlice<'b, C>, Output = Polynomial<C>>
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, other: PolynomialSliceIn<'b, Var, C>) -> Self::Output {
        assert_eq!(self.var(), other.var());
        self.poly.mul(other.poly).in_var(self.var.clone())
    }
}

impl<'a, Var, C: Coeff> Mul<PolynomialIn<Var, C>>
    for PolynomialSliceIn<'a, Var, C>
where
    for<'b> PolynomialSliceIn<'a, Var, C>:
        Mul<PolynomialSliceIn<'b, Var, C>, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, other: PolynomialIn<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b PolynomialIn<Var, C>>
    for PolynomialSliceIn<'a, Var, C>
where
    PolynomialSliceIn<'a, Var, C>:
        Mul<PolynomialSliceIn<'b, Var, C>, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, other: &'b PolynomialIn<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var: Clone, C: Coeff> Mul<C> for PolynomialSliceIn<'a, Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, scalar: C) -> Self::Output {
        self.poly.mul(scalar).in_var(self.var().clone())
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Mul<&'b C> for PolynomialSliceIn<'a, Var, C>
where
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, scalar: &'b C) -> Self::Output {
        self.poly.mul(scalar).in_var(self.var().clone())
    }
}

impl<'a, Var: Clone, C: Coeff> Div<C> for PolynomialSliceIn<'a, Var, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn div(self, scalar: C) -> Self::Output {
        self.poly.div(scalar).in_var(self.var().clone())
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Div<&'b C> for PolynomialSliceIn<'a, Var, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn div(self, scalar: &'b C) -> Self::Output {
        self.poly.div(scalar).in_var(self.var().clone())
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<PolynomialSliceIn<'b, Var, C>>
    for PolynomialSliceIn<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    PolynomialSlice<'a, C>: KaratsubaMul<PolynomialSlice<'b, C>, Output = Polynomial<C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: PolynomialSliceIn<'b, Var, C>,
        min_size: usize,
    ) -> Self::Output {
        assert_eq!(self.var, rhs.var);
        self.poly.karatsuba_mul(rhs.poly, min_size).in_var(self.var().clone())
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<&'b PolynomialIn<Var, C>>
    for PolynomialSliceIn<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    PolynomialSlice<'a, C>: KaratsubaMul<&'b Polynomial<C>, Output = Polynomial<C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: &'b PolynomialIn<Var, C>,
        min_size: usize,
    ) -> Self::Output {
        self.poly.karatsuba_mul(&rhs.poly, min_size).in_var(self.var.clone())
    }
}
