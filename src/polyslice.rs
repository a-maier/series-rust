use crate::poly::MulHelper;
use crate::traits::{AsSlice, KaratsubaMul};
use crate::util::{trim_slice_end, trim_slice_start};
use crate::{Coeff, Iter, Polynomial};

use std::fmt;
use std::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

/// View into a Laurent polynomial
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct PolynomialSlice<'a, Var, C: Coeff> {
    pub(crate) var: &'a Var,
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: &'a [C],
    pub(crate) zero: &'a C,
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, Var, C: Coeff> std::marker::Copy for PolynomialSlice<'a, Var, C> {}

impl<'a, Var, C: Coeff> std::clone::Clone for PolynomialSlice<'a, Var, C> {
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
        let mut res = PolynomialSlice {
            var,
            min_pow: Some(min_pow),
            coeffs,
            zero,
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
            if let Some(min_pow) = self.min_pow.as_mut() {
                *min_pow += removed as isize
            }
        }
    }

    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).min_pow(), Some(-1));
    /// assert_eq!(p.as_slice(0..).min_pow(), Some(0));
    /// let p = series::Polynomial::new("x", -1, vec![0]);
    /// assert_eq!(p.as_slice(..).min_pow(), None);
    /// ```
    pub fn min_pow(&self) -> Option<isize> {
        self.min_pow
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
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).max_pow(), Some(1));
    /// assert_eq!(p.as_slice(..1).max_pow(), Some(0));
    /// let p = series::Polynomial::new("x", -1, vec![0]);
    /// assert_eq!(p.max_pow(), None);
    /// ```
    pub fn max_pow(&self) -> Option<isize> {
        self.min_pow.map(|c| c - 1 + self.coeffs.len() as isize)
    }

    /// Get the difference between the highest and the lowest power of
    /// the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).len(), 3);
    /// assert_eq!(p.as_slice(0..2).len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Check if the polynomial is zero
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert!(!p.as_slice(..).is_empty());
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(0));
    /// assert!(p.as_slice(..).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn coeff(&self, pow: isize) -> &C {
        if let Some(min_pow) = self.min_pow() {
            let idx = pow - min_pow;
            if idx < 0 || idx > self.len() as isize {
                self.zero
            } else {
                &self.coeffs[idx as usize]
            }
        } else {
            self.zero
        }
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// let slice = p.as_slice(..);
    /// let mut iter = slice.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.iter())
    }

    /// Split a polynomial slice into two at the given power
    /// of the polynomial variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// let (lower, upper) = p.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), Some(-1));
    /// assert_eq!(upper.min_pow(), Some(0));
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos - self.min_pow().unwrap()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = PolynomialSlice {
            var: self.var,
            min_pow: self.min_pow(),
            coeffs: lower,
            zero: self.zero,
        };
        let upper = PolynomialSlice {
            var: self.var,
            min_pow: Some(pos),
            coeffs: upper,
            zero: self.zero,
        };
        (lower, upper)
    }
}

impl<'a, Var, C: Coeff> Index<isize> for PolynomialSlice<'a, Var, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow.unwrap()) as usize]
    }
}

impl<'a, Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for PolynomialSlice<'a, Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(min_pow) = self.min_pow() {
            for (i, coeff) in self.coeffs.iter().enumerate() {
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

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs>
    for PolynomialSlice<'a, Var, C>
where
    Polynomial<Var, C>: AddAssign<Rhs>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = Polynomial::from(self);
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
        let mut res = Polynomial::from(self);
        res -= other;
        res
    }
}

const MIN_KARATSUBA_SIZE: usize = 8;

impl<'a, 'b, Var, C: Coeff> Mul<PolynomialSlice<'b, Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<Var, C>:
        AddAssign<&'c Polynomial<Var, C>> + SubAssign<&'c Polynomial<Var, C>>,
    Polynomial<Var, C>:
        AddAssign<Polynomial<Var, C>> + SubAssign<Polynomial<Var, C>>,
    for<'c> PolynomialSlice<'c, Var, C>: Add<Output = Polynomial<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: PolynomialSlice<'b, Var, C>) -> Self::Output {
        assert_eq!(self.var, other.var);
        self.karatsuba_mul(other, MIN_KARATSUBA_SIZE)
    }
}

impl<'a, Var, C: Coeff> Mul<Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    for<'b> PolynomialSlice<'a, Var, C>:
        Mul<PolynomialSlice<'b, Var, C>, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: Polynomial<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b Polynomial<Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    PolynomialSlice<'a, Var, C>:
        Mul<PolynomialSlice<'b, Var, C>, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: &'b Polynomial<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var: Clone, C: Coeff> Mul<C> for PolynomialSlice<'a, Var, C>
where
    for<'b> &'b C: Mul<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Mul<&'b C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, Var: Clone, C: Coeff> Div<C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn div(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / &scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> Div<&'b C> for PolynomialSlice<'a, Var, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn div(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / scalar).collect();
        Polynomial::new(self.var.clone(), self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<PolynomialSlice<'b, Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<Var, C>:
        AddAssign<&'c Polynomial<Var, C>> + SubAssign<&'c Polynomial<Var, C>>,
    Polynomial<Var, C>:
        AddAssign<Polynomial<Var, C>> + SubAssign<Polynomial<Var, C>>,
    for<'c> PolynomialSlice<'c, Var, C>: Add<Output = Polynomial<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: PolynomialSlice<'b, Var, C>,
        min_size: usize,
    ) -> Self::Output {
        let mut result = Polynomial {
            var: self.var.clone(),
            min_pow: None,
            coeffs: vec![],
            zero: C::from(0),
        };
        result.add_prod(self, rhs, min_size);
        result
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<&'b Polynomial<Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<Var, C>:
        AddAssign<&'c Polynomial<Var, C>> + SubAssign<&'c Polynomial<Var, C>>,
    Polynomial<Var, C>:
        AddAssign<Polynomial<Var, C>> + SubAssign<Polynomial<Var, C>>,
    for<'c> PolynomialSlice<'c, Var, C>: Add<Output = Polynomial<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: &'b Polynomial<Var, C>,
        min_size: usize,
    ) -> Self::Output {
        self.karatsuba_mul(rhs.as_slice(..), min_size)
    }
}
