use crate::poly::MulHelper;
use crate::traits::{AsSlice, KaratsubaMul};
use crate::util::{trim_slice_end, trim_slice_start};
use crate::zero_ref::zero_ref;
use crate::{Coeff, Iter, Polynomial, PolynomialSliceIn};

use std::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

/// View into a Laurent polynomial
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub struct PolynomialSlice<'a, C: Coeff> {
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: &'a [C],
}

// needs manual implementation,
// #[derive(Copy)] can't deal with lifetimes in rust 1.36
impl<'a, C: Coeff> std::marker::Copy for PolynomialSlice<'a, C> {}

impl<'a, C: Coeff> std::clone::Clone for PolynomialSlice<'a, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C: Coeff> PolynomialSlice<'a, C> {
    pub(super) fn new(min_pow: isize, coeffs: &'a [C]) -> Self {
        let mut res = PolynomialSlice {
            min_pow: Some(min_pow),
            coeffs,
        };
        res.trim();
        res
    }

    fn trim(&mut self) {
        let (coeffs, _removed) = trim_slice_end(self.coeffs, &C::zero());
        self.coeffs = coeffs;
        if self.coeffs.is_empty() {
            self.min_pow = None;
        } else {
            let (coeffs, removed) = trim_slice_start(self.coeffs, &C::zero());
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).min_pow(), Some(-1));
    /// assert_eq!(p.as_slice(0..).min_pow(), Some(0));
    /// let p = series::Polynomial::new(-1, vec![0]);
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.as_slice(..).max_pow(), Some(1));
    /// assert_eq!(p.as_slice(..1).max_pow(), Some(0));
    /// let p = series::Polynomial::new(-1, vec![0]);
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert!(!p.as_slice(..).is_empty());
    ///
    /// let p = series::Polynomial::new(-1, vec!(0));
    /// assert!(p.as_slice(..).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// let (lower, upper) = p.as_slice(..).split_at(0);
    /// assert_eq!(lower.min_pow(), Some(-1));
    /// assert_eq!(upper.min_pow(), Some(0));
    /// ```
    pub fn split_at(&self, pos: isize) -> (Self, Self) {
        let upos = (pos - self.min_pow().unwrap()) as usize;
        let (lower, upper) = self.coeffs.split_at(upos);
        let lower = PolynomialSlice {
            min_pow: self.min_pow(),
            coeffs: lower,
        };
        let upper = PolynomialSlice {
            min_pow: Some(pos),
            coeffs: upper,
        };
        (lower, upper)
    }

    /// Turn into a slice with a named variable
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::AsSlice;
    ///
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// let p = p.as_slice(..).in_var(&"x");
    /// assert_eq!(p.var(), &"x");
    /// ```
    pub fn in_var<Var>(self, var: &'a Var) -> PolynomialSliceIn<Var, C> {
        PolynomialSliceIn {
            var,
            poly: self
        }
    }

    pub fn try_coeff(&self, pow: isize) -> Option<&'a C> {
        let min_pow = self.min_pow()?;
        self.coeffs.get((pow - min_pow) as usize)
    }
}

impl<'a, C: 'static + Coeff + Send + Sync> PolynomialSlice<'a, C> {
    pub fn coeff(&self, pow: isize) -> &'a C {
        if let Some(min_pow) = self.min_pow() {
            let idx = pow - min_pow;
            if idx < 0 || idx >= self.len() as isize {
                zero_ref()
            } else {
                &self.coeffs[idx as usize]
            }
        } else {
            zero_ref()
        }
    }
}

impl<'a, C: Coeff> Index<isize> for PolynomialSlice<'a, C> {
    type Output = C;

    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow.unwrap()) as usize]
    }
}

impl<'a, C: Coeff> Neg for PolynomialSlice<'a, C>
where
    for<'c> &'c C: Neg<Output = C>,
{
    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), neg_coeff)
    }
}

impl<'a, C: Coeff + Clone, Rhs> Add<Rhs> for PolynomialSlice<'a, C>
where
    Polynomial<C>: AddAssign<Rhs>,
{
    type Output = Polynomial<C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = Polynomial::from(self);
        res += other;
        res
    }
}

impl<'a, C: Coeff, T> Sub<T> for PolynomialSlice<'a, C>
where
    C: Clone,
    Polynomial<C>: SubAssign<T>,
{
    type Output = Polynomial<C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = Polynomial::from(self);
        res -= other;
        res
    }
}

const MIN_KARATSUBA_SIZE: usize = 8;

impl<'a, 'b, C: Coeff> Mul<PolynomialSlice<'b, C>> for PolynomialSlice<'a, C>
where
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<C>:
        AddAssign<&'c Polynomial<C>> + SubAssign<&'c Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>> + SubAssign<Polynomial<C>>,
    for<'c> PolynomialSlice<'c, C>: Add<Output = Polynomial<C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<C>;

    fn mul(self, other: PolynomialSlice<'b, C>) -> Self::Output {
        self.karatsuba_mul(other, MIN_KARATSUBA_SIZE)
    }
}

impl<'a, C: Coeff> Mul<Polynomial<C>> for PolynomialSlice<'a, C>
where
    for<'b> PolynomialSlice<'a, C>:
        Mul<PolynomialSlice<'b, C>, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn mul(self, other: Polynomial<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b Polynomial<C>> for PolynomialSlice<'a, C>
where
    PolynomialSlice<'a, C>: Mul<PolynomialSlice<'b, C>, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn mul(self, other: &'b Polynomial<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, C: Coeff> Mul<C> for PolynomialSlice<'a, C>
where
    for<'b> &'b C: Mul<Output = C>,
{
    type Output = Polynomial<C>;

    fn mul(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * &scalar).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b C> for PolynomialSlice<'a, C>
where
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<C>;

    fn mul(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c * scalar).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, C: Coeff> Div<C> for PolynomialSlice<'a, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = Polynomial<C>;

    fn div(self, scalar: C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / &scalar).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, 'b, C: Coeff> Div<&'b C> for PolynomialSlice<'a, C>
where
    for<'c> &'c C: Div<Output = C>,
{
    type Output = Polynomial<C>;

    fn div(self, scalar: &'b C) -> Self::Output {
        let coeffs = self.coeffs.iter().map(|c| c / scalar).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), coeffs)
    }
}

impl<'a, C: Coeff, Var> From<PolynomialSliceIn<'a, Var, C>>
    for PolynomialSlice<'a, C>
{
    fn from(source: PolynomialSliceIn<'a, Var, C>) -> Self {
        source.poly
    }
}

impl<'a, 'b, C: Coeff> KaratsubaMul<PolynomialSlice<'b, C>>
    for PolynomialSlice<'a, C>
where
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<C>:
        AddAssign<&'c Polynomial<C>> + SubAssign<&'c Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>> + SubAssign<Polynomial<C>>,
    for<'c> PolynomialSlice<'c, C>: Add<Output = Polynomial<C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<C>;

    fn karatsuba_mul(
        self,
        rhs: PolynomialSlice<'b, C>,
        min_size: usize,
    ) -> Self::Output {
        let mut result = Polynomial {
            min_pow: None,
            coeffs: vec![],
        };
        result.add_prod(self, rhs, min_size);
        result
    }
}

impl<'a, 'b, C: Coeff> KaratsubaMul<&'b Polynomial<C>>
    for PolynomialSlice<'a, C>
where
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<C>:
        AddAssign<&'c Polynomial<C>> + SubAssign<&'c Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>> + SubAssign<Polynomial<C>>,
    for<'c> PolynomialSlice<'c, C>: Add<Output = Polynomial<C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = Polynomial<C>;

    fn karatsuba_mul(
        self,
        rhs: &'b Polynomial<C>,
        min_size: usize,
    ) -> Self::Output {
        self.karatsuba_mul(rhs.as_slice(..), min_size)
    }
}
