use crate::traits::{AsSlice, KaratsubaMul};
use crate::util::{trim_end, trim_start};
use crate::{Coeff, IntoIter, Iter, PolynomialIn, PolynomialInParts};
use crate::{PolynomialSlice, Series};

use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};
use std::{convert, iter};

use log::trace;
use num_traits::{One, Zero};

/// Laurent polynomial in a single variable
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct Polynomial<C: Coeff> {
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: Vec<C>,
    pub(crate) zero: C, // TODO: evil hack, learn how to do this better
}

impl<C: Coeff> Polynomial<C> {
    /// Create a new Laurent polynomial
    ///
    /// # Example
    ///
    /// This creates a Laurent polynomial starting at power -1 with
    /// coefficients 1, 2, 3. In other words, the polynomial x^-1 + 2
    /// + 3*x.
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// ```
    pub fn new(min_pow: isize, coeffs: Vec<C>) -> Polynomial<C> {
        let mut res = Polynomial {
            min_pow: Some(min_pow),
            coeffs,
            zero: C::zero(),
        };
        res.trim();
        res
    }

    /// Turn into a polynomial in a named variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Polynomial::new(-1, vec!(1,2,3)).in_var("x");
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn in_var<Var>(self, var: Var) -> PolynomialIn<Var, C> {
        let Self {
            min_pow,
            coeffs,
            zero,
        } = self;
        PolynomialIn {
            var,
            min_pow,
            coeffs,
            zero,
        }
    }

    /// Get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.coeff(-5), &0);
    /// assert_eq!(p.coeff(-2), &0);
    /// assert_eq!(p.coeff(-1), &1);
    /// assert_eq!(p.coeff(0), &2);
    /// assert_eq!(p.coeff(1), &3);
    /// assert_eq!(p.coeff(2), &0);
    /// assert_eq!(p.coeff(5), &0);
    /// ```
    pub fn coeff(&self, pow: isize) -> &C {
        self.as_slice(..).coeff(pow)
    }

    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.min_pow(), Some(-1));
    /// let p = series::Polynomial::new(-1, vec![0]);
    /// assert_eq!(p.min_pow(), None);
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
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.max_pow(), Some(1));
    /// let p = series::Polynomial::new(-1, vec![0]);
    /// assert_eq!(p.max_pow(), None);
    /// ```
    pub fn max_pow(&self) -> Option<isize> {
        self.min_pow.map(|c| c - 1 + self.len() as isize)
    }

    /// Get the difference between the highest and the lowest power of
    /// the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p.len(), 3);
    /// ```
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Check if the polynomial is zero
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert!(!p.is_empty());
    ///
    /// let p = series::Polynomial::new(-1, vec!(0));
    /// assert!(p.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// let mut iter = p.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.iter())
    }

    /// Turn a polynomial into a series with the given cutoff
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// let s = series::Series::with_cutoff(-1..5, vec!(1,2,3));
    /// assert_eq!(p.cutoff_at(5), s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the cutoff power is lower than the starting power
    ///
    /// ```
    pub fn cutoff_at(self, cutoff_pow: isize) -> Series<C> {
        Series::with_cutoff(
            self.min_pow.unwrap_or(cutoff_pow)..cutoff_pow,
            self.coeffs,
        )
    }

    fn as_empty_slice(&self) -> PolynomialSlice<'_, C> {
        PolynomialSlice::new(0, &self.coeffs[self.len()..], &self.zero)
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, Range<isize>> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = (r.start - min_pow) as usize;
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(r.start, &self.coeffs[start..end], &self.zero)
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeInclusive<isize>> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let (start, end) = r.into_inner();
            let ustart = (start - min_pow) as usize;
            let end = (end - min_pow) as usize;
            PolynomialSlice::new(start, &self.coeffs[ustart..=end], &self.zero)
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(min_pow, &self.coeffs[..=end], &self.zero)
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = r.start - min_pow;
            PolynomialSlice::new(
                r.start,
                &self.coeffs[(start as usize)..],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(min_pow, &self.coeffs[..end], &self.zero)
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, C: 'a + Coeff> AsSlice<'a, RangeFull> for Polynomial<C> {
    type Output = PolynomialSlice<'a, C>;

    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            PolynomialSlice::new(min_pow, &self.coeffs[r], &self.zero)
        } else {
            self.as_empty_slice()
        }
    }
}

impl<C: Coeff> convert::From<Series<C>> for Polynomial<C> {
    fn from(s: Series<C>) -> Self {
        Polynomial::new(s.min_pow, s.coeffs)
    }
}

impl<C: Coeff> Index<isize> for Polynomial<C> {
    type Output = C;

    /// Get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// # Panics
    ///
    /// Panics if the index is smaller than the leading power or
    /// bigger than the highest power
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// assert_eq!(p[-1], 1);
    /// assert_eq!(p[0], 2);
    /// assert_eq!(p[1], 3);
    /// assert!(std::panic::catch_unwind(|| p[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| p[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index - self.min_pow.unwrap()) as usize]
    }
}

impl<C: Coeff> std::iter::IntoIterator for Polynomial<C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-1, vec!(1,2,3));
    /// let mut iter = p.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs)
    }
}

impl<C: Coeff> Polynomial<C> {
    fn trim(&mut self) {
        trim_end(&mut self.coeffs, &self.zero);
        if self.coeffs.is_empty() {
            self.min_pow = None;
        } else {
            let min_pow_shift = trim_start(&mut self.coeffs, &self.zero);
            if let Some(min_pow) = self.min_pow.as_mut() {
                *min_pow += min_pow_shift as isize
            }
        }
    }

    fn extend_min(&mut self, extend: usize) {
        debug_assert!(self.min_pow.is_some());
        let to_insert = iter::repeat_with(C::zero).take(extend);
        self.coeffs.splice(0..0, to_insert);
        if let Some(min_pow) = self.min_pow.as_mut() {
            *min_pow -= extend as isize
        }
    }

    fn extend_max(&mut self, extend: usize) {
        let to_insert = iter::repeat_with(C::zero).take(extend);
        self.coeffs.extend(to_insert);
    }
}

impl<C: Coeff + Neg<Output = C>> Neg for Polynomial<C> {
    type Output = Polynomial<C>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let minus_p = series::Polynomial::new(-3, vec!(-1.,0.,3.));
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        Polynomial::new(self.min_pow.unwrap_or(0), neg_coeff)
    }
}

impl<'a, C: Coeff> Neg for &'a Polynomial<C>
where
    PolynomialSlice<'a, C>: Neg,
{
    type Output = <PolynomialSlice<'a, C> as Neg>::Output;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let minus_p = series::Polynomial::new(-3, vec!(-1.,0.,3.));
    /// assert_eq!(-&p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, C: Coeff + Clone> AddAssign<&'a Polynomial<C>> for Polynomial<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set p = p + q for two Laurent polynomials p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let q = Polynomial::new(-1, vec!(3., 4., 5.));
    /// let res = Polynomial::new(-3, vec!(1.,0.,0.,4.,5.));
    /// p += &q;
    /// assert_eq!(res, p);
    /// ```
    ///
    fn add_assign(&mut self, other: &'a Polynomial<C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, C: Coeff + Clone> AddAssign<PolynomialSlice<'a, C>> for Polynomial<C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSlice<'a, C>) {
        if other.min_pow().is_none() {
            return;
        }
        if self.min_pow().is_none() {
            self.coeffs = other.coeffs.to_owned();
            self.min_pow = other.min_pow();
            return;
        }
        let min_pow = self.min_pow().unwrap();
        let other_min_pow = other.min_pow().unwrap();
        if other_min_pow < min_pow {
            self.extend_min((min_pow - other_min_pow) as usize);
        }
        let max_pow = self.max_pow().unwrap();
        let other_max_pow = other.max_pow().unwrap();
        if other_max_pow > max_pow {
            self.extend_max((other_max_pow - max_pow) as usize);
        }
        debug_assert!(self.min_pow().unwrap() <= other_min_pow);
        debug_assert!(self.max_pow().unwrap() >= other_max_pow);
        let min_pow = self.min_pow().unwrap();
        for (pow, coeff) in other.iter() {
            self.coeffs[(pow - min_pow) as usize] += coeff;
        }
        self.trim();
    }
}

impl<C: Coeff> AddAssign<Polynomial<C>> for Polynomial<C>
where
    for<'c> C: AddAssign<&'c C>,
    C: Clone + AddAssign,
{
    /// Set p = p + q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let q = Polynomial::new(-1, vec!(3., 4., 5.));
    /// let res = Polynomial::new(-3, vec!(1.,0.,0.,4.,5.));
    /// p += q;
    /// assert_eq!(res, p);
    /// ```
    fn add_assign(&mut self, other: Polynomial<C>) {
        //TODO: code duplication with AddAssign<PolynomialSlice>
        if other.min_pow().is_none() {
            return;
        }
        if self.min_pow().is_none() {
            self.coeffs = other.coeffs.to_owned();
            self.min_pow = other.min_pow();
            return;
        }
        let min_pow = self.min_pow().unwrap();
        let other_min_pow = other.min_pow().unwrap();
        if other_min_pow < min_pow {
            self.extend_min((min_pow - other_min_pow) as usize);
        }
        let max_pow = self.max_pow().unwrap();
        let other_max_pow = other.max_pow().unwrap();
        if other_max_pow > max_pow {
            self.extend_max((other_max_pow - max_pow) as usize);
        }
        debug_assert!(self.min_pow() <= other.min_pow());
        debug_assert!(self.max_pow() >= other.max_pow());
        let min_pow = self.min_pow().unwrap();
        for (pow, coeff) in other.into_iter() {
            self.coeffs[(pow - min_pow) as usize] += coeff;
        }
        self.trim();
    }
}

impl<C: Coeff + Clone, Rhs> Add<Rhs> for Polynomial<C>
where
    Polynomial<C>: AddAssign<Rhs>,
{
    type Output = Polynomial<C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, C: Coeff> SubAssign<&'a Polynomial<C>> for Polynomial<C>
where
    for<'c> &'c Polynomial<C>: Neg<Output = Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let res = Polynomial::new(0, vec!());
    /// p -= &p.clone();
    /// assert_eq!(res, p);
    /// ```
    fn sub_assign(&mut self, other: &'a Polynomial<C>) {
        *self += -other;
    }
}

impl<'a, C: Coeff> SubAssign<PolynomialSlice<'a, C>> for Polynomial<C>
where
    for<'c> PolynomialSlice<'c, C>: Neg<Output = Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>>,
{
    fn sub_assign(&mut self, other: PolynomialSlice<'a, C>) {
        *self += -other;
    }
}

impl<C: Coeff> SubAssign<Polynomial<C>> for Polynomial<C>
where
    Polynomial<C>: AddAssign + Neg<Output = Polynomial<C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// let res = Polynomial::new(0, vec!());
    /// p -= p.clone();
    /// assert_eq!(res, p);
    /// ```
    fn sub_assign(&mut self, other: Polynomial<C>) {
        *self += -other;
    }
}

impl<C: Coeff, T> Sub<T> for Polynomial<C>
where
    Polynomial<C>: SubAssign<T>,
{
    type Output = Polynomial<C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, C: Coeff + Clone + AddAssign> MulAssign<&'a Polynomial<C>>
    for Polynomial<C>
where
    Polynomial<C>: MulAssign<PolynomialSlice<'a, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = Polynomial::new(-6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a Polynomial<C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, C: Coeff> MulAssign<PolynomialSlice<'a, C>> for Polynomial<C>
where
    for<'b> PolynomialSlice<'b, C>:
        Mul<PolynomialSlice<'a, C>, Output = Polynomial<C>>,
{
    fn mul_assign(&mut self, other: PolynomialSlice<'a, C>) {
        let prod = self.as_slice(..) * other;
        *self = prod;
    }
}

impl<C: Coeff> MulAssign for Polynomial<C>
where
    for<'a> Polynomial<C>: MulAssign<&'a Polynomial<C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = Polynomial::new(-6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: Polynomial<C>) {
        *self *= &other
    }
}

impl<C: Coeff> MulAssign<C> for Polynomial<C>
where
    for<'a> C: MulAssign<&'a C>,
{
    /// Multiply each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p *= 2.;
    /// let res = Polynomial::new(-3, vec!(2.,0.,-6.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: C) {
        self.mul_assign(&other)
    }
}

impl<'a, C: Coeff> MulAssign<&'a C> for Polynomial<C>
where
    C: MulAssign<&'a C>,
{
    /// Multiply each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p *= &2.;
    /// let res = Polynomial::new(-3, vec!(2.,0.,-6.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a C) {
        for c in &mut self.coeffs {
            *c *= other
        }
    }
}

impl<C: Coeff> DivAssign<C> for Polynomial<C>
where
    for<'a> C: DivAssign<&'a C>,
{
    /// Divide each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p /= 2.;
    /// let res = Polynomial::new(-3, vec!(0.5,0.,-1.5));
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: C) {
        self.div_assign(&other)
    }
}

impl<'a, C: Coeff> DivAssign<&'a C> for Polynomial<C>
where
    C: DivAssign<&'a C>,
{
    /// Divide each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new(-3, vec!(1.,0.,-3.));
    /// p /= &2.;
    /// let res = Polynomial::new(-3, vec!(0.5,0.,-1.5));
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: &'a C) {
        for c in &mut self.coeffs {
            *c /= other
        }
    }
}

impl<C: Coeff> AddAssign<C> for Polynomial<C>
where
    for<'c> Polynomial<C>: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: C) {
        self.add_assign(&other)
    }
}

impl<'a, C: Coeff> AddAssign<&'a C> for Polynomial<C>
where
    Polynomial<C>: AddAssign,
    C: Clone + AddAssign<&'a C>,
{
    fn add_assign(&mut self, other: &'a C) {
        match self.min_pow() {
            None => {
                self.min_pow = Some(0);
                self.coeffs.push(other.clone())
            }
            Some(min_pow) => {
                if min_pow > 0 {
                    self.extend_min(min_pow as usize);
                } else if self.max_pow().unwrap() < 0 {
                    self.extend_max((-self.max_pow().unwrap()) as usize);
                }
                debug_assert!(self.min_pow().unwrap() <= 0);
                debug_assert!(0 <= self.max_pow().unwrap());
                let min_pow = self.min_pow().unwrap();
                self.coeffs[-min_pow as usize] += other;
            }
        }
        self.trim();
    }
}

impl<C: Coeff> SubAssign<C> for Polynomial<C>
where
    Polynomial<C>: AddAssign<C>,
    C: Neg<Output = C> + AddAssign,
{
    fn sub_assign(&mut self, other: C) {
        self.add_assign(-other);
    }
}

impl<'c, C: Coeff> SubAssign<&'c C> for Polynomial<C>
where
    Polynomial<C>: AddAssign<C>,
    C: AddAssign,
    &'c C: Neg<Output = C>,
{
    fn sub_assign(&mut self, other: &'c C) {
        self.add_assign(-other);
    }
}

impl<C: Coeff> Mul for Polynomial<C>
where
    Polynomial<C>: MulAssign,
{
    type Output = Polynomial<C>;

    fn mul(mut self, other: Polynomial<C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a Polynomial<C>> for Polynomial<C>
where
    Polynomial<C>: MulAssign<PolynomialSlice<'a, C>>,
{
    type Output = Polynomial<C>;

    fn mul(self, other: &'a Polynomial<C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, C: Coeff> Mul<PolynomialSlice<'a, C>> for Polynomial<C>
where
    Polynomial<C>: MulAssign<PolynomialSlice<'a, C>>,
{
    type Output = Polynomial<C>;

    fn mul(mut self, other: PolynomialSlice<'a, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff> Mul<C> for Polynomial<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, C: Coeff> Mul<&'a C> for Polynomial<C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff> Div<C> for Polynomial<C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<C>;

    fn div(mut self, other: C) -> Self::Output {
        self /= &other;
        self
    }
}

impl<'a, C: Coeff> Div<&'a C> for Polynomial<C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<C>;

    fn div(mut self, other: &'a C) -> Self::Output {
        self /= other;
        self
    }
}

impl<'a, C: Coeff, T> Mul<T> for &'a Polynomial<C>
where
    PolynomialSlice<'a, C>: Mul<T, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, C: Coeff, T> Div<T> for &'a Polynomial<C>
where
    PolynomialSlice<'a, C>: Div<T, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn div(self, other: T) -> Self::Output {
        self.as_slice(..) / other
    }
}

impl<'a, C: Coeff, T> Add<T> for &'a Polynomial<C>
where
    PolynomialSlice<'a, C>: Add<T, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn add(self, other: T) -> Self::Output {
        self.as_slice(..) + other
    }
}

impl<'a, C: Coeff, T> Sub<T> for &'a Polynomial<C>
where
    PolynomialSlice<'a, C>: Sub<T, Output = Polynomial<C>>,
{
    type Output = Polynomial<C>;

    fn sub(self, other: T) -> Self::Output {
        self.as_slice(..) - other
    }
}

impl<'a, C: Coeff + Clone> From<PolynomialSlice<'a, C>> for Polynomial<C> {
    fn from(s: PolynomialSlice<'a, C>) -> Polynomial<C> {
        Polynomial::new(s.min_pow.unwrap_or(0), s.coeffs.to_vec())
    }
}

impl<C: Coeff> Zero for Polynomial<C>
where
    Polynomial<C>: Add<Output = Polynomial<C>>,
{
    fn zero() -> Self {
        Self {
            min_pow: None,
            coeffs: vec![],
            zero: C::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }
}

impl<C: Coeff> One for Polynomial<C>
where
    Polynomial<C>: Mul<Output = Polynomial<C>>,
{
    fn one() -> Self {
        Self {
            min_pow: Some(0),
            coeffs: vec![C::one()],
            zero: C::zero(),
        }
    }

    fn is_one(&self) -> bool {
        self.min_pow == Some(0)
            && self.coeffs.len() == 1
            && self.coeff(0).is_one()
    }
}

impl<C: Coeff, Var> From<PolynomialIn<Var, C>> for Polynomial<C> {
    fn from(source: PolynomialIn<Var, C>) -> Self {
        let PolynomialInParts { var: _, min_pow, coeffs  } = source.into();
        Self { min_pow, coeffs, zero: C::zero() }
    }
}

/// Data parts of a polynomial
///
/// # Example
///
/// ```rust
/// // destructure a polynomial
/// let p = series::Polynomial::new(-1, vec![1,2,3]);
/// let series::PolynomialParts{min_pow, coeffs} = p.into();
/// assert_eq!(min_pow, Some(-1));
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct PolynomialParts<C> {
    pub min_pow: Option<isize>,
    pub coeffs: Vec<C>,
}

impl<C: Coeff> From<Polynomial<C>> for PolynomialParts<C> {
    fn from(p: Polynomial<C>) -> Self {
        PolynomialParts {
            min_pow: p.min_pow,
            coeffs: p.coeffs,
        }
    }
}

impl<'a, 'b, C: Coeff> KaratsubaMul<&'b Polynomial<C>> for &'a Polynomial<C>
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
        self.as_slice(..).karatsuba_mul(rhs.as_slice(..), min_size)
    }
}

impl<'a, 'b, C: Coeff> KaratsubaMul<PolynomialSlice<'b, C>>
    for &'a Polynomial<C>
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
        self.as_slice(..).karatsuba_mul(rhs, min_size)
    }
}

// dubious helpers trait that only serve to prevent obscure
// compiler errors (rust 1.36.0)
pub(crate) trait MulHelper<'a, 'b, C: Coeff> {
    fn add_prod(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    );

    fn add_prod_naive(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
    );

    fn add_prod_karatsuba(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    );

    fn add_prod_unchecked(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    );

    fn resize_to_fit(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
    );
}

impl<'a, 'b, C: Coeff> MulHelper<'a, 'b, C> for Polynomial<C>
where
    C: 'a + 'b + Clone,
    for<'c> C: AddAssign,
    for<'c> Polynomial<C>:
        AddAssign<&'c Polynomial<C>> + SubAssign<&'c Polynomial<C>>,
    Polynomial<C>: AddAssign<Polynomial<C>> + SubAssign<Polynomial<C>>,
    for<'c> PolynomialSlice<'c, C>:
        Add<Output = Polynomial<C>> + Mul<Output = Polynomial<C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    fn add_prod(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    ) {
        if a.min_pow().is_none() || b.min_pow().is_none() {
            return;
        }
        self.resize_to_fit(a, b);
        self.add_prod_unchecked(a, b, min_karatsuba_size);
        self.trim();
    }

    fn resize_to_fit(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
    ) {
        debug_assert_ne!(a.min_pow(), None);
        debug_assert_ne!(b.min_pow(), None);
        let prod_min_pow = a.min_pow().unwrap() + b.min_pow().unwrap();
        match self.min_pow() {
            Some(min_pow) => {
                if min_pow > prod_min_pow {
                    let num_missing = (min_pow - prod_min_pow) as usize;
                    let to_prepend = iter::repeat(C::zero()).take(num_missing);
                    self.coeffs.splice(0..0, to_prepend);
                    self.min_pow = Some(prod_min_pow);
                }
            }
            None => {
                self.min_pow = Some(prod_min_pow);
            }
        }
        let prod_len = a.len() + b.len() - 1;
        if self.len() < prod_len {
            self.coeffs.resize(prod_len, C::zero());
        }
    }

    fn add_prod_unchecked(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    ) {
        if std::cmp::min(a.len(), b.len()) < min_karatsuba_size {
            trace!(
                "Cauchy product of polynomials with lengths {}, {}",
                a.len(),
                b.len()
            );
            self.add_prod_naive(a, b);
        } else {
            // TODO: split a or b if it's too long?
            trace!(
                "Karatsuba product of polynomials with lengths {}, {}",
                a.len(),
                b.len()
            );
            self.add_prod_karatsuba(a, b, min_karatsuba_size);
        }
    }

    fn add_prod_karatsuba(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
        min_karatsuba_size: usize,
    ) {
        let mid = ((std::cmp::min(a.len(), b.len()) + 1) / 2) as isize;
        let (a_low, mut a_high) = a.split_at(a.min_pow().unwrap() + mid);
        let (b_low, mut b_high) = b.split_at(b.min_pow().unwrap() + mid);
        if let Some(min_pow) = a_high.min_pow.as_mut() {
            *min_pow -= mid
        }
        if let Some(min_pow) = b_high.min_pow.as_mut() {
            *min_pow -= mid
        }
        let t = a_low + a_high;
        let mut u = b_low + b_high;
        if let Some(min_pow) = u.min_pow.as_mut() {
            *min_pow += mid
        }
        self.add_prod_unchecked(
            t.as_slice(..),
            u.as_slice(..),
            min_karatsuba_size,
        );
        let mut t = a_low * b_low;
        *self += &t;
        if let Some(min_pow) = t.min_pow.as_mut() {
            *min_pow += mid
        }
        *self -= t;
        let mut t = a_high * b_high;
        if let Some(min_pow) = t.min_pow.as_mut() {
            *min_pow += mid
        }
        *self -= &t;
        if let Some(min_pow) = t.min_pow.as_mut() {
            *min_pow += mid
        }
        *self += t;
    }

    fn add_prod_naive(
        &mut self,
        a: PolynomialSlice<'a, C>,
        b: PolynomialSlice<'b, C>,
    ) {
        let min_pow = self.min_pow().unwrap();
        for (i, a) in a.iter() {
            for (j, b) in b.iter() {
                self.coeffs[(i + j - min_pow) as usize] += a * b;
            }
        }
    }
}
