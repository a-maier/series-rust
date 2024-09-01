use crate::traits::{AsSlice, KaratsubaMul};
use crate::util::{trim_end, trim_start};
use crate::PolynomialSliceIn;
use crate::SeriesIn;
use crate::{Coeff, IntoIter, Iter};

use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};
use std::{convert, fmt, iter};

use log::trace;

/// Laurent polynomial in a single variable
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct PolynomialIn<Var, C: Coeff> {
    pub(crate) var: Var,
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: Vec<C>,
}

impl<Var, C: Coeff> PolynomialIn<Var, C> {
    /// Create a new Laurent polynomial
    ///
    /// # Example
    ///
    /// This creates a Laurent polynomial in the variable "x", starting
    /// at "x"^-1 with coefficients 1, 2, 3. In other words, the polynomial
    /// x^-1 + 2 + 3*x.
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// ```
    pub fn new(
        var: Var,
        min_pow: isize,
        coeffs: Vec<C>,
    ) -> PolynomialIn<Var, C> {
        let mut res = PolynomialIn {
            var,
            min_pow: Some(min_pow),
            coeffs,
        };
        res.trim();
        res
    }

    /// Get the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        &self.var
    }

    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.min_pow(), Some(-1));
    /// let p = series::PolynomialIn::new("x", -1, vec![0]);
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.max_pow(), Some(1));
    /// let p = series::PolynomialIn::new("x", -1, vec![0]);
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// assert!(!p.is_empty());
    ///
    /// let p = series::PolynomialIn::new("x", -1, vec!(0));
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
    /// let s = series::SeriesIn::with_cutoff("x", -1, 5, vec!(1,2,3));
    /// assert_eq!(p.cutoff_at(5), s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the cutoff power is lower than the starting power
    ///
    pub fn cutoff_at(self, cutoff_pow: isize) -> SeriesIn<Var, C> {
        SeriesIn::with_cutoff(
            self.var,
            self.min_pow.unwrap_or(cutoff_pow),
            cutoff_pow,
            self.coeffs,
        )
    }

    fn as_empty_slice(&self) -> PolynomialSliceIn<'_, Var, C> {
        PolynomialSliceIn::new(
            &self.var,
            0,
            &self.coeffs[self.len()..],
        )
    }

    /// Try to get the coefficient of the polynomial variable to the
    /// given power. Returns [None] if `pow` is less than [min_pow] or
    /// greater than [max_pow].
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec![1,0,3]);
    /// assert_eq!(p.try_coeff(-5), None);
    /// assert_eq!(p.try_coeff(-2), None);
    /// assert_eq!(p.try_coeff(-1), Some(&1));
    /// assert_eq!(p.try_coeff(0), Some(&0));
    /// assert_eq!(p.try_coeff(1), Some(&3));
    /// assert_eq!(p.try_coeff(2), None);
    /// assert_eq!(p.try_coeff(5), None);
    /// ```
    pub fn try_coeff(&self, pow: isize) -> Option<&C> {
        self.as_slice(..).try_coeff(pow)
    }
}

impl<Var, C: 'static + Coeff + Send + Sync> PolynomialIn<Var, C> {
    /// Get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
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
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, Range<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = (r.start - min_pow) as usize;
            let end = (r.end - min_pow) as usize;
            PolynomialSliceIn::new(
                &self.var,
                r.start,
                &self.coeffs[start..end],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeInclusive<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let (start, end) = r.into_inner();
            let ustart = (start - min_pow) as usize;
            let end = (end - min_pow) as usize;
            PolynomialSliceIn::new(
                &self.var,
                start,
                &self.coeffs[ustart..=end],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSliceIn::new(
                &self.var,
                min_pow,
                &self.coeffs[..=end],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = r.start - min_pow;
            PolynomialSliceIn::new(
                &self.var,
                r.start,
                &self.coeffs[(start as usize)..],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSliceIn::new(
                &self.var,
                min_pow,
                &self.coeffs[..end],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFull>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            PolynomialSliceIn::new(
                &self.var,
                min_pow,
                &self.coeffs[r],
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<Var, C: Coeff> convert::From<SeriesIn<Var, C>> for PolynomialIn<Var, C> {
    fn from(s: SeriesIn<Var, C>) -> Self {
        PolynomialIn::new(s.var, s.min_pow, s.coeffs)
    }
}

impl<Var, C: Coeff> Index<isize> for PolynomialIn<Var, C> {
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
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
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

impl<Var, C: Coeff> std::iter::IntoIterator for PolynomialIn<Var, C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -1, vec!(1,2,3));
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

impl<Var, C: Coeff> PolynomialIn<Var, C> {
    fn trim(&mut self) {
        let zero = C::zero();
        trim_end(&mut self.coeffs, &zero);
        if self.coeffs.is_empty() {
            self.min_pow = None;
        } else {
            let min_pow_shift = trim_start(&mut self.coeffs, &zero);
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

impl<Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for PolynomialIn<Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.as_slice(..).fmt(f)
    }
}

impl<Var, C: Coeff + Neg<Output = C>> Neg for PolynomialIn<Var, C> {
    type Output = PolynomialIn<Var, C>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_p = series::PolynomialIn::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        PolynomialIn::new(self.var, self.min_pow.unwrap_or(0), neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a PolynomialIn<Var, C>
where
    PolynomialSliceIn<'a, Var, C>: Neg,
{
    type Output = <PolynomialSliceIn<'a, Var, C> as Neg>::Output;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_p = series::PolynomialIn::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-&p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<&'a PolynomialIn<Var, C>> for PolynomialIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set p = p + q for two Laurent polynomials p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let q = PolynomialIn::new("x", -1, vec!(3., 4., 5.));
    /// let res = PolynomialIn::new("x", -3, vec!(1.,0.,0.,4.,5.));
    /// p += &q;
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn add_assign(&mut self, other: &'a PolynomialIn<Var, C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<PolynomialSliceIn<'a, Var, C>> for PolynomialIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSliceIn<'a, Var, C>) {
        assert_eq!(self.var, *other.var);
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

impl<Var, C: Coeff> AddAssign<PolynomialIn<Var, C>> for PolynomialIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
    Var: PartialEq + fmt::Debug,
    C: Clone + AddAssign,
{
    /// Set p = p + q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let q = PolynomialIn::new("x", -1, vec!(3., 4., 5.));
    /// let res = PolynomialIn::new("x", -3, vec!(1.,0.,0.,4.,5.));
    /// p += q;
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn add_assign(&mut self, other: PolynomialIn<Var, C>) {
        //TODO: code duplication with AddAssign<PolynomialSlice>
        assert_eq!(self.var, other.var);
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

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: AddAssign<Rhs>,
{
    type Output = PolynomialIn<Var, C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, Var, C: Coeff> SubAssign<&'a PolynomialIn<Var, C>>
    for PolynomialIn<Var, C>
where
    for<'c> &'c PolynomialIn<Var, C>: Neg<Output = PolynomialIn<Var, C>>,
    PolynomialIn<Var, C>: AddAssign<PolynomialIn<Var, C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let res = PolynomialIn::new("x", 0, vec!());
    /// p -= &p.clone();
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn sub_assign(&mut self, other: &'a PolynomialIn<Var, C>) {
        *self += -other;
    }
}

impl<'a, Var, C: Coeff> SubAssign<PolynomialSliceIn<'a, Var, C>>
    for PolynomialIn<Var, C>
where
    for<'c> PolynomialSliceIn<'c, Var, C>: Neg<Output = PolynomialIn<Var, C>>,
    PolynomialIn<Var, C>: AddAssign<PolynomialIn<Var, C>>,
{
    fn sub_assign(&mut self, other: PolynomialSliceIn<'a, Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff> SubAssign<PolynomialIn<Var, C>> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: AddAssign + Neg<Output = PolynomialIn<Var, C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// let res = PolynomialIn::new("x", 0, vec!());
    /// p -= p.clone();
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomial have different expansion variables.
    fn sub_assign(&mut self, other: PolynomialIn<Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff, T> Sub<T> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: SubAssign<T>,
{
    type Output = PolynomialIn<Var, C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone + AddAssign>
    MulAssign<&'a PolynomialIn<Var, C>> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: MulAssign<PolynomialSliceIn<'a, Var, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = PolynomialIn::new("x", -6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn mul_assign(&mut self, other: &'a PolynomialIn<Var, C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, Var, C: Coeff> MulAssign<PolynomialSliceIn<'a, Var, C>>
    for PolynomialIn<Var, C>
where
    for<'b> PolynomialSliceIn<'b, Var, C>:
        Mul<PolynomialSliceIn<'a, Var, C>, Output = PolynomialIn<Var, C>>,
{
    fn mul_assign(&mut self, other: PolynomialSliceIn<'a, Var, C>) {
        let prod = self.as_slice(..) * other;
        *self = prod;
    }
}

impl<Var, C: Coeff> MulAssign for PolynomialIn<Var, C>
where
    for<'a> PolynomialIn<Var, C>: MulAssign<&'a PolynomialIn<Var, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = PolynomialIn::new("x", -6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: PolynomialIn<Var, C>) {
        *self *= &other
    }
}

impl<Var, C: Coeff> MulAssign<C> for PolynomialIn<Var, C>
where
    for<'a> C: MulAssign<&'a C>,
{
    /// Multiply each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p *= 2.;
    /// let res = PolynomialIn::new("x", -3, vec!(2.,0.,-6.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: C) {
        self.mul_assign(&other)
    }
}

impl<'a, Var, C: Coeff> MulAssign<&'a C> for PolynomialIn<Var, C>
where
    C: MulAssign<&'a C>,
{
    /// Multiply each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p *= &2.;
    /// let res = PolynomialIn::new("x", -3, vec!(2.,0.,-6.));
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a C) {
        for c in &mut self.coeffs {
            *c *= other
        }
    }
}

impl<Var, C: Coeff> DivAssign<C> for PolynomialIn<Var, C>
where
    for<'a> C: DivAssign<&'a C>,
{
    /// Divide each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p /= 2.;
    /// let res = PolynomialIn::new("x", -3, vec!(0.5,0.,-1.5));
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: C) {
        self.div_assign(&other)
    }
}

impl<'a, Var, C: Coeff> DivAssign<&'a C> for PolynomialIn<Var, C>
where
    C: DivAssign<&'a C>,
{
    /// Divide each monomial by a factor
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::PolynomialIn;
    /// let mut p = PolynomialIn::new("x", -3, vec!(1.,0.,-3.));
    /// p /= &2.;
    /// let res = PolynomialIn::new("x", -3, vec!(0.5,0.,-1.5));
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: &'a C) {
        for c in &mut self.coeffs {
            *c /= other
        }
    }
}

impl<Var, C: Coeff> AddAssign<C> for PolynomialIn<Var, C>
where
    for<'c> PolynomialIn<Var, C>: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: C) {
        self.add_assign(&other)
    }
}

impl<'a, Var, C: Coeff> AddAssign<&'a C> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: AddAssign,
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

impl<Var, C: Coeff> SubAssign<C> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: AddAssign<C>,
    C: Neg<Output = C> + AddAssign,
{
    fn sub_assign(&mut self, other: C) {
        self.add_assign(-other);
    }
}

impl<'c, Var, C: Coeff> SubAssign<&'c C> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: AddAssign<C>,
    C: AddAssign,
    &'c C: Neg<Output = C>,
{
    fn sub_assign(&mut self, other: &'c C) {
        self.add_assign(-other);
    }
}

impl<Var, C: Coeff> Mul for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: MulAssign,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(mut self, other: PolynomialIn<Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a PolynomialIn<Var, C>> for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: MulAssign<PolynomialSliceIn<'a, Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, other: &'a PolynomialIn<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var, C: Coeff> Mul<PolynomialSliceIn<'a, Var, C>>
    for PolynomialIn<Var, C>
where
    PolynomialIn<Var, C>: MulAssign<PolynomialSliceIn<'a, Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(mut self, other: PolynomialSliceIn<'a, Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<Var, C: Coeff> Mul<C> for PolynomialIn<Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a C> for PolynomialIn<Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<Var, C: Coeff> Div<C> for PolynomialIn<Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = PolynomialIn<Var, C>;

    fn div(mut self, other: C) -> Self::Output {
        self /= &other;
        self
    }
}

impl<'a, Var, C: Coeff> Div<&'a C> for PolynomialIn<Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = PolynomialIn<Var, C>;

    fn div(mut self, other: &'a C) -> Self::Output {
        self /= other;
        self
    }
}

impl<'a, Var, C: Coeff, T> Mul<T> for &'a PolynomialIn<Var, C>
where
    PolynomialSliceIn<'a, Var, C>: Mul<T, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, Var, C: Coeff, T> Div<T> for &'a PolynomialIn<Var, C>
where
    PolynomialSliceIn<'a, Var, C>: Div<T, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn div(self, other: T) -> Self::Output {
        self.as_slice(..) / other
    }
}

impl<'a, Var, C: Coeff, T> Add<T> for &'a PolynomialIn<Var, C>
where
    PolynomialSliceIn<'a, Var, C>: Add<T, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn add(self, other: T) -> Self::Output {
        self.as_slice(..) + other
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for &'a PolynomialIn<Var, C>
where
    PolynomialSliceIn<'a, Var, C>: Sub<T, Output = PolynomialIn<Var, C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        self.as_slice(..) - other
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> From<PolynomialSliceIn<'a, Var, C>>
    for PolynomialIn<Var, C>
{
    fn from(s: PolynomialSliceIn<'a, Var, C>) -> PolynomialIn<Var, C> {
        PolynomialIn::new(
            s.var.clone(),
            s.min_pow.unwrap_or(0),
            s.coeffs.to_vec(),
        )
    }
}

/// Data parts of a polynomial
///
/// # Example
///
/// ```rust
/// // destructure a polynomial
/// let p = series::PolynomialIn::new("x", -1, vec![1,2,3]);
/// let series::PolynomialInParts{var, min_pow, coeffs} = p.into();
/// assert_eq!(var, "x");
/// assert_eq!(min_pow, Some(-1));
/// assert_eq!(coeffs, vec![1,2,3]);
/// ```
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct PolynomialInParts<Var, C> {
    pub var: Var,
    pub min_pow: Option<isize>,
    pub coeffs: Vec<C>,
}

impl<Var, C: Coeff> From<PolynomialIn<Var, C>> for PolynomialInParts<Var, C> {
    fn from(p: PolynomialIn<Var, C>) -> Self {
        PolynomialInParts {
            var: p.var,
            min_pow: p.min_pow,
            coeffs: p.coeffs,
        }
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<&'b PolynomialIn<Var, C>>
    for &'a PolynomialIn<Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> PolynomialIn<Var, C>: AddAssign<&'c PolynomialIn<Var, C>>
        + SubAssign<&'c PolynomialIn<Var, C>>,
    PolynomialIn<Var, C>:
        AddAssign<PolynomialIn<Var, C>> + SubAssign<PolynomialIn<Var, C>>,
    for<'c> PolynomialSliceIn<'c, Var, C>: Add<Output = PolynomialIn<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: &'b PolynomialIn<Var, C>,
        min_size: usize,
    ) -> Self::Output {
        self.as_slice(..).karatsuba_mul(rhs.as_slice(..), min_size)
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<PolynomialSliceIn<'b, Var, C>>
    for &'a PolynomialIn<Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    C: Clone,
    for<'c> C: AddAssign,
    for<'c> PolynomialIn<Var, C>: AddAssign<&'c PolynomialIn<Var, C>>
        + SubAssign<&'c PolynomialIn<Var, C>>,
    PolynomialIn<Var, C>:
        AddAssign<PolynomialIn<Var, C>> + SubAssign<PolynomialIn<Var, C>>,
    for<'c> PolynomialSliceIn<'c, Var, C>: Add<Output = PolynomialIn<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: PolynomialSliceIn<'b, Var, C>,
        min_size: usize,
    ) -> Self::Output {
        self.as_slice(..).karatsuba_mul(rhs, min_size)
    }
}

// dubious helpers trait that only serve to prevent obscure
// compiler errors (rust 1.36.0)
pub(crate) trait MulHelper<'a, 'b, Var, C: Coeff> {
    fn add_prod(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
        min_karatsuba_size: usize,
    );

    fn add_prod_naive(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
    );

    fn add_prod_karatsuba(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
        min_karatsuba_size: usize,
    );

    fn add_prod_unchecked(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
        min_karatsuba_size: usize,
    );

    fn resize_to_fit(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
    );
}

impl<'a, 'b, Var, C: Coeff> MulHelper<'a, 'b, Var, C> for PolynomialIn<Var, C>
where
    Var: 'a + 'b,
    C: 'a + 'b + Clone,
    for<'c> C: AddAssign,
    for<'c> PolynomialIn<Var, C>: AddAssign<&'c PolynomialIn<Var, C>>
        + SubAssign<&'c PolynomialIn<Var, C>>,
    PolynomialIn<Var, C>:
        AddAssign<PolynomialIn<Var, C>> + SubAssign<PolynomialIn<Var, C>>,
    for<'c> PolynomialSliceIn<'c, Var, C>:
        Add<Output = PolynomialIn<Var, C>> + Mul<Output = PolynomialIn<Var, C>>,
    for<'c> &'c C: Mul<Output = C>,
{
    fn add_prod(
        &mut self,
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
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
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
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
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
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
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
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
        a: PolynomialSliceIn<'a, Var, C>,
        b: PolynomialSliceIn<'b, Var, C>,
    ) {
        let min_pow = self.min_pow().unwrap();
        for (i, a) in a.iter() {
            for (j, b) in b.iter() {
                self.coeffs[(i + j - min_pow) as usize] += a * b;
            }
        }
    }
}
