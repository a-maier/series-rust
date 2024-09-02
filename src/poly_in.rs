use crate::traits::{AsSlice, KaratsubaMul};
use crate::{Polynomial, PolynomialSlice, PolynomialSliceIn, SeriesInParts};
use crate::SeriesIn;
use crate::{Coeff, IntoIter, Iter};

use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};
use std::{convert, fmt};

/// Laurent polynomial in a single variable
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct PolynomialIn<Var, C: Coeff> {
    pub(crate) var: Var,
    pub(crate) poly: Polynomial<C>
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
        let poly = Polynomial::new(min_pow, coeffs);
        Self {
            var,
            poly,
        }
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
        self.poly.min_pow()
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
        self.poly.max_pow()
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
        self.poly.len()
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
        self.poly.is_empty()
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
        self.poly.iter()
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
        self.poly.cutoff_at(cutoff_pow).in_var(self.var)
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
        self.poly.try_coeff(pow)
    }

    /// Transform all coefficients
    ///
    /// `f(p, c)` is applied to each monomial, where `p` is the power
    /// of the variable and `c` a mutable reference to the
    /// coefficient. `p` takes all values in the range
    /// `min_pow()..=max_pow()`.
    ///
    /// # Example
    ///
    /// Replace each coefficient by its square
    /// ```rust
    /// let mut s = series::PolynomialIn::new("x", -1, vec!(1,2,3,4));
    /// s.for_each(|_, c| *c *= *c);
    /// assert_eq!(s.coeff(-1), &1);
    /// assert_eq!(s.coeff(0), &4);
    /// assert_eq!(s.coeff(1), &9);
    /// assert_eq!(s.coeff(2), &16);
    /// ```
    pub fn for_each<F>(&mut self, f: F)
    where
        F: FnMut(isize, &mut C)
    {
        self.poly.for_each(f)
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
        self.poly.coeff(pow)
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, Range<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeInclusive<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFull>
    for PolynomialIn<Var, C>
{
    type Output = PolynomialSliceIn<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        self.poly.as_slice(r).in_var(self.var())
    }
}

impl<Var, C: Coeff> convert::From<SeriesIn<Var, C>> for PolynomialIn<Var, C> {
    fn from(s: SeriesIn<Var, C>) -> Self {
        let SeriesInParts { var, min_pow, coeffs } = s.into();
        PolynomialIn::new(var, min_pow, coeffs)
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
        &self.poly[index]
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
        self.poly.into_iter()
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
        (-self.poly).in_var(self.var)
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
        assert_eq!(self.var(), other.var());
        self.poly.add_assign(&other.poly)
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<PolynomialSliceIn<'a, Var, C>> for PolynomialIn<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSliceIn<'a, Var, C>) {
        assert_eq!(self.var(), other.var());
        self.poly.add_assign(other.poly)
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
        assert_eq!(self.var(), other.var());
        self.poly.add_assign(other.poly)
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
    Polynomial<C>: MulAssign<&'a Polynomial<C>>,
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
        assert_eq!(self.var, other.var);
        self.poly.mul_assign(&other.poly)
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
    fn mul_assign(&mut self, rhs: &'a C) {
        self.poly.mul_assign(rhs)
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
    fn div_assign(&mut self, rhs: &'a C) {
        self.poly.div_assign(rhs)
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
    C: Clone + AddAssign<&'a C>,
{
    fn add_assign(&mut self, rhs: &'a C) {
        self.poly.add_assign(rhs)
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
        Polynomial::from(s.poly).in_var(s.var.clone())
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
        let PolynomialIn { var, poly } = p;
        PolynomialInParts {
            var,
            min_pow: poly.min_pow,
            coeffs: poly.coeffs,
        }
    }
}

impl<'a, 'b, Var, C: Coeff> KaratsubaMul<&'b PolynomialIn<Var, C>>
    for &'a PolynomialIn<Var, C>
where
    Var: Clone + fmt::Debug + PartialEq,
    &'a Polynomial<C>: KaratsubaMul<&'b Polynomial<C>, Output = Polynomial<C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: &'b PolynomialIn<Var, C>,
        min_size: usize,
    ) -> Self::Output {
        assert_eq!(self.var(), rhs.var());
        self.poly.karatsuba_mul(&rhs.poly, min_size).in_var(self.var.clone())
    }
}

impl<'a, 'b, Var: Clone, C: Coeff> KaratsubaMul<PolynomialSliceIn<'b, Var, C>>
    for &'a PolynomialIn<Var, C>
where
    Var: Clone + PartialEq + fmt::Debug,
    &'a Polynomial<C>: KaratsubaMul<PolynomialSlice<'b, C>, Output = Polynomial<C>>,
{
    type Output = PolynomialIn<Var, C>;

    fn karatsuba_mul(
        self,
        rhs: PolynomialSliceIn<'b, Var, C>,
        min_size: usize,
    ) -> Self::Output {
        assert_eq!(self.var(), rhs.var());
        self.poly.karatsuba_mul(rhs.poly, min_size).in_var(self.var.clone())
    }
}
