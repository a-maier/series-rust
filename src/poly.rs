use crate::polyslice::PolynomialSlice;
use crate::{Coeff, Iter, IterMut, IntoIter};
use crate::traits::AsSlice;
use crate::util::{trim_start, trim_end};

use std::fmt;
use std::iter;
use std::ops::{
    Add, AddAssign, Mul, MulAssign, Neg, Range,
    RangeInclusive, RangeToInclusive, RangeFrom, RangeFull, RangeTo,
    Sub, SubAssign, Index, IndexMut
};

/// Laurent polynomial in a single variable
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct Polynomial<Var, C: Coeff> {
    pub(crate) var: Var,
    pub(crate) min_pow: Option<isize>,
    pub(crate) coeffs: Vec<C>,
    pub(crate) zero: C, // TODO: evil hack, learn how to do this better
}

impl<Var, C: Coeff> Polynomial<Var, C> {
    /// Create a new Laurent polynomial
    ///
    /// # Example
    ///
    /// This creates a Laurent polynomial in the variable "x", starting
    /// at "x"^-1 with coefficients 1, 2, 3. In other words, the polynomial
    /// x^-1 + 2 + 3*x.
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// ```
    pub fn new(var: Var, min_pow: isize, coeffs: Vec<C>) -> Polynomial<Var, C> {
        let mut res = Polynomial{
            var,
            min_pow: Some(min_pow),
            coeffs,
            zero: C::from(0),
        };
        res.trim();
        res
    }

    /// Get the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        &self.var
    }

    /// Get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.coeff(-5), &0);
    /// assert_eq!(p.coeff(-2), &0);
    /// assert_eq!(p.coeff(-1), &1);
    /// assert_eq!(p.coeff(0), &2);
    /// assert_eq!(p.coeff(1), &3);
    /// assert_eq!(p.coeff(2), &0);
    /// assert_eq!(p.coeff(5), &0);
    /// ```
    pub fn coeff(& self, pow: isize) -> & C {
        //TODO: code duplication with PolynomialSlice
        // (can't use a temporary slice because it doesn't live long enough)
        if let Some(min_pow) = self.min_pow() {
            let idx = pow - min_pow;
            if idx < 0 || idx >= self.len() as isize {
                &self.zero
            }
            else {
                &self.coeffs[idx as usize]
            }
        } else {
            &self.zero
        }
    }

    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.min_pow(), Some(-1));
    /// let p = series::Polynomial::new("x", -1, vec![0]);
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
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.max_pow(), Some(2));
    /// let p = series::Polynomial::new("x", -1, vec![0]);
    /// assert_eq!(p.max_pow(), None);
    /// ```
    pub fn max_pow(&self) -> Option<isize> {
        self.min_pow.map(|c| c + self.len() as isize)
    }

    /// Get the number of coefficients
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p.len(), 3);
    /// ```
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// let mut iter = p.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.iter())
    }

    /// An iterator that allows modifying each polynomial coefficient.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// for (_pow, coeff) in p.iter_mut() {
    ///     *coeff += 1
    /// }
    /// let inc = series::Polynomial::new("x", -1, vec!(2,3,4));
    /// assert_eq!(p, inc);
    /// ```
    pub fn iter_mut(&mut self) -> IterMut<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.iter_mut())
    }

    fn as_empty_slice<'a>(&'a self) -> PolynomialSlice<'a, Var, C> {
        PolynomialSlice::new(
            &self.var,
            0,
            &self.coeffs[self.len()..],
            &self.zero,
        )
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, Range<isize>> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = (r.start - min_pow) as usize;
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(
                &self.var,
                r.start,
                &self.coeffs[start..end],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeInclusive<isize>> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let (start, end) = r.into_inner();
            let ustart = (start - min_pow) as usize;
            let end = (end - min_pow) as usize;
            PolynomialSlice::new(
                &self.var,
                start,
                &self.coeffs[ustart..=end],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeToInclusive<isize>> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(
                &self.var,
                min_pow,
                &self.coeffs[..=end],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeFrom<isize>> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let start = r.start - min_pow;
            PolynomialSlice::new(
                &self.var,
                start,
                &self.coeffs[(start as usize)..],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeTo<isize>> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            let end = (r.end - min_pow) as usize;
            PolynomialSlice::new(
                &self.var,
                min_pow,
                &self.coeffs[..end],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff>  AsSlice<'a, RangeFull> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFull) -> Self::Output {
        if let Some(min_pow) = self.min_pow() {
            PolynomialSlice::new(
                &self.var,
                min_pow,
                &self.coeffs[r],
                &self.zero,
            )
        } else {
            self.as_empty_slice()
        }
    }
}

impl<Var, C: Coeff> Index<isize> for Polynomial<Var, C> {
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
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// assert_eq!(p[-1], 1);
    /// assert_eq!(p[0], 2);
    /// assert_eq!(p[1], 3);
    /// assert!(std::panic::catch_unwind(|| p[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| p[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        &self.coeffs[(index-self.min_pow.unwrap()) as usize]
    }
}

impl<Var, C: Coeff> IndexMut<isize> for Polynomial<Var, C> {
    /// Access the (mutable) coefficient of the polynomial
    /// variable to the given power.
    ///
    /// # Panics
    ///
    /// Panics if the index is smaller than the leading power or
    /// bigger than the highest power
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// p[-1] = 0;
    /// assert_eq!(p[-1], 0);
    /// assert_eq!(p[0], 2);
    /// assert_eq!(p[1], 3);
    /// assert!(std::panic::catch_unwind(|| p[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| p[2]).is_err());
    /// ```
    fn index_mut(&mut self, index: isize) -> &mut Self::Output {
        &mut self.coeffs[(index-self.min_pow.unwrap()) as usize]
    }
}

impl<Var, C: Coeff> std::iter::IntoIterator for Polynomial<Var, C> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -1, vec!(1,2,3));
    /// let mut iter = p.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        (self.min_pow().unwrap_or(0)..).zip(self.coeffs.into_iter())
    }
}

impl<Var, C: Coeff> Polynomial<Var, C> {
    fn trim(&mut self) {
        trim_end(&mut self.coeffs, &self.zero);
        if self.coeffs.is_empty() {
            self.min_pow = None;
        }
        else {
            let min_pow_shift = trim_start(&mut self.coeffs, &self.zero);
            self.min_pow.as_mut().map(|m| *m += min_pow_shift as isize);
        }
    }

    fn extend_min(&mut self, extend: usize) {
        debug_assert!(self.min_pow != None);
        let to_insert = iter::repeat_with(|| C::from(0)).take(extend);
        self.coeffs.splice(0..0, to_insert);
        self.min_pow.as_mut().map(|m| *m -= extend as isize);
    }

    fn extend_max(&mut self, extend: usize) {
        let to_insert = iter::repeat_with(|| C::from(0)).take(extend);
        self.coeffs.extend(to_insert);
    }
}

impl<Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display
    for Polynomial<Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.as_slice(..).fmt(f)
    }
}

impl<Var, C: Coeff + Neg<Output = C>> Neg for Polynomial<Var, C> {
    type Output = Polynomial<Var, C>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_p = series::Polynomial::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        Polynomial::new(self.var, self.min_pow.unwrap_or(0), neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a Polynomial<Var, C>
where
   PolynomialSlice<'a, Var, C>: Neg
{
    type Output = <PolynomialSlice<'a, Var, C> as Neg>::Output;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = series::Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_p = series::Polynomial::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-&p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<&'a Polynomial<Var, C>> for Polynomial<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    /// Set p = p + q for two Laurent polynomials p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let q = Polynomial::new("x", -1, vec!(3., 4., 5.));
    /// let res = Polynomial::new("x", -3, vec!(1.,0.,0.,4.,5.));
    /// p += &q;
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn add_assign(&mut self, other: &'a Polynomial<Var, C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<PolynomialSlice<'a, Var, C>> for Polynomial<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        assert_eq!(self.var, *other.var);
        if other.min_pow() == None {
            return;
        }
        if self.min_pow() == None {
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
        for (pow, coeff) in other.iter() {
            self[pow] += coeff;
        }
        self.trim();
    }
}

impl<Var, C: Coeff> AddAssign<Polynomial<Var, C>> for Polynomial<Var, C>
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
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let q = Polynomial::new("x", -1, vec!(3., 4., 5.));
    /// let res = Polynomial::new("x", -3, vec!(1.,0.,0.,4.,5.));
    /// p += q;
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn add_assign(&mut self, other: Polynomial<Var, C>) {
        //TODO: code duplication with AddAssign<PolynomialSlice>
        assert_eq!(self.var, other.var);
        if other.min_pow() == None {
            return;
        }
        if self.min_pow() == None {
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
        for (pow, coeff) in other.into_iter() {
            self[pow] += coeff;
        }
        self.trim();
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for &'a Polynomial<Var, C>
where
    Polynomial<Var, C>: AddAssign<Rhs>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs> for Polynomial<Var, C>
where
    Polynomial<Var, C>: AddAssign<Rhs>,
{
    type Output = Polynomial<Var, C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, Var, C: Coeff> SubAssign<&'a Polynomial<Var, C>> for Polynomial<Var, C>
where
    for<'c> &'c Polynomial<Var, C>: Neg<Output = Polynomial<Var, C>>,
    Polynomial<Var, C>: AddAssign<Polynomial<Var, C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Polynomial::new("x", 0, vec!());
    /// p -= &p.clone();
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn sub_assign(&mut self, other: &'a Polynomial<Var, C>) {
        *self += -other;
    }
}

impl<'a, Var, C: Coeff> SubAssign<PolynomialSlice<'a, Var, C>> for Polynomial<Var, C>
where
    for<'c> PolynomialSlice<'c, Var, C>: Neg<Output = Polynomial<Var, C>>,
    Polynomial<Var, C>: AddAssign<Polynomial<Var, C>>,
{
    fn sub_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff> SubAssign<Polynomial<Var, C>> for Polynomial<Var, C>
where
    Polynomial<Var, C>: AddAssign + Neg<Output = Polynomial<Var, C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Polynomial::new("x", 0, vec!());
    /// p -= p.clone();
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomial have different expansion variables.
    fn sub_assign(&mut self, other: Polynomial<Var, C>) {
        *self += -other;
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for &'a Polynomial<Var, C>
where
    Polynomial<Var, C>: Clone + SubAssign<T>,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        let mut res = self.clone();
        res -= other;
        res
    }
}

impl<Var, C: Coeff, T> Sub<T> for Polynomial<Var, C>
where
    Polynomial<Var, C>: SubAssign<T>,
{
    type Output = Polynomial<Var, C>;

    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone + AddAssign>
    MulAssign<&'a Polynomial<Var, C>> for Polynomial<Var, C>
where
    Polynomial<Var, C>: MulAssign<PolynomialSlice<'a, Var, C>>
{
    /// Set p = p * q for two polynomial p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials have different expansion variables.
    fn mul_assign(&mut self, other: &'a Polynomial<Var, C>) {
        self.mul_assign(other.as_slice(..))
    }
}

impl<'a, Var, C: Coeff> MulAssign<PolynomialSlice<'a, Var, C>> for Polynomial<Var, C>
where
    for <'b> PolynomialSlice<'b, Var, C>: Mul<PolynomialSlice<'a, Var, C>, Output=Polynomial<Var,C>>
{
    fn mul_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        let prod = self.as_slice(..) * other;
        *self = prod;
    }
}

impl<Var, C: Coeff> MulAssign for Polynomial<Var, C>
where
    for<'a> Polynomial<Var, C>: MulAssign<&'a Polynomial<Var, C>>,
{
    /// Set p = p * q for two series p,q
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1.,0.,-3.));
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec!(1.,0.,-6.,0.,9.));
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(&mut self, other: Polynomial<Var, C>) {
        *self *= &other
    }
}

const MIN_KARATSUBA_SIZE: usize = 8;

// dubious helpers trait that only serve to prevent obscure
// compiler errors (rust 1.36.0)
pub(crate) trait MulHelper<'a, 'b, Var, C: Coeff> {

    fn add_prod(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>);

    fn add_prod_naive(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>);

    fn add_prod_karatsuba(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>);

    fn add_prod_unchecked(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>);

    fn resize_to_fit(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>);
}

impl<'a, 'b, Var, C: Coeff> MulHelper<'a, 'b, Var, C> for Polynomial<Var, C>
where
    Var: 'a + 'b,
    C: 'a + 'b + Clone,
    for <'c> C: AddAssign,
    for <'c> Polynomial<Var, C>: AddAssign<&'c Polynomial<Var, C>> + SubAssign<&'c Polynomial<Var, C>>,
    Polynomial<Var, C>: AddAssign<Polynomial<Var, C>> + SubAssign<Polynomial<Var, C>>,
    for<'c> PolynomialSlice<'c, Var, C>: Add<Output = Polynomial<Var, C>> + Mul<Output = Polynomial<Var, C>>,
    for <'c> &'c C: Mul<Output=C>
{

    fn add_prod(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>) {
        if a.min_pow() == None || b.min_pow() == None {
            return;
        }
        self.resize_to_fit(a, b);
        self.add_prod_unchecked(a, b);
        self.trim();
    }

    fn resize_to_fit(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>) {
        debug_assert_ne!(a.min_pow(), None);
        debug_assert_ne!(b.min_pow(), None);
        let prod_min_pow = a.min_pow().unwrap() + b.min_pow().unwrap();
        match self.min_pow() {
            Some(min_pow) => {
                if min_pow > prod_min_pow {
                    let num_missing = (min_pow - prod_min_pow) as usize;
                    let to_prepend = iter::repeat(C::from(0)).take(num_missing);
                    self.coeffs.splice(0..0, to_prepend);
                    self.min_pow = Some(prod_min_pow);
                }
            },
            None => {
                self.min_pow = Some(prod_min_pow);
            }
        }
        let prod_len = a.len() + b.len() - 1;
        if self.len() < prod_len {
            self.coeffs.resize(prod_len, C::from(0));
        }
    }

    fn add_prod_unchecked(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>) {
        if a.len() > b.len() {
            return self.add_prod_unchecked(b, a)
        }
        if a.len() < MIN_KARATSUBA_SIZE {
            self.add_prod_naive(a, b);
        }
        else {
            // TODO: split b if it's too long?
            self.add_prod_karatsuba(a, b);
        }
    }

    fn add_prod_karatsuba(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>) {
        debug_assert!(a.len() <= b.len());
        let mid = ((a.len() + 1)/2) as isize;
        let (a_low, mut a_high) = a.split_at(a.min_pow().unwrap() + mid);
        let (b_low, mut b_high) = b.split_at(b.min_pow().unwrap() + mid);
        a_high.min_pow.as_mut().map(|m| *m -= mid);
        b_high.min_pow.as_mut().map(|m| *m -= mid);
        let t = a_low + a_high;
        let mut u = b_low + b_high;
        u.min_pow.as_mut().map(|m| *m += mid);
        self.add_prod_unchecked(t.as_slice(..), u.as_slice(..));
        let mut t = a_low * b_low;
        *self += &t;
        t.min_pow.as_mut().map(|m| *m += mid);
        *self -= t;
        let mut t = a_high * b_high;
        t.min_pow.as_mut().map(|m| *m += mid);
        *self -= &t;
        t.min_pow.as_mut().map(|m| *m += mid);
        *self += t;
    }

    fn add_prod_naive(&mut self, a: PolynomialSlice<'a, Var, C>, b: PolynomialSlice<'b, Var, C>) {
        for (i, a) in a.iter() {
            for (j, b) in b.iter() {
                self[i + j] += a*b;
            }
        }
    }
}

impl<'a, Var, C: Coeff, T> Mul<T> for &'a Polynomial<Var, C>
where
    PolynomialSlice<'a, Var, C>: Mul<T>,
{
    type Output = <PolynomialSlice<'a, Var, C> as Mul<T>>::Output;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<Var, C: Coeff, T> Mul<T> for Polynomial<Var, C>
where
    for<'a> PolynomialSlice<'a, Var, C>: Mul<T, Output=Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}
