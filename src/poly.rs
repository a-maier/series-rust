use crate::traits::AsSlice;
use crate::util::{NumDisplay, trim_slice_zero, trim_zero};
use crate::zero_ref::zero_ref;
use crate::{Coeff, IntoIter, Series, SeriesParts};

use core::slice;
use std::fmt::Display;
use std::iter::FusedIterator;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Range,
    RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive, Sub,
    SubAssign,
};
use std::{fmt::Debug, iter};

use num_traits::{One, Zero};

/// Laurent polynomial in a single variable
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub enum Polynomial<Var, C> {
    Const(C),
    Poly(NonConstPoly<Var, C>),
}

/// A non-constant polynomial
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub struct NonConstPoly<Var, C> {
    min_pow: isize,
    coeffs: Vec<C>,
    var: Var,
}

impl<Var, C: Coeff> NonConstPoly<Var, C> {
    /// Turn a polynomial into a series with the given cutoff
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{Polynomial, Series};
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let Polynomial::Poly(p) = p else {
    ///    unreachable!("Polynomial is not a constant")
    /// };
    /// let s = Series::with_cutoff("x", -1..5, vec![1, 2, 3]);
    /// assert_eq!(p.cutoff_at(5), s);
    /// ```
    pub fn cutoff_at(self, cutoff_pow: isize) -> Series<Var, C> {
        let Self {
            min_pow,
            coeffs,
            var,
        } = self;
        Series::with_cutoff(var, min_pow..cutoff_pow, coeffs)
    }

    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    pub fn var(&self) -> &Var {
        &self.var
    }
}

impl<Var, C: Coeff> Display for Polynomial<Var, C>
where
    for<'c> PolynomialSlice<'c, Var, C>: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_slice(..).fmt(f)
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
/// Data parts of a polynomial
///
/// # Example
///
/// ```rust
/// // destructure a polynomial
/// # use series::{Polynomial, PolynomialParts};
/// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
/// let Polynomial::Poly(p) = p else {
///    unreachable!("Polynomial is not a constant")
/// };
/// let PolynomialParts{min_pow, coeffs, var} = p.into();
/// assert_eq!(min_pow, -1);
/// assert_eq!(coeffs, vec![1, 2, 3]);
/// assert_eq!(var, "x");
/// ```
pub struct PolynomialParts<Var, C> {
    pub min_pow: isize,
    pub coeffs: Vec<C>,
    pub var: Var,
}

impl<Var, C> From<NonConstPoly<Var, C>> for PolynomialParts<Var, C> {
    fn from(value: NonConstPoly<Var, C>) -> Self {
        let NonConstPoly {
            min_pow,
            coeffs,
            var,
        } = value;
        Self {
            min_pow,
            coeffs,
            var,
        }
    }
}

impl<Var, C: Coeff> Polynomial<Var, C> {
    /// Create a new Laurent polynomial
    ///
    /// # Example
    ///
    /// This creates a Laurent polynomial in x, starting at power -1 with
    /// coefficients 1, 2, 3. In other words, the polynomial x^-1 + 2 + 3*x.
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// ```
    pub fn new(var: Var, min_pow: isize, coeffs: Vec<C>) -> Polynomial<Var, C> {
        let mut res = Self::Poly(NonConstPoly {
            min_pow,
            coeffs,
            var,
        });
        res.trim();
        res
    }

    fn trim(&mut self) {
        if let Self::Poly(NonConstPoly {
            min_pow, coeffs, ..
        }) = self
        {
            let removed_from_start = trim_zero(coeffs);
            *min_pow += removed_from_start as isize;
            if coeffs.len() == 1 && *min_pow == 0 {
                *self = Self::Const(coeffs.pop().unwrap());
            } else if coeffs.is_empty() {
                *self = Self::zero();
            };
        }
    }

    /// Create a constant polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let two: Polynomial<(), _> = Polynomial::from_const(2);
    /// ```
    pub const fn from_const(coeff: C) -> Polynomial<Var, C> {
        Self::Const(coeff)
    }

    /// Create the zero polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let zero: Polynomial<(), _> = Polynomial::zero();
    /// assert_eq!(zero, Polynomial::from_const(0));
    /// ```
    pub fn zero() -> Self {
        Self::Const(C::zero())
    }

    /// Check if the polynomial is zero
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// type IntPoly = Polynomial<(), i32>;
    /// assert!(IntPoly::zero().is_zero());
    /// assert!(!IntPoly::one().is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        if let Self::Const(c) = self
            && c.is_zero()
        {
            true
        } else {
            false
        }
    }

    /// Create the unit polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let one: Polynomial<(), _> = Polynomial::one();
    /// assert_eq!(one, Polynomial::from_const(1));
    /// ```
    pub fn one() -> Self {
        Self::Const(C::one())
    }

    /// Check if the polynomial is one
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// type IntPoly = Polynomial<(), i32>;
    /// assert!(IntPoly::one().is_one());
    /// assert!(!IntPoly::zero().is_one());
    /// ```
    pub fn is_one(&self) -> bool {
        if let Self::Const(c) = self
            && c.is_one()
        {
            true
        } else {
            false
        }
    }

    /// Get the leading power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.min_pow(), Some(-1));
    ///
    /// type IntPoly = Polynomial<(), i32>;
    /// let p = IntPoly::from_const(1);
    /// assert_eq!(p.min_pow(), Some(0));
    ///
    /// let p = IntPoly::zero();
    /// assert_eq!(p.min_pow(), None);
    /// ```
    pub fn min_pow(&self) -> Option<isize> {
        self.as_slice(..).min_pow()
    }

    /// Get the highest power of the polynomial variable
    ///
    /// For vanishing polynomials `None` is returned
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.max_pow(), Some(1));
    ///
    /// let p: Polynomial<(), i32> = Polynomial::zero();
    /// assert_eq!(p.max_pow(), None);
    /// ```
    pub fn max_pow(&self) -> Option<isize> {
        self.as_slice(..).max_pow()
    }

    /// Get the difference between the highest and the lowest power of
    /// the polynomial variable
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1 ,2, 3]);
    /// assert_eq!(p.len(), 3);
    ///
    /// let p: Polynomial<(), i32> = Polynomial::zero();
    /// assert_eq!(p.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.as_slice(..).len()
    }

    /// Check if the polynomial is zero
    ///
    /// See [is_zero].
    pub fn is_empty(&self) -> bool {
        self.is_zero()
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let mut iter = p.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    ///
    /// let p: Polynomial<(), i32> = Polynomial::zero();
    /// ```
    pub fn iter(&self) -> Iter<'_, C> {
        self.as_slice(..).iter()
    }

    /// Try to get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// Similar to [coeff](Self::coeff), but returns [None] if `pow`
    /// is less than [min_pow](Self::min_pow) or greater than
    /// [max_pow](Self::max_pow).
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 0, 3]);
    /// assert_eq!(p.get_coeff(-5), None);
    /// assert_eq!(p.get_coeff(-2), None);
    /// assert_eq!(p.get_coeff(-1), Some(&1));
    /// assert_eq!(p.get_coeff(0), Some(&0));
    /// assert_eq!(p.get_coeff(1), Some(&3));
    /// assert_eq!(p.get_coeff(2), None);
    /// assert_eq!(p.get_coeff(5), None);
    /// ```
    pub fn get_coeff(&self, pow: isize) -> Option<&C> {
        self.as_slice(..).get_coeff(pow)
    }

    /// Transform all coefficients
    ///
    /// `f(p, c)` is applied to each monomial, where `p` is the power
    /// of the variable and `c` the coefficient. `p` takes all values
    /// in the range `min_pow()..=max_pow()`.
    ///
    /// # Example
    ///
    /// Replace each coefficient by its square
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3, 4]);
    /// let p = p.map(|_, c| c * c);
    /// assert_eq!(p.coeff(-1), &1);
    /// assert_eq!(p.coeff(0), &4);
    /// assert_eq!(p.coeff(1), &9);
    /// assert_eq!(p.coeff(2), &16);
    /// ```
    pub fn map<D, F>(self, mut f: F) -> Polynomial<Var, D>
    where
        F: FnMut(isize, C) -> D,
        D: Coeff,
    {
        match self {
            Polynomial::Const(c) => Polynomial::Const(f(0, c)),
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let coeffs = coeffs
                    .into_iter()
                    .enumerate()
                    .map(|(n, c)| f(min_pow + n as isize, c))
                    .collect();
                Polynomial::new(var, min_pow, coeffs)
            }
        }
    }

    /// Get the polynomial variable
    ///
    /// Returns `None` if the polynomial is constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.var(), Some(&"x"));
    /// let p: Polynomial<&str, _> = Polynomial::from_const(2);
    /// assert!(p.var().is_none());
    /// ```
    pub fn var(&self) -> Option<&Var> {
        self.as_slice(..).var()
    }

    /// Replace the polynomial variable
    ///
    /// Returns the new polynomial and, if the polynomial was not a
    /// constant, the old variable.
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let (p, var) = p.replace_var("y");
    /// assert_eq!(p.var(), Some(&"y"));
    /// assert_eq!(var, Some("x"));
    ///
    /// let p: Polynomial<&str, _> = Polynomial::from_const(2);
    /// let (p, var) = p.replace_var("y");
    /// assert!(var.is_none());
    /// assert!(p.var().is_none());
    /// ```
    pub fn replace_var<W>(self, new_var: W) -> (Polynomial<W, C>, Option<Var>) {
        match self {
            Polynomial::Const(c) => (Polynomial::Const(c), None),
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => (
                Polynomial::Poly(NonConstPoly {
                    min_pow,
                    coeffs,
                    var: new_var,
                }),
                Some(var),
            ),
        }
    }

    /// Check if the polynomial is constant
    pub const fn is_const(&self) -> bool {
        matches!(self, Polynomial::Const(..))
    }
}

impl<Var: Clone + Debug + PartialEq, C: Coeff> Polynomial<Var, C> {
    /// Turn a polynomial into a series with the given cutoff
    ///
    /// Since constant polynomials do not store the expansion variable
    /// it has to be specified. See [NonConstPoly::cutoff_at] for the
    /// case where we know that the polynomial is not a constant.
    ///
    /// # Panics
    ///
    /// Panics if the passed expansion variable does not agree with
    /// the polynomial variable.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::{Polynomial, Series};
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let Polynomial::Poly(p) = p else {
    ///    unreachable!("Polynomial is not a constant")
    /// };
    /// let s = Series::with_cutoff("x", -1..5, vec![1, 2, 3]);
    /// assert_eq!(p.cutoff_at(5), s);
    /// ```
    pub fn cutoff_at(self, var: &Var, cutoff_pow: isize) -> Series<Var, C> {
        match self {
            Polynomial::Const(c) => {
                Series::with_cutoff(var.to_owned(), 0..cutoff_pow, vec![c])
            }
            Polynomial::Poly(poly) => {
                assert_eq!(var, poly.var());
                poly.cutoff_at(cutoff_pow)
            }
        }
    }
}

impl<Var, C: 'static + Coeff + Send + Sync> Polynomial<Var, C> {
    /// Get the coefficient of the polynomial variable to the
    /// given power.
    ///
    /// Returns a reference to zero if `pow` is less than
    /// [min_pow](Self::min_pow) or greater than
    /// [max_pow](Self::max_pow).
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
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

impl<Var, C: Coeff> Default for Polynomial<Var, C> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<'a, Var: 'a, C: 'static + Coeff + Send + Sync> AsSlice<'a, Range<isize>>
    for Polynomial<Var, C>
{
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: Range<isize>) -> Self::Output {
        match self {
            Polynomial::Const(c) => {
                if r.is_empty() {
                    PolynomialSlice::zero()
                } else if r.start != 0 || r.end != 0 || c.is_zero() {
                    panic!("index out of bounds");
                } else {
                    PolynomialSlice::Const(c)
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let [start, end] =
                    [r.start, r.end].map(|c| (c - *min_pow) as usize);
                PolynomialSlice::new(*min_pow, &coeffs[start..end], var)
            }
        }
    }
}

impl<'a, Var: 'a, C: 'static + Coeff + Send + Sync>
    AsSlice<'a, RangeInclusive<isize>> for Polynomial<Var, C>
{
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeInclusive<isize>) -> Self::Output {
        match self {
            Polynomial::Const(c) => {
                if r.is_empty() {
                    PolynomialSlice::zero()
                } else if *r.start() != 0 || *r.end() != 0 || c.is_zero() {
                    panic!("index out of bounds");
                } else {
                    PolynomialSlice::Const(c)
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let [start, end] =
                    [r.start(), r.end()].map(|c| (c - *min_pow) as usize);
                PolynomialSlice::new(*min_pow, &coeffs[start..=end], var)
            }
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeToInclusive<isize>>
    for Polynomial<Var, C>
{
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeToInclusive<isize>) -> Self::Output {
        match self {
            Polynomial::Const(c) => {
                if r.end != 0 || c.is_zero() {
                    panic!("index out of bounds");
                } else {
                    PolynomialSlice::Const(c)
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let r = ..=(r.end - *min_pow) as usize;
                PolynomialSlice::new(*min_pow, &coeffs[r], var)
            }
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFrom<isize>>
    for Polynomial<Var, C>
{
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeFrom<isize>) -> Self::Output {
        match self {
            Polynomial::Const(c) => {
                if r.start != 0 || c.is_zero() {
                    panic!("index out of bounds");
                } else {
                    PolynomialSlice::Const(c)
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let r = (r.start - *min_pow) as usize..;
                PolynomialSlice::new(*min_pow, &coeffs[r], var)
            }
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeTo<isize>>
    for Polynomial<Var, C>
{
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, r: RangeTo<isize>) -> Self::Output {
        match self {
            Polynomial::Const(c) => {
                if r.end != 1 || c.is_zero() {
                    panic!("index out of bounds");
                } else {
                    PolynomialSlice::Const(c)
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => {
                let r = ..(r.end - *min_pow) as usize;
                PolynomialSlice::new(*min_pow, &coeffs[r], var)
            }
        }
    }
}

impl<'a, Var: 'a, C: 'a + Coeff> AsSlice<'a, RangeFull> for Polynomial<Var, C> {
    type Output = PolynomialSlice<'a, Var, C>;

    fn as_slice(&'a self, _: RangeFull) -> Self::Output {
        match self {
            Polynomial::Const(c) => PolynomialSlice::Const(c),
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var,
            }) => PolynomialSlice::Poly {
                min_pow: *min_pow,
                coeffs,
                var,
            },
        }
    }
}

impl<Var, C: Coeff> From<Series<Var, C>> for Polynomial<Var, C> {
    fn from(s: Series<Var, C>) -> Self {
        let SeriesParts {
            var,
            min_pow,
            coeffs,
        } = s.into();
        Polynomial::new(var, min_pow, coeffs)
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
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p[-1], 1);
    /// assert_eq!(p[0], 2);
    /// assert_eq!(p[1], 3);
    /// assert!(std::panic::catch_unwind(|| p[-2]).is_err());
    /// assert!(std::panic::catch_unwind(|| p[2]).is_err());
    /// ```
    fn index(&self, index: isize) -> &Self::Output {
        match self {
            Polynomial::Const(c) => {
                if c.is_zero() || index == 0 {
                    panic!("index {index} out of bounds");
                } else {
                    c
                }
            }
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var: _,
            }) => &coeffs[(index - min_pow) as usize],
        }
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
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let mut iter = p.into_iter();
    /// assert_eq!(iter.next(), Some((-1, 1)));
    /// assert_eq!(iter.next(), Some((0, 2)));
    /// assert_eq!(iter.next(), Some((1, 3)));
    /// assert_eq!(iter.next(), None);
    /// ```
    fn into_iter(self) -> IntoIter<C> {
        match self {
            // TODO: avoid the vec allocation?
            Polynomial::Const(c) => (0..).zip(vec![c]),
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var: _,
            }) => (min_pow..).zip(coeffs),
        }
    }
}

fn extend_to_range<C: Coeff>(
    coeffs: &mut Vec<C>,
    min_pow: &mut isize,
    pow_range: Range<isize>,
) {
    let min_pow_diff = *min_pow - pow_range.start;
    if min_pow_diff > 0 {
        extend_min(coeffs, min_pow, min_pow_diff as usize);
    }
    let max_pow = *min_pow + coeffs.len() as isize;
    let max_pow_diff = pow_range.end - max_pow;
    if max_pow_diff > 0 {
        extend_max(coeffs, max_pow_diff as usize);
    }
    debug_assert!(*min_pow <= pow_range.start);
    debug_assert!(*min_pow + coeffs.len() as isize >= pow_range.end);
}

fn extend_min<C: Coeff>(
    coeffs: &mut Vec<C>,
    min_pow: &mut isize,
    extend: usize,
) {
    let to_insert = iter::repeat_with(C::zero).take(extend);
    coeffs.splice(0..0, to_insert);
    *min_pow -= extend as isize
}

fn extend_max<C: Coeff>(coeffs: &mut Vec<C>, extend: usize) {
    let to_insert = iter::repeat_with(C::zero).take(extend);
    coeffs.extend(to_insert);
}

impl<Var, C: Coeff + Neg> Neg for Polynomial<Var, C>
where
    <C as Neg>::Output: Coeff,
{
    type Output = Polynomial<Var, <C as Neg>::Output>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let minus_p = Polynomial::new("x", -3, vec![-1, 0, 3]);
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.map(|_, c| -c)
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a Polynomial<Var, C>
where
    &'a C: Neg,
    <&'a C as Neg>::Output: Coeff,
{
    type Output = Polynomial<Var, <&'a C as Neg>::Output>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let minus_p = Polynomial::new("x", -3, vec![-1, 0, 3]);
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.as_slice(..).neg()
    }
}

macro_rules! impl_add_assign_const {
    ($t:ty) => {
        impl<'a, Var, C: Coeff> AddAssign<$t> for Polynomial<Var, C>
        where C: AddAssign<$t>
        {
            /// Add a constant to the polynomial
            ///
            /// # Example
            ///
            /// ```rust
            /// # use series::Polynomial;
            /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
            /// p += 1;
            /// p += &1;
            /// let res = Polynomial::new("x", -3, vec![1, 0, -3, 2]);
            /// assert_eq!(res, p);
            /// ```
            fn add_assign(&mut self, other: $t) {
                if other.is_zero() {
                    return;
                }
                match self {
                    Polynomial::Const(c) => c.add_assign(other),
                    Polynomial::Poly(NonConstPoly {
                        min_pow,
                        coeffs,
                        var: _,
                    }) => {
                        extend_to_range(coeffs, min_pow, 0..1);
                        let pos = (-*min_pow) as usize;
                        coeffs[pos].add_assign(other);
                        self.trim();
                    }
                }
            }
        }
    };
}

impl_add_assign_const!(C);
impl_add_assign_const!(&'a C);

macro_rules! impl_sub_assign_const {
    ($t:ty) => {
        impl<'a, Var, C: Coeff> SubAssign<$t> for Polynomial<Var, C>
        where C: SubAssign<$t>
        {
            /// Subtrac a constant from the polynomial
            ///
            /// # Example
            ///
            /// ```rust
            /// # use series::Polynomial;
            /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
            /// p -= 1;
            /// p -= &1;
            /// let res = Polynomial::new("x", -3, vec![1, 0, -3, -2]);
            /// assert_eq!(res, p);
            /// ```
            fn sub_assign(&mut self, other: $t) {
                if other.is_zero() {
                    return;
                }
                match self {
                    Polynomial::Const(c) => c.sub_assign(other),
                    Polynomial::Poly(NonConstPoly {
                        min_pow,
                        coeffs,
                        var: _,
                    }) => {
                        extend_to_range(coeffs, min_pow, 0..1);
                        let pos = (-*min_pow) as usize;
                        coeffs[pos].sub_assign(other);
                        self.trim();
                    }
                }
            }
        }
    };
}

impl_sub_assign_const!(C);
impl_sub_assign_const!(&'a C);


impl<'a, Var, C> AddAssign<&'a Polynomial<Var, C>> for Polynomial<Var, C>
where
    C: Coeff + Clone,
    Var: Clone + Debug + PartialEq,
    for<'c> C: AddAssign<&'c C>,
{
    /// Set p = p + q for two Laurent polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec!(1., 0., -3.));
    /// let q = Polynomial::new("x", -1, vec!(3., 4., 5.));
    /// let res = Polynomial::new("x", -3, vec!(1., 0., 0., 4., 5.));
    /// p += &q;
    /// assert_eq!(res, p);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the polynomials are non-constant and have different variables.
    fn add_assign(&mut self, other: &'a Polynomial<Var, C>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, Var: Clone + Debug + PartialEq, C: Coeff + Clone>
    AddAssign<PolynomialSlice<'a, Var, C>> for Polynomial<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        match other {
            PolynomialSlice::Const(c) => self.add_assign(c),
            PolynomialSlice::Poly {
                min_pow: other_min_pow,
                coeffs: other_coeffs,
                var: other_var,
            } => match self {
                Polynomial::Const(c) => {
                    let mut res = Polynomial::from(other);
                    res.add_assign(&*c);
                    *self = res;
                }
                Polynomial::Poly(NonConstPoly {
                    min_pow,
                    coeffs,
                    var,
                }) => {
                    assert_eq!(var, other_var);
                    let other_max_pow =
                        other_min_pow + other_coeffs.len() as isize;
                    let pow_range = other_min_pow..other_max_pow;
                    extend_to_range(coeffs, min_pow, pow_range);
                    for (pow, coeff) in other.iter() {
                        coeffs[(pow - *min_pow) as usize].add_assign(coeff);
                    }
                    self.trim();
                }
            },
        }
    }
}

impl<Var, C: Coeff> AddAssign for Polynomial<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
    C: AddAssign,
    Var: Debug + PartialEq,
{
    /// Set p = p + q for two Laurent polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let q = Polynomial::new("x", -1, vec![3, 4, 5]);
    /// let res = Polynomial::new("x", -3, vec![1, 0, 0, 4, 5]);
    /// p += q;
    /// assert_eq!(res, p);
    /// ```
    fn add_assign(&mut self, other: Polynomial<Var, C>) {
        match (&mut *self, other) {
            (Polynomial::Const(c), Polynomial::Const(d)) => c.add_assign(d),
            (
                Polynomial::Const(c),
                mut other @ Polynomial::Poly(NonConstPoly { .. }),
            ) => {
                other.add_assign(std::mem::replace(c, C::zero()));
                *self = other;
            }
            (Polynomial::Poly(NonConstPoly { .. }), Polynomial::Const(d)) => {
                self.add_assign(d)
            }
            (
                Polynomial::Poly(NonConstPoly {
                    min_pow,
                    coeffs,
                    var,
                }),
                Polynomial::Poly(NonConstPoly {
                    min_pow: mut other_min_pow,
                    coeffs: mut other_coeffs,
                    var: other_var,
                }),
            ) => {
                assert_eq!(*var, other_var);
                if other_coeffs.len() > coeffs.len() {
                    std::mem::swap(coeffs, &mut other_coeffs);
                    std::mem::swap(min_pow, &mut other_min_pow);
                }
                let other_max_pow = other_min_pow + other_coeffs.len() as isize;
                let pow_range = other_min_pow..other_max_pow;
                extend_to_range(coeffs, min_pow, pow_range);
                let lhs = &mut coeffs[(other_min_pow - *min_pow) as usize..];
                for (lhs, rhs) in lhs.iter_mut().zip(other_coeffs) {
                    lhs.add_assign(rhs);
                }
                self.trim();
            }
        }
    }
}

macro_rules! impl_add_via_add_assign {
    ($($t:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Add<$t> for Polynomial<Var, C>
            where
                Polynomial<Var, C>: AddAssign<$t>,
            {
                type Output = Polynomial<Var, C>;

                fn add(mut self, other: $t) -> Self::Output {
                    self.add_assign(other);
                    self
                }
            }
        )*
    };
}

impl_add_via_add_assign!(Self, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, C, &'a C);

// TODO: avoid potentially costly clone
impl<'a, Var: Clone, C: Coeff + Clone> SubAssign<PolynomialSlice<'a, Var, C>>
    for Polynomial<Var, C>
where
    Polynomial<Var, C>: SubAssign,
{
    /// Set p = p - q for two polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    fn sub_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        self.sub_assign(Polynomial::from(other));
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> SubAssign<&'a Polynomial<Var, C>>
    for Polynomial<Var, C>
where
    Polynomial<Var, C>: SubAssign<PolynomialSlice<'a, Var, C>>,
{
    /// Set p = p - q for two polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    fn sub_assign(&mut self, other: &'a Polynomial<Var, C>) {
        self.sub_assign(other.as_slice(..));
    }
}

impl<Var, C: Coeff> SubAssign for Polynomial<Var, C>
where
    Polynomial<Var, C>: AddAssign + Neg<Output = Polynomial<Var, C>>,
{
    /// Set p = p - q for two polynomial p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// let res = Polynomial::zero();
    /// p -= p.clone();
    /// assert_eq!(res, p);
    /// ```
    fn sub_assign(&mut self, other: Polynomial<Var, C>) {
        *self += -other;
    }
}

macro_rules! impl_sub_via_sub_assign {
    ($($t:ty), *) => {
        $(
            impl<'a, Var, C: Coeff> Sub<$t> for Polynomial<Var, C>
            where
                Polynomial<Var, C>: SubAssign<$t>,
            {
                type Output = Polynomial<Var, C>;

                fn sub(mut self, other: $t) -> Self::Output {
                    self.sub_assign(other);
                    self
                }
            }
        )*
    };
}

impl_sub_via_sub_assign!(Self, &'a Polynomial<Var, C>, PolynomialSlice<'a, Var, C>, C, &'a C);

impl<'a, Var, C: Coeff + Clone + AddAssign> MulAssign<&'a Polynomial<Var, C>>
    for Polynomial<Var, C>
where
    Polynomial<Var, C>: MulAssign<PolynomialSlice<'a, Var, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec![1., 0., -6., 0., 9.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a Polynomial<Var, C>) {
        self.mul_assign(other.as_slice(..))
    }
}

// TODO: pass `var` in `Mul` so it does not have to be cloned

impl<'a, Var, C: Coeff> MulAssign<PolynomialSlice<'a, Var, C>>
    for Polynomial<Var, C>
where
    for<'b> PolynomialSlice<'b, Var, C>:
        Mul<PolynomialSlice<'a, Var, C>, Output = Polynomial<Var, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    fn mul_assign(&mut self, other: PolynomialSlice<'a, Var, C>) {
        let prod = self.as_slice(..) * other;
        *self = prod;
    }
}

impl<Var, C: Coeff> MulAssign for Polynomial<Var, C>
where
    for<'a> Polynomial<Var, C>: MulAssign<&'a Polynomial<Var, C>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec![1., 0., -6. ,0. ,9.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: Polynomial<Var, C>) {
        *self *= &other
    }
}

impl<Var, C: Coeff> MulAssign<C> for Polynomial<Var, C>
where
    for<'a> C: MulAssign<&'a C>,
{
    /// Multiply by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= 2.;
    /// let res = Polynomial::new("x", -3, vec![2., 0., -6.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: C) {
        self.mul_assign(&other)
    }
}

impl<'a, Var, C: Coeff> MulAssign<&'a C> for Polynomial<Var, C>
where
    C: MulAssign<&'a C>,
{
    /// Multiply by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= &2.;
    /// let res = Polynomial::new("x", -3, vec![2., 0., -6.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a C) {
        match self {
            Polynomial::Const(c) => c.mul_assign(other),
            Polynomial::Poly(NonConstPoly { coeffs, .. }) => {
                for c in coeffs {
                    c.mul_assign(other)
                }
                // `other` might be zero or a zero divisor in some ring
                self.trim();
            }
        }
    }
}

impl<Var, C: Coeff> DivAssign<C> for Polynomial<Var, C>
where
    for<'a> C: DivAssign<&'a C>,
{
    /// Divide by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p /= 2.;
    /// let res = Polynomial::new("x", -3, vec![0.5, 0., -1.5]);
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: C) {
        self.div_assign(&other)
    }
}

impl<'a, Var, C: Coeff> DivAssign<&'a C> for Polynomial<Var, C>
where
    C: DivAssign<&'a C>,
{
    /// Divide by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p /= &2.;
    /// let res = Polynomial::new("x", -3, vec![0.5, 0., -1.5]);
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: &'a C) {
        match self {
            Polynomial::Const(c) => c.div_assign(other),
            Polynomial::Poly(NonConstPoly { coeffs, .. }) => {
                for c in coeffs {
                    c.div_assign(other)
                }
                // Division can result in zero, e.g. over integers
                self.trim();
            }
        }
    }
}

impl<Var, C: Coeff> Mul for Polynomial<Var, C>
where
    Polynomial<Var, C>: MulAssign,
{
    type Output = Polynomial<Var, C>;

    fn mul(mut self, other: Polynomial<Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a Polynomial<Var, C>> for Polynomial<Var, C>
where
    Polynomial<Var, C>: MulAssign<PolynomialSlice<'a, Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: &'a Polynomial<Var, C>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, Var, C: Coeff> Mul<PolynomialSlice<'a, Var, C>> for Polynomial<Var, C>
where
    Polynomial<Var, C>: MulAssign<PolynomialSlice<'a, Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(mut self, other: PolynomialSlice<'a, Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<Var, C: Coeff> Mul<C> for Polynomial<Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a C> for Polynomial<Var, C>
where
    C: MulAssign<&'a C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<Var, C: Coeff> Div<C> for Polynomial<Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn div(mut self, other: C) -> Self::Output {
        self /= &other;
        self
    }
}

impl<'a, Var, C: Coeff> Div<&'a C> for Polynomial<Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn div(mut self, other: &'a C) -> Self::Output {
        self /= other;
        self
    }
}

impl<'a, Var, C: Coeff, T> Mul<T> for &'a Polynomial<Var, C>
where
    PolynomialSlice<'a, Var, C>: Mul<T, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, Var, C: Coeff, T> Div<T> for &'a Polynomial<Var, C>
where
    PolynomialSlice<'a, Var, C>: Div<T, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn div(self, other: T) -> Self::Output {
        self.as_slice(..) / other
    }
}

impl<'a, Var, C: Coeff, T> Add<T> for &'a Polynomial<Var, C>
where
    PolynomialSlice<'a, Var, C>: Add<T, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, other: T) -> Self::Output {
        self.as_slice(..) + other
    }
}

impl<'a, Var, C: Coeff, T> Sub<T> for &'a Polynomial<Var, C>
where
    PolynomialSlice<'a, Var, C>: Sub<T, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, other: T) -> Self::Output {
        self.as_slice(..) - other
    }
}

impl<Var, C: Coeff> Zero for Polynomial<Var, C>
where
    Polynomial<Var, C>: Add<Output = Polynomial<Var, C>>,
{
    fn zero() -> Self {
        Polynomial::zero()
    }

    fn is_zero(&self) -> bool {
        Polynomial::is_zero(self)
    }
}

impl<Var, C: AddAssign + Coeff + Clone> One for Polynomial<Var, C>
where
    Polynomial<Var, C>: Add<Output = Polynomial<Var, C>>,
    Polynomial<Var, C>: Mul<Output = Polynomial<Var, C>>,
{
    fn one() -> Self {
        Polynomial::one()
    }

    fn is_one(&self) -> bool {
        Polynomial::is_one(self)
    }
}

/// View into a Laurent polynomial
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub enum PolynomialSlice<'a, Var, C> {
    Const(&'a C),
    Poly {
        min_pow: isize,
        coeffs: &'a [C],
        var: &'a Var,
    },
}

impl<Var, C: Coeff> std::marker::Copy for PolynomialSlice<'_, Var, C> {}

impl<Var, C: Coeff> std::clone::Clone for PolynomialSlice<'_, Var, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C: Coeff + Clone, Var: Clone> From<PolynomialSlice<'a, Var, C>>
    for Polynomial<Var, C>
{
    fn from(value: PolynomialSlice<'a, Var, C>) -> Self {
        match value {
            PolynomialSlice::Const(c) => Polynomial::Const(c.clone()),
            PolynomialSlice::Poly {
                min_pow,
                coeffs,
                var,
            } => Polynomial::new(var.clone(), min_pow, coeffs.to_owned()),
        }
    }
}

impl<'a, Var: 'a, C: Coeff + 'a> PolynomialSlice<'a, Var, C> {
    /// Get the leading power of the polynomial variable
    ///
    /// See [Polynomial::min_pow] for details.
    pub fn min_pow(self) -> Option<isize> {
        match self {
            PolynomialSlice::Const(c) => {
                if c.is_zero() {
                    None
                } else {
                    Some(0)
                }
            }
            PolynomialSlice::Poly { min_pow, .. } => Some(min_pow),
        }
    }

    /// Get the highest power of the polynomial variable
    ///
    /// See [Polynomial::max_pow] for details.
    pub fn max_pow(self) -> Option<isize> {
        self.min_pow().map(|c| c + (self.len() - 1) as isize)
    }

    /// Get the difference between the highest and the lowest power of
    /// the polynomial variable
    ///
    /// See [Polynomial::len] for details.
    pub fn len(self) -> usize {
        match self {
            PolynomialSlice::Const(c) => {
                if c.is_zero() {
                    0
                } else {
                    1
                }
            }
            PolynomialSlice::Poly { coeffs, .. } => coeffs.len(),
        }
    }

    /// Check if the polynomial is zero
    ///
    /// See [Polynomial::is_zero].
    pub fn is_empty(self) -> bool {
        match self {
            PolynomialSlice::Const(c) => c.is_zero(),
            PolynomialSlice::Poly { coeffs, .. } => coeffs.is_empty(),
        }
    }

    /// Iterator over the polynomial powers and coefficients.
    ///
    /// See [Polynomial::iter] for details.
    pub fn iter(self) -> Iter<'a, C> {
        match self {
            PolynomialSlice::Const(c) => Iter {
                pow: 0,
                coeffs: if c.is_zero() { &[] } else { slice::from_ref(c) },
            },
            PolynomialSlice::Poly {
                min_pow,
                coeffs,
                var: _,
            } => Iter {
                pow: min_pow,
                coeffs,
            },
        }
    }

    /// Try to get the coefficient of the polynomial variable to the given power
    ///
    /// See [Polynomial::get_coeff] for details.
    pub fn get_coeff(self, pow: isize) -> Option<&'a C> {
        match self {
            PolynomialSlice::Const(c) => {
                if pow != 0 || c.is_zero() {
                    None
                } else {
                    Some(c)
                }
            }
            PolynomialSlice::Poly {
                min_pow,
                coeffs,
                var: _,
            } => coeffs.get((pow - min_pow) as usize),
        }
    }

    /// Get the polynomial variable
    ///
    /// See [Polynomial::var] for details.
    pub fn var(&self) -> Option<&'a Var> {
        if let PolynomialSlice::Poly { var, .. } = self {
            Some(var)
        } else {
            None
        }
    }

    /// Check if the polynomial is constant
    pub const fn is_const(&self) -> bool {
        matches!(self, PolynomialSlice::Const(..))
    }

    pub fn new(mut min_pow: isize, mut coeffs: &'a [C], var: &'a Var) -> Self {
        min_pow += trim_slice_zero(&mut coeffs) as isize;
        Self::Poly {
            min_pow,
            coeffs,
            var,
        }
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for PolynomialSlice<'a, Var, C>
where
    &'a C: Neg,
    <&'a C as Neg>::Output: Coeff,
{
    type Output = Polynomial<Var, <&'a C as Neg>::Output>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// # use series::Polynomial;
    /// let p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let minus_p = Polynomial::new("x", -3, vec![-1, 0, 3]);
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        match self {
            PolynomialSlice::Const(c) => Polynomial::from_const(-c),
            PolynomialSlice::Poly {
                min_pow,
                coeffs,
                var,
            } => {
                let coeffs = coeffs.iter().map(Neg::neg).collect();
                Polynomial::new(var.clone(), min_pow, coeffs)
            }
        }
    }
}

impl<'a, Var: Clone, C> Add<C> for PolynomialSlice<'a, Var, C>
where
    C: Coeff + Clone + AddAssign,
{
    type Output = Polynomial<Var, C>;

    fn add(self, rhs: C) -> Self::Output {
        Polynomial::from(self).add(rhs)
    }
}

impl<'a, 'b, Var: Clone, C> Add<&'b C> for PolynomialSlice<'a, Var, C>
where
    C: Coeff + Clone + AddAssign<&'b C>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, rhs: &'b C) -> Self::Output {
        Polynomial::from(self).add(rhs)
    }
}

impl<'a, Var: Clone, C> Sub<C> for PolynomialSlice<'a, Var, C>
where
    C: Coeff + Clone + SubAssign,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, rhs: C) -> Self::Output {
        Polynomial::from(self).sub(rhs)
    }
}

impl<'a, 'b, Var: Clone, C> Sub<&'b C> for PolynomialSlice<'a, Var, C>
where
    C: Coeff + Clone + SubAssign<&'b C>,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, rhs: &'b C) -> Self::Output {
        Polynomial::from(self).sub(rhs)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> Mul<C> for PolynomialSlice<'a, Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, rhs: C) -> Self::Output {
        Polynomial::from(self).mul(rhs)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> Mul<&C> for PolynomialSlice<'a, Var, C>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, rhs: &C) -> Self::Output {
        Polynomial::from(self).mul(rhs)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> Div<C> for PolynomialSlice<'a, Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn div(self, rhs: C) -> Self::Output {
        Polynomial::from(self).div(rhs)
    }
}

impl<'a, Var: Clone, C: Coeff + Clone> Div<&C> for PolynomialSlice<'a, Var, C>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<Var, C>;

    fn div(self, rhs: &C) -> Self::Output {
        Polynomial::from(self).div(rhs)
    }
}

impl<'a, Var, C> Add for PolynomialSlice<'a, Var, C>
where
    Var: Clone + Debug + PartialEq,
    C: Coeff + Clone,
    for<'c> C: AddAssign<&'c C>,
    &'a C: Add<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (PolynomialSlice::Const(c), PolynomialSlice::Const(d)) => {
                Polynomial::Const(c + d)
            }
            (PolynomialSlice::Const(c), PolynomialSlice::Poly { .. }) => {
                Polynomial::from(rhs) + c
            }
            (PolynomialSlice::Poly { .. }, PolynomialSlice::Const(c)) => {
                Polynomial::from(self) + c
            }
            (
                PolynomialSlice::Poly { coeffs, .. },
                PolynomialSlice::Poly {
                    coeffs: rhs_coeffs, ..
                },
            ) => {
                if coeffs.len() >= rhs_coeffs.len() {
                    Polynomial::from(self) + rhs
                } else {
                    Polynomial::from(rhs) + self
                }
            }
        }
    }
}

impl<'a, Var, C: Coeff> Add<&'a Polynomial<Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Self: Add,
{
    type Output = <Self as Add>::Output;

    fn add(self, rhs: &'a Polynomial<Var, C>) -> Self::Output {
        self.add(rhs.as_slice(..))
    }
}

impl<'a, Var, C: Coeff> Add<Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    Polynomial<Var, C>: Add<Self>,
{
    type Output = <Polynomial<Var, C> as Add<Self>>::Output;

    fn add(self, rhs: Polynomial<Var, C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, Var, C> Sub for PolynomialSlice<'a, Var, C>
where
    Var: Clone + Debug + PartialEq,
    C: Coeff + Clone + Neg<Output = C> + AddAssign,
    for<'c> C: AddAssign<&'c C> + SubAssign<&'c C>,
    &'a C: Sub<Output = C>,
{
    type Output = Polynomial<Var, C>;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (PolynomialSlice::Const(c), PolynomialSlice::Const(d)) => {
                Polynomial::Const(c - d)
            }
            (PolynomialSlice::Const(c), PolynomialSlice::Poly { .. }) => {
                -(Polynomial::from(rhs) - c)
            }
            (PolynomialSlice::Poly { .. }, PolynomialSlice::Const(c)) => {
                Polynomial::from(self) - c
            }
            (
                PolynomialSlice::Poly { coeffs, .. },
                PolynomialSlice::Poly {
                    coeffs: rhs_coeffs, ..
                },
            ) => {
                if coeffs.len() >= rhs_coeffs.len() {
                    Polynomial::from(self) - rhs
                } else {
                    -(Polynomial::from(rhs) - self)
                }
            }
        }
    }
}

impl<'a, Var, C: Coeff> Sub<&'a Polynomial<Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Self: Sub,
{
    type Output = <Self as Sub>::Output;

    fn sub(self, rhs: &'a Polynomial<Var, C>) -> Self::Output {
        self.sub(rhs.as_slice(..))
    }
}

impl<'a, Var, C: Coeff> Sub<Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    Polynomial<Var, C>: Sub<Self, Output = Polynomial<Var, C>>
        + Neg<Output = Polynomial<Var, C>>,
{
    type Output = <Polynomial<Var, C> as Sub<Self>>::Output;

    fn sub(self, rhs: Polynomial<Var, C>) -> Self::Output {
        -rhs.sub(self)
    }
}

impl<'a, 'b, Var, C: Coeff> Mul<&'b Polynomial<Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    Self: Mul<PolynomialSlice<'b, Var, C>, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, rhs: &'b Polynomial<Var, C>) -> Self::Output {
        self.mul(rhs.as_slice(..))
    }
}

impl<'a, Var, C: Coeff> Mul<Polynomial<Var, C>> for PolynomialSlice<'a, Var, C>
where
    for<'c> Self: Mul<PolynomialSlice<'c, Var, C>, Output = Polynomial<Var, C>>,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, rhs: Polynomial<Var, C>) -> Self::Output {
        self.mul(rhs.as_slice(..))
    }
}

impl<'a, 'b, Var, C> Mul<PolynomialSlice<'b, Var, C>>
    for PolynomialSlice<'a, Var, C>
where
    C: Coeff + Clone + AddAssign,
    for<'c> &'c C: Mul<Output = C>,
    for<'c> C: MulAssign<&'c C>,
    Var: Debug + Clone + PartialEq,
{
    type Output = Polynomial<Var, C>;

    fn mul(self, rhs: PolynomialSlice<'b, Var, C>) -> Self::Output {
        match (self, rhs) {
            (PolynomialSlice::Const(c), PolynomialSlice::Const(d)) => {
                Polynomial::Const(c * d)
            }
            (PolynomialSlice::Const(c), PolynomialSlice::Poly { .. }) => {
                Polynomial::from(rhs) * c
            }
            (PolynomialSlice::Poly { .. }, PolynomialSlice::Const(c)) => {
                Polynomial::from(self) * c
            }
            (
                PolynomialSlice::Poly {
                    min_pow,
                    coeffs,
                    var,
                },
                PolynomialSlice::Poly {
                    min_pow: other_min_pow,
                    coeffs: other_coeffs,
                    var: other_var,
                },
            ) => {
                assert_eq!(var, other_var);
                let res_min_pow = min_pow + other_min_pow;
                let res_len = coeffs.len() + other_coeffs.len();
                let mut res_coeffs = Vec::with_capacity(res_len);
                for n in 0..res_len {
                    let mut c = C::zero();
                    let imin = 1 + n - std::cmp::min(other_coeffs.len(), 1 + n);
                    let imax = std::cmp::min(n + 1, coeffs.len());
                    for i in imin..imax {
                        c += &coeffs[i] * &other_coeffs[n - i];
                    }
                    res_coeffs.push(c)
                }
                Polynomial::new(var.clone(), res_min_pow, res_coeffs)
            }
        }
    }
}

impl<'a, Var, C: 'static + Coeff + Send + Sync> PolynomialSlice<'a, Var, C> {
    pub fn zero() -> Self {
        Self::Const(zero_ref())
    }

    pub fn coeff(self, pow: isize) -> &'a C {
        self.get_coeff(pow).unwrap_or(zero_ref())
    }
}

macro_rules! impl_num_display {
    ($($t:ty), *) => {
        $(
            impl<'a, Var: Display> Display for PolynomialSlice<'a, Var, $t> {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    match self {
                        PolynomialSlice::Const(c) => return write!(f, "{c}"),
                        PolynomialSlice::Poly { min_pow, coeffs, var } => {
                            if coeffs.is_empty() {
                                return write!(f, "0");
                            }
                            let terms = coeffs.iter()
                                .enumerate()
                                .filter_map(|(n, c)| if c.is_zero() {
                                    None
                                } else {
                                    Some((*min_pow + n as isize, *c))
                                });
                            let mut first = true;
                            for (pow, mut c) in terms {
                                if !first {
                                    if c.starts_with_minus() {
                                        c = c.abs();
                                        write!(f, " - ")?;
                                    } else {
                                        write!(f, " + ")?;
                                    }
                                }
                                first = false;
                                if pow == 0 {
                                    write!(f, "{c}")?;
                                } else {
                                    if !c.is_one() {
                                        write!(f, "{c}*")?;
                                    }
                                    write!(f, "{var}")?;
                                    if pow != 1 {
                                        write!(f, "^{pow}")?;
                                    }
                                }
                            }
                            Ok(())
                        },
                    }
                }
            }
        )*
    };
}

impl_num_display!(
    i8, i16, i32, i64, i128, isize, f32, f64, u8, u16, u32, u64, u128, usize
);

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Iter<'a, C> {
    pow: isize,
    coeffs: &'a [C],
}

impl<C: Coeff> std::marker::Copy for Iter<'_, C> {}

impl<C: Coeff> std::clone::Clone for Iter<'_, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C> ExactSizeIterator for Iter<'a, C> {}

impl<'a, C> DoubleEndedIterator for Iter<'a, C> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let (coeff, rest) = self.coeffs.split_last()?;
        self.coeffs = rest;
        Some((self.pow + self.coeffs.len() as isize, coeff))
    }
}

impl<'a, C> FusedIterator for Iter<'a, C> {}

impl<'a, C> Iterator for Iter<'a, C> {
    type Item = (isize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        let (coeff, rest) = self.coeffs.split_first()?;
        self.coeffs = rest;
        let pow = self.pow;
        self.pow += 1;
        Some((pow, coeff))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.coeffs.len();
        (len, Some(len))
    }

    fn count(self) -> usize {
        self.len()
    }

    fn last(self) -> Option<Self::Item> {
        let Self { pow, coeffs: coeff } = self;
        let last = coeff.last()?;
        Some((pow + coeff.len() as isize, last))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.pow += n as isize;
        self.coeffs = &self.coeffs[n..];
        self.next()
    }
}
