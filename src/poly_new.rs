use crate::traits::AsSlice;
use crate::util_new::trim_zero;
use crate::zero_ref::zero_ref;
use crate::{Coeff, IntoIter};

use core::slice;
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
pub enum Polynomial<C, V> {
    Const(C),
    Poly(NonConstPoly<C, V>),
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
/// A non-constant polynomial
pub struct NonConstPoly<C, V> {
    min_pow: isize,
    coeffs: Vec<C>,
    var: V,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
/// Data parts of a polynomial
///
/// # Example
///
/// ```rust
/// // destructure a polynomial
/// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
/// let Polynomial::Poly(p) = p else {
///    unreachable!("Polynomial is not a constant")
/// };
/// let PolynomialParts{min_pow, coeffs, var} = p.into();
/// assert_eq!(min_pow, Some(-1));
/// assert_eq!(coeffs, vec![1, 2, 3]);
/// assert_eq!(var, "x");
/// ```
pub struct PolynomialParts<C, V> {
    pub min_pow: isize,
    pub coeffs: Vec<C>,
    pub var: V,
}

impl<C, V> From<NonConstPoly<C, V>> for PolynomialParts<C, V> {
    fn from(value: NonConstPoly<C, V>) -> Self {
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

impl<C: Coeff, V> Polynomial<C, V> {
    /// Create a new Laurent polynomial
    ///
    /// # Example
    ///
    /// This creates a Laurent polynomial in x, starting at power -1 with
    /// coefficients 1, 2, 3. In other words, the polynomial x^-1 + 2 + 3*x.
    /// ```rust
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// ```
    pub fn new(var: V, min_pow: isize, coeffs: Vec<C>) -> Polynomial<C, V> {
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
    /// let two: Polynomial<_, ()> = Polynomial::from_const(2);
    /// ```
    pub const fn from_const(coeff: C) -> Polynomial<C, V> {
        Self::Const(coeff)
    }

    /// Create the zero polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// let zero: Polynomial<_, ()> = Polynomial::zero();
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
    /// type IntPoly = Polynomial<i32, ()>;
    /// assert!(IntPoly::zero.is_zero());
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
    /// let one: Polynomial<_, ()> = Polynomial::one();
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
    /// type IntPoly = Polynomial<i32, ()>;
    /// assert!(IntPoly::one.is_one());
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
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.min_pow(), Some(-1));
    ///
    /// type IntPoly = Polynomial<i32, ()>;
    /// let p = IntPoly::from_coeff(1);
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
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.max_pow(), Some(1));
    ///
    /// let p: Polynomial<i32, ()> = Polynomial::zero();
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
    /// let p = Polynomial::new("x", -1, vec![1 ,2, 3]);
    /// assert_eq!(p.len(), 3);
    ///
    /// let p: Polynomial<i32, ()> = Polynomial::zero();
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
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let mut iter = p.iter();
    /// assert_eq!(iter.next(), Some((-1, &1)));
    /// assert_eq!(iter.next(), Some((0, &2)));
    /// assert_eq!(iter.next(), Some((1, &3)));
    /// assert_eq!(iter.next(), None);
    ///
    /// let p: Polynomial<i32, ()> = Polynomial::zero();
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
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3, 4]);
    /// let p = p.map(|_, c| c * c);
    /// assert_eq!(p.coeff(-1), &1);
    /// assert_eq!(p.coeff(0), &4);
    /// assert_eq!(p.coeff(1), &9);
    /// assert_eq!(p.coeff(2), &16);
    /// ```
    pub fn map<D, F>(self, mut f: F) -> Polynomial<D, V>
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
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// assert_eq!(p.var(), Some(&"x"));
    /// let p = Polynomial::from_const(2);
    /// assert!(p.var().is_none());
    /// ```
    pub fn var(&self) -> Option<&V> {
        self.as_slice(..).var()
    }

    /// Replace the polynomial variable
    ///
    /// Returns the new polynomial and, if the polynomial was not a
    /// constant, the old variable.
    /// # Example
    ///
    /// ```rust
    /// let p = Polynomial::new("x", -1, vec![1, 2, 3]);
    /// let (p, var) = p.replace_var("y");
    /// assert_eq!(p.var(), Some(&"y"));
    /// assert_eq!(var, Some("x"));
    ///
    /// let p = Polynomial::from_const(2);
    /// let (p, var) = p.replace_var("y");
    /// assert!(var.is_none());
    /// assert!(p.var().is_none());
    /// ```
    pub fn replace_var<W>(self, new_var: W) -> (Polynomial<C, W>, Option<V>) {
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

//     /// Turn a polynomial into a series with the given cutoff
//     ///
//     /// # Example
//     ///
//     /// ```rust
//     /// let p = Polynomial::new(-1, vec!(1,2,3));
//     /// let s = with_cutoff(-1..5, vec!(1,2,3));
//     /// assert_eq!(p.cutoff_at(5), s);
//     /// ```
//     ///
//     /// # Panics
//     ///
//     /// Panics if the cutoff power is lower than the starting power
//     ///
//     pub fn cutoff_at(self, cutoff_pow: isize) -> Series<C> {
//         with_cutoff(
//             self.min_pow.unwrap_or(cutoff_pow)..cutoff_pow,
//             self.coeffs,
//         )
//     }

impl<C: 'static + Coeff + Send + Sync, V> Polynomial<C, V> {
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

impl<C: Coeff, V> Default for Polynomial<C, V> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<'a, C: 'static + Coeff + Send + Sync, V: 'a> AsSlice<'a, Range<isize>>
    for Polynomial<C, V>
{
    type Output = PolynomialSlice<'a, C, V>;

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
                PolynomialSlice::Poly {
                    min_pow: *min_pow,
                    coeffs: &coeffs[start..end],
                    var,
                }
            }
        }
    }
}

impl<'a, C: 'static + Coeff + Send + Sync, V: 'a>
    AsSlice<'a, RangeInclusive<isize>> for Polynomial<C, V>
{
    type Output = PolynomialSlice<'a, C, V>;

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
                PolynomialSlice::Poly {
                    min_pow: *min_pow,
                    coeffs: &coeffs[start..=end],
                    var,
                }
            }
        }
    }
}

impl<'a, C: 'a + Coeff, V: 'a> AsSlice<'a, RangeToInclusive<isize>>
    for Polynomial<C, V>
{
    type Output = PolynomialSlice<'a, C, V>;

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
                PolynomialSlice::Poly {
                    min_pow: *min_pow,
                    coeffs: &coeffs[r],
                    var,
                }
            }
        }
    }
}

impl<'a, C: 'a + Coeff, V: 'a> AsSlice<'a, RangeFrom<isize>>
    for Polynomial<C, V>
{
    type Output = PolynomialSlice<'a, C, V>;

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
                PolynomialSlice::Poly {
                    min_pow: *min_pow,
                    coeffs: &coeffs[r],
                    var,
                }
            }
        }
    }
}

impl<'a, C: 'a + Coeff, V: 'a> AsSlice<'a, RangeTo<isize>>
    for Polynomial<C, V>
{
    type Output = PolynomialSlice<'a, C, V>;

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
                PolynomialSlice::Poly {
                    min_pow: *min_pow,
                    coeffs: &coeffs[r],
                    var,
                }
            }
        }
    }
}

impl<'a, C: 'a + Coeff, V: 'a> AsSlice<'a, RangeFull> for Polynomial<C, V> {
    type Output = PolynomialSlice<'a, C, V>;

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

// impl<C: Coeff> convert::From<Series<C>> for Polynomial<C> {
//     fn from(s: Series<C>) -> Self {
//         Polynomial::new(s.min_pow, s.coeffs)
//     }
// }

impl<C: Coeff, V> Index<isize> for Polynomial<C, V> {
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

impl<C: Coeff, V> std::iter::IntoIterator for Polynomial<C, V> {
    type Item = (isize, C);
    type IntoIter = crate::IntoIter<C>;

    /// Consuming iterator over the polynomial powers and coefficients.
    ///
    /// # Example
    ///
    /// ```rust
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
    debug_assert!(max_pow >= pow_range.end);
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

impl<C: Coeff + Neg, V> Neg for Polynomial<C, V>
where
    <C as Neg>::Output: Coeff,
{
    type Output = Polynomial<<C as Neg>::Output, V>;

    /// Compute -p for a Laurent polynomial p
    ///
    /// # Example
    ///
    /// ```rust
    /// let p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let minus_p = Polynomial::new("x", -3, vec![-1, 0, 3]);
    /// assert_eq!(-p, minus_p);
    /// ```
    fn neg(self) -> Self::Output {
        self.map(|_, c| -c)
    }
}

impl<'a, C: AddAssign + Coeff, V> AddAssign<C> for Polynomial<C, V> {
    /// Add a constant to the polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// p += 1;
    /// let res = Polynomial::new("x", -3, vec![1, 0, -3, 1]);
    /// assert_eq!(res, p);
    /// ```
    fn add_assign(&mut self, other: C) {
        if other.is_zero() {
            return;
        }
        match self {
            Polynomial::Const(c) => *c += other,
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var: _,
            }) => {
                extend_to_range(coeffs, min_pow, 0..1);
                let pos = (-*min_pow) as usize;
                coeffs[pos] += other;
                self.trim();
            }
        }
    }
}

// TODO: code duplication with AddAssign<C>
impl<'a, C: Coeff, V> AddAssign<&'a C> for Polynomial<C, V>
where
    C: AddAssign<&'a C>,
{
    /// Add a constant to the polynomial
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// p += &1;
    /// let res = Polynomial::new("x", -3, vec![1, 0, -3, 1]);
    /// assert_eq!(res, p);
    /// ```
    fn add_assign(&mut self, other: &'a C) {
        if other.is_zero() {
            return;
        }
        match self {
            Polynomial::Const(c) => *c += other,
            Polynomial::Poly(NonConstPoly {
                min_pow,
                coeffs,
                var: _,
            }) => {
                extend_to_range(coeffs, min_pow, 0..1);
                let pos = (-*min_pow) as usize;
                coeffs[pos] += other;
                self.trim();
            }
        }
    }
}

impl<'a, C, V> AddAssign<&'a Polynomial<C, V>> for Polynomial<C, V>
where
    C: Coeff + Clone,
    V: Clone + Debug + PartialEq,
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
    fn add_assign(&mut self, other: &'a Polynomial<C, V>) {
        self.add_assign(other.as_slice(..))
    }
}

impl<'a, C: Coeff + Clone, V: Clone + Debug + PartialEq>
    AddAssign<PolynomialSlice<'a, C, V>> for Polynomial<C, V>
where
    for<'c> C: AddAssign<&'c C>,
{
    fn add_assign(&mut self, other: PolynomialSlice<'a, C, V>) {
        match other {
            PolynomialSlice::Const(c) => *self += c,
            PolynomialSlice::Poly {
                min_pow: other_min_pow,
                coeffs: other_coeffs,
                var: other_var,
            } => match self {
                Polynomial::Const(c) => {
                    let mut res = Polynomial::from(other);
                    res += &*c;
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
                        coeffs[(pow - *min_pow) as usize] += coeff;
                    }
                    self.trim();
                }
            },
        }
    }
}

impl<C: Coeff, V> AddAssign<Polynomial<C, V>> for Polynomial<C, V>
where
    for<'c> C: AddAssign<&'c C>,
    C: Clone + AddAssign,
    V: Debug + PartialEq,
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
    /// let mut p = Polynomial::new("x", -3, vec![1, 0, -3]);
    /// let q = Polynomial::new("x", -1, vec![3, 4, 5]);
    /// let res = Polynomial::new("x", -3, vec![1, 0, 0, 4, 5]);
    /// p += q;
    /// assert_eq!(res, p);
    /// ```
    fn add_assign(&mut self, other: Polynomial<C, V>) {
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

impl<C: Coeff + Clone, V, Rhs> Add<Rhs> for Polynomial<C, V>
where
    Polynomial<C, V>: AddAssign<Rhs>,
{
    type Output = Polynomial<C, V>;

    /// Add two Laurent polynomials
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

// TODO: avoid potentially costly clone
impl<'a, C: Coeff, V> SubAssign<&'a Polynomial<C, V>> for Polynomial<C, V>
where
    Polynomial<C, V>: Clone + SubAssign,
{
    /// Set p = p - q for two polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// let res = Polynomial::zero();
    /// p -= &p.clone();
    /// assert_eq!(res, p);
    /// ```
    fn sub_assign(&mut self, other: &'a Polynomial<C, V>) {
        *self -= other.to_owned();
    }
}

// TODO: avoid potentially costly clone
impl<'a, C: Coeff + Clone, V: Clone> SubAssign<PolynomialSlice<'a, C, V>>
    for Polynomial<C, V>
where
    Polynomial<C, V>: SubAssign,
{
    /// Set p = p - q for two polynomials p and q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    fn sub_assign(&mut self, other: PolynomialSlice<'a, C, V>) {
        *self -= Polynomial::from(other);
    }
}

impl<C: Coeff, V> SubAssign<Polynomial<C, V>> for Polynomial<C, V>
where
    Polynomial<C, V>: AddAssign + Neg<Output = Polynomial<C, V>>,
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
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// let res = Polynomial::zero();
    /// p -= p.clone();
    /// assert_eq!(res, p);
    /// ```
    fn sub_assign(&mut self, other: Polynomial<C, V>) {
        *self += -other;
    }
}

impl<C: Coeff, V, T> Sub<T> for Polynomial<C, V>
where
    Polynomial<C, V>: SubAssign<T>,
{
    type Output = Polynomial<C, V>;

    /// Subtract two Laurent polynomials
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<'a, C: Coeff + Clone + AddAssign, V> MulAssign<&'a Polynomial<C, V>>
    for Polynomial<C, V>
where
    Polynomial<C, V>: MulAssign<PolynomialSlice<'a, C, V>>,
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
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec![1., 0., -6., 0., 9.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: &'a Polynomial<C, V>) {
        self.mul_assign(other.as_slice(..))
    }
}

// TODO: pass `var` in `Mul` so it does not have to be cloned

impl<'a, C: Coeff, V> MulAssign<PolynomialSlice<'a, C, V>> for Polynomial<C, V>
where
    for<'b> PolynomialSlice<'b, C, V>:
        Mul<PolynomialSlice<'a, C, V>, Output = Polynomial<C, V>>,
{
    /// Set p = p * q for two polynomials p,q
    ///
    /// # Panics
    ///
    /// Panics if the polynomial variables differ and neither of the
    /// polynomials is a constant.
    ///
    fn mul_assign(&mut self, other: PolynomialSlice<'a, C, V>) {
        let prod = self.as_slice(..) * other;
        *self = prod;
    }
}

impl<C: Coeff, V> MulAssign for Polynomial<C, V>
where
    for<'a> Polynomial<C, V>: MulAssign<&'a Polynomial<C, V>>,
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
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= &p.clone();
    /// let res = Polynomial::new("x", -6, vec![1., 0., -6. ,0. ,9.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: Polynomial<C, V>) {
        *self *= &other
    }
}

impl<C: Coeff, V> MulAssign<C> for Polynomial<C, V>
where
    for<'a> C: MulAssign<&'a C>,
{
    /// Multiply by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p *= 2.;
    /// let res = Polynomial::new("x", -3, vec![2., 0., -6.]);
    /// assert_eq!(res, p);
    /// ```
    fn mul_assign(&mut self, other: C) {
        self.mul_assign(&other)
    }
}

impl<'a, C: Coeff, V> MulAssign<&'a C> for Polynomial<C, V>
where
    C: MulAssign<&'a C>,
{
    /// Multiply by a constant
    ///
    /// # Example
    ///
    /// ```rust
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

impl<C: Coeff, V> DivAssign<C> for Polynomial<C, V>
where
    for<'a> C: DivAssign<&'a C>,
{
    /// Divide by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// use Polynomial;
    /// let mut p = Polynomial::new("x", -3, vec![1., 0., -3.]);
    /// p /= 2.;
    /// let res = Polynomial::new("x", -3, vec![0.5, 0., -1.5]);
    /// assert_eq!(res, p);
    /// ```
    fn div_assign(&mut self, other: C) {
        self.div_assign(&other)
    }
}

impl<'a, C: Coeff, V> DivAssign<&'a C> for Polynomial<C, V>
where
    C: DivAssign<&'a C>,
{
    /// Divide by a constant
    ///
    /// # Example
    ///
    /// ```rust
    /// use Polynomial;
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

impl<C: Coeff, V> Mul for Polynomial<C, V>
where
    Polynomial<C, V>: MulAssign,
{
    type Output = Polynomial<C, V>;

    fn mul(mut self, other: Polynomial<C, V>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, C: Coeff, V> Mul<&'a Polynomial<C, V>> for Polynomial<C, V>
where
    Polynomial<C, V>: MulAssign<PolynomialSlice<'a, C, V>>,
{
    type Output = Polynomial<C, V>;

    fn mul(self, other: &'a Polynomial<C, V>) -> Self::Output {
        self * other.as_slice(..)
    }
}

impl<'a, C: Coeff, V> Mul<PolynomialSlice<'a, C, V>> for Polynomial<C, V>
where
    Polynomial<C, V>: MulAssign<PolynomialSlice<'a, C, V>>,
{
    type Output = Polynomial<C, V>;

    fn mul(mut self, other: PolynomialSlice<'a, C, V>) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff, V> Mul<C> for Polynomial<C, V>
where
    for<'c> C: MulAssign<&'c C>,
{
    type Output = Polynomial<C, V>;

    fn mul(mut self, other: C) -> Self::Output {
        self *= &other;
        self
    }
}

impl<'a, C: Coeff, V> Mul<&'a C> for Polynomial<C, V>
where
    C: MulAssign<&'a C>,
{
    type Output = Polynomial<C, V>;

    fn mul(mut self, other: &'a C) -> Self::Output {
        self *= other;
        self
    }
}

impl<C: Coeff, V> Div<C> for Polynomial<C, V>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<C, V>;

    fn div(mut self, other: C) -> Self::Output {
        self /= &other;
        self
    }
}

impl<'a, C: Coeff, V> Div<&'a C> for Polynomial<C, V>
where
    for<'c> C: DivAssign<&'c C>,
{
    type Output = Polynomial<C, V>;

    fn div(mut self, other: &'a C) -> Self::Output {
        self /= other;
        self
    }
}

impl<'a, C: Coeff, V, T> Mul<T> for &'a Polynomial<C, V>
where
    PolynomialSlice<'a, C, V>: Mul<T, Output = Polynomial<C, V>>,
{
    type Output = Polynomial<C, V>;

    fn mul(self, other: T) -> Self::Output {
        self.as_slice(..) * other
    }
}

impl<'a, C: Coeff, V, T> Div<T> for &'a Polynomial<C, V>
where
    PolynomialSlice<'a, C, V>: Div<T, Output = Polynomial<C, V>>,
{
    type Output = Polynomial<C, V>;

    fn div(self, other: T) -> Self::Output {
        self.as_slice(..) / other
    }
}

impl<'a, C: Coeff, V, T> Add<T> for &'a Polynomial<C, V>
where
    PolynomialSlice<'a, C, V>: Add<T, Output = Polynomial<C, V>>,
{
    type Output = Polynomial<C, V>;

    fn add(self, other: T) -> Self::Output {
        self.as_slice(..) + other
    }
}

impl<'a, C: Coeff, V, T> Sub<T> for &'a Polynomial<C, V>
where
    PolynomialSlice<'a, C, V>: Sub<T, Output = Polynomial<C, V>>,
{
    type Output = Polynomial<C, V>;

    fn sub(self, other: T) -> Self::Output {
        self.as_slice(..) - other
    }
}

impl<C: Coeff, V> Zero for Polynomial<C, V>
where
    Polynomial<C, V>: Add<Output = Polynomial<C, V>>,
{
    fn zero() -> Self {
        Polynomial::zero()
    }

    fn is_zero(&self) -> bool {
        Polynomial::is_zero(&self)
    }
}

impl<C: AddAssign + Coeff + Clone, V> One for Polynomial<C, V>
where
    Polynomial<C, V>: Add<Output = Polynomial<C, V>>,
    Polynomial<C, V>: Mul<Output = Polynomial<C, V>>,
{
    fn one() -> Self {
        Polynomial::one()
    }

    fn is_one(&self) -> bool {
        Polynomial::is_one(&self)
    }
}

/// View into a Laurent polynomial
#[derive(PartialEq, Eq, Debug, Hash, Ord, PartialOrd)]
pub enum PolynomialSlice<'a, C, V> {
    Const(&'a C),
    Poly {
        min_pow: isize,
        coeffs: &'a [C],
        var: &'a V,
    },
}

impl<C: Coeff, V> std::marker::Copy for PolynomialSlice<'_, C, V> {}

impl<C: Coeff, V> std::clone::Clone for PolynomialSlice<'_, C, V> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C: Coeff + Clone, V: Clone> From<PolynomialSlice<'a, C, V>>
    for Polynomial<C, V>
{
    fn from(value: PolynomialSlice<'a, C, V>) -> Self {
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

impl<'a, C: Coeff + 'a, V: 'a> PolynomialSlice<'a, C, V> {
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
    pub fn var(&self) -> Option<&'a V> {
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
}

impl<'a, C, V> Mul for PolynomialSlice<'a, C, V>
where
    C: Coeff + Clone + MulAssign<&'a C> + AddAssign,
    &'a C: Mul<Output = C>,
    V: Debug + Clone + PartialEq,
{
    type Output = Polynomial<C, V>;

    fn mul(self, rhs: Self) -> Self::Output {
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
                    for i in 0..=n {
                        c += &coeffs[i] * &other_coeffs[n - i];
                    }
                    res_coeffs.push(c)
                }
                Polynomial::new(var.clone(), res_min_pow, res_coeffs)
            }
        }
    }
}

impl<'a, C: 'static + Coeff + Send + Sync, V> PolynomialSlice<'a, C, V> {
    pub fn zero() -> Self {
        Self::Const(zero_ref())
    }

    pub fn coeff(self, pow: isize) -> &'a C {
        self.get_coeff(pow).unwrap_or(zero_ref())
    }
}

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
