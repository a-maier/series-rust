/// Laurent series with support for usual arithmetic operations and some
/// common functions
#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
use std::fmt;
use std::ops::{Add,AddAssign,Sub,SubAssign,Mul,MulAssign,Div,DivAssign,Neg};
use std::convert::From;
use std::cmp::min;

pub mod ops;
use self::ops::{Ln,Exp,Pow};

/// Minimum requirements on series coefficients
pub trait Coeff: From<i32> + PartialEq {}
impl<T: From<i32> + PartialEq> Coeff for T {}

/// Laurent series in a single variable up to some power
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq,Eq,Debug,Clone,Hash,Ord,PartialOrd)]
pub struct Series<Var, C: Coeff> {
    var: Var,
    min_pow: isize,
    coeffs: Vec<C>,
    zero: C // TODO: evil hack, learn how to do this better
}

impl<Var, C: Coeff> Series<Var, C> {
    /// Create a new series
    ///
    /// # Example
    ///
    /// This creates a series in the variable "x", starting at "x"^-1
    /// with coefficients 1, 2, 3. In other words, the series x^-1 + 2 +
    /// 3*x + O(x^2).
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// ```
    pub fn new(var: Var, min_pow: isize, coeffs: Vec<C>) -> Series<Var, C> {
        let mut res = Series{var, min_pow, coeffs, zero: C::from(0)};
        res.trim();
        res
    }

    /// Get the expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.var(), &"x");
    /// ```
    pub fn var(&self) -> &Var {
        &self.var
    }

    /// Get the series coefficient of the expansion variable to the
    /// given power.
    ///
    /// Returns None if the requested power is above the highest known
    /// power. Coefficients below the leading power are zero.
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.coeff(-5), Some(&0));
    /// assert_eq!(s.coeff(-2), Some(&0));
    /// assert_eq!(s.coeff(-1), Some(&1));
    /// assert_eq!(s.coeff(0), Some(&2));
    /// assert_eq!(s.coeff(1), Some(&3));
    /// assert_eq!(s.coeff(2), None);
    /// assert_eq!(s.coeff(5), None);
    /// ```
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        if pow < self.min_pow() {
            return Some(&self.zero) // TODO this is a bad hack
        }
        if pow >= self.max_pow() {
            return None
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(&self.coeffs[idx])
    }

    /// Get the leading power of the series expansion variable
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.min_pow(), -1);
    /// ```
    pub fn min_pow(&self) -> isize {
        self.min_pow
    }

    /// Get the power of the expansion variable where the series is
    /// truncated, that is \\( N \\) in a series of the form \\(
    /// \sum_{n=n_0}^{N-1} a_n x^n + O(x^N), \\)
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -1, vec!(1,2,3));
    /// assert_eq!(s.max_pow(), 2);
    /// ```
    pub fn max_pow(&self) -> isize {
        self.min_pow + (self.coeffs.len() as isize)
    }
}

/// Multiplicative inverse
pub trait MulInverse {
    type Output;

    fn mul_inverse(self) -> Self::Output;
}

impl<'a, Var: Clone, C: Coeff + SubAssign> MulInverse for &'a Series<Var, C>
where
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>
{
    type Output = Series<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = (&s).mul_inverse();
    /// let one = series::Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        let inv_min_pow = -self.min_pow;
        if self.coeffs.is_empty() {
            return Series::new(self.var.clone(), inv_min_pow, vec!())
        }
        let a : Vec<_> = self.coeffs.iter().map(|c| c/&self.coeffs[0]).collect();
        let mut b = Vec::with_capacity(a.len());
        b.push(C::from(1));
        for n in 1..a.len() {
            let mut b_n = C::from(0);
            for i in 0..n {
                b_n -= &a[n-i] * &b[i];
            }
            b.push(b_n);
        }
        let inv_coeffs : Vec<_> = b.iter().map(|b| b/&self.coeffs[0]).collect();
        Series::new(self.var.clone(), inv_min_pow, inv_coeffs)
    }
}

impl<Var: Clone, C: Coeff + SubAssign> MulInverse for Series<Var, C>
where
    for<'a> &'a Series<Var, C>: MulInverse<Output=Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Compute 1/s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::MulInverse;
    /// let s = series::Series::new("x", -1, vec!(1.,2.,3.));
    /// let s_inv = s.clone().mul_inverse();
    /// let one = series::Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(s * s_inv, one);
    /// ```
    fn mul_inverse(self) -> Self::Output {
        (&self).mul_inverse()
    }
}

impl<Var, C: Coeff> Series<Var, C> {
    fn trim(& mut self) {
        let idx = {
            let non_zero = self.coeffs.iter().enumerate().
                find(|&(_, c)| *c != C::from(0));
            match non_zero {
                Some((idx, _)) => Some(idx),
                _ => None
            }
        };
        match idx {
            Some(idx) => {
                if idx > 0 {
                    self.min_pow += idx as isize;
                    self.coeffs.drain(0..idx);
                }
            },
            None => {
                self.min_pow += self.coeffs.len() as isize;
                self.coeffs = vec!();
            }
        }
    }
}

impl<Var: fmt::Display, C: Coeff + fmt::Display> fmt::Display for Series<Var, C>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if *coeff == C::from(0) { continue; }
            let cur_pow = self.min_pow() + i as isize;
            if i > 0 { write!(f, " + ")?; }
            write!(f, "({})", coeff)?;
            if cur_pow != 0 {
                write!(f, "*{}", self.var)?;
                if cur_pow != 1 { write!(f, "^{}", cur_pow)?; }
            }
        }
        if ! self.coeffs.is_empty() { write!(f, " + ")?; }
        write!(f, "O({}^{})", self.var, self.max_pow())
    }
}

impl<Var, C: Coeff + Neg<Output=C>> Neg for Series<Var, C> {
    type Output = Series<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.into_iter().map(|c| -c).collect();
        Series::new(self.var, self.min_pow, neg_coeff)
    }
}

impl<'a, Var: Clone, C: Coeff> Neg for &'a Series<Var, C>
    where for<'c> &'c C: Neg<Output=C>
{
    type Output = Series<Var, C>;

    /// Compute -s for a series s
    ///
    /// # Example
    ///
    /// ```rust
    /// let s = series::Series::new("x", -3, vec!(1.,0.,-3.));
    /// let minus_s = series::Series::new("x", -3, vec!(-1.,0.,3.));
    /// assert_eq!(-&s, minus_s);
    /// ```
    fn neg(self) -> Self::Output {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Series::new(self.var.clone(), self.min_pow, neg_coeff)
    }
}

// Spurious trait, needed for rust 1.28
// direct implementation used to work in 1.24
trait AddAssignHelper {
    fn truncate_max_pow(&mut self, other: &Self);
    fn add_overlap(&mut self, other: &Self);
    fn num_leading(&mut self, other: &Self) -> usize;
}

impl<Var, C: Coeff> AddAssignHelper for Series<Var, C>
where for<'c> C: AddAssign<&'c C>
{
    fn truncate_max_pow(&mut self, other: &Self){
        if other.max_pow() < self.max_pow() {
            let to_remove = min(
                (self.max_pow() - other.max_pow()) as usize,
                self.coeffs.len()
            );
            let new_size = self.coeffs.len() - to_remove;
            self.coeffs.truncate(new_size);
            debug_assert!(
                self.coeffs.is_empty() || other.max_pow() == self.max_pow()
            );
        }
    }

    fn add_overlap(&mut self, other: &Self){
        let offset = self.min_pow();
        for (i, c) in self.coeffs.iter_mut().enumerate() {
            let power = offset + i as isize;
            *c += other.coeff(power).unwrap();
        }
    }

    fn num_leading(&mut self, other: &Self) -> usize {
        min(
            (self.min_pow() - other.min_pow()) as usize,
            other.coeffs.len()
        )
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone>
    AddAssign<&'a Series<Var, C>>
    for Series<Var, C>
    where for<'c> C: AddAssign<&'c C>
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = Series::new("x", -1, vec!(3., 4., 5.));
    /// let res = Series::new("x", -3, vec!(1.,0.,0.));
    /// s += &t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(& mut self, other: &'a Series<Var, C>) {
        assert_eq!(self.var, other.var);
        self.truncate_max_pow(other);
        self.add_overlap(other);
        if other.min_pow() < self.min_pow() {
            let num_leading = self.num_leading(other);
            let leading_coeff = other.coeffs[0..num_leading].iter().cloned();
            self.coeffs.splice(0..0, leading_coeff);
            self.min_pow = other.min_pow;
        }
        debug_assert!(other.max_pow() >= self.max_pow());
        self.trim();
    }
}

impl<Var, C: Coeff>
    AddAssign<Series<Var, C>> for Series<Var, C>
where
    for<'c> C: AddAssign<&'c C>,
    Var: PartialEq + fmt::Debug
{
    /// Set s = s + t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = Series::new("x", -1, vec!(3., 4., 5.));
    /// let res = Series::new("x", -3, vec!(1.,0.,0.));
    /// s += t;
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn add_assign(& mut self, mut other: Series<Var, C>) {
        assert_eq!(self.var, other.var);
        self.truncate_max_pow(&other);
        self.add_overlap(&other);
        if other.min_pow() < self.min_pow() {
            let num_leading = self.num_leading(&other);
            let leading_coeff = other.coeffs.drain(0..num_leading);
            self.coeffs.splice(0..0, leading_coeff);
            self.min_pow = other.min_pow;
        }
        self.trim();
    }
}

impl<'a, Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs>
    for &'a Series<Var, C>
    where Series<Var, C>: AddAssign<Rhs>
{
    type Output = Series<Var, C>;

    fn add(self, other: Rhs) -> Self::Output {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<Var: Clone, C: Coeff + Clone, Rhs> Add<Rhs>
    for Series<Var, C>
    where Series<Var, C>: AddAssign<Rhs>
{
    type Output = Series<Var, C>;

    fn add(mut self, other: Rhs) -> Self::Output {
        self += other;
        self
    }
}

impl<'a, Var, C: Coeff>
    SubAssign<&'a Series<Var, C>>
    for Series<Var, C>
where
    for<'c> &'c Series<Var, C>: Neg<Output=Series<Var, C>>,
    Series<Var, C>: AddAssign<Series<Var, C>>
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!());
    /// s -= &s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(& mut self, other: &'a Series<Var, C>) {
        *self += -other;
    }
}

impl<Var, C: Coeff>
    SubAssign<Series<Var, C>> for Series<Var, C>
where
    Series<Var, C>: AddAssign + Neg<Output=Series<Var, C>>
{
    /// Set s = s - t for two series s and t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!());
    /// s -= s.clone();
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub_assign(& mut self, other: Series<Var, C>) {
        *self += -other;
    }
}

impl<'a, 'b, Var, C: Coeff>
    Sub<&'b Series<Var, C>> for &'a Series<Var, C>
where
    for<'c> &'c Series<Var, C>:
    Add<Series<Var, C>, Output = Series<Var, C>>
    + Neg<Output = Series<Var,C>>
{
    type Output = Series<Var, C>;

    /// Subtract two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", 0, vec!());
    /// assert_eq!(res, &s - &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub(self, other: &'b Series<Var, C>) -> Self::Output {
        self + (-other)
    }
}

impl<'a, Var, C: Coeff>
    Sub<Series<Var, C>> for &'a Series<Var, C>
where
    for<'b> &'b Series<Var, C>:
    Add<Series<Var, C>, Output = Series<Var, C>>,
    Series<Var, C>: Neg<Output = Series<Var,C>>
{
    type Output = Series<Var, C>;

    /// Subtract two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", 0, vec!());
    /// assert_eq!(res, &s - &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub(self, other: Series<Var, C>) -> Self::Output {
        self + (-other)
    }
}

impl<'a, Var, C: Coeff>
    Sub<&'a Series<Var, C>> for Series<Var, C>
where
    for<'b> Series<Var, C>: SubAssign<&'b Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Subtract two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", 0, vec!());
    /// assert_eq!(res, &s - &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub(mut self, other: &'a Series<Var, C>) -> Self::Output {
        self -= other;
        self
    }
}

impl<Var, C: Coeff> Sub for Series<Var, C>
   where Series<Var, C>: SubAssign<Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Subtract two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", 0, vec!());
    /// assert_eq!(res, &s - &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn sub(mut self, other: Series<Var, C>) -> Series<Var, C> {
        self -= other;
        self
    }
}

impl<'a, Var: PartialEq + fmt::Debug, C: Coeff + Clone + AddAssign>
    MulAssign<&'a Series<Var, C>>
    for Series<Var, C>
where
    for<'b> &'b C: Mul<Output=C>,
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(& mut self, other: &'a Series<Var, C>) {
        assert_eq!(self.var, other.var);
        self.min_pow += other.min_pow();
        let num_coeffs = min(self.coeffs.len(), other.coeffs.len());
        self.coeffs.truncate(num_coeffs);
        // compute Cauchy product
        // cloning here is inefficient, but the borrow checker won't allow us
        // to modify self.coeffs directly
        let mut c: Vec<_> = self.coeffs.iter()
            .map(|c| c * &other.coeffs[0])
            .collect();
        for (k,c) in c.iter_mut().enumerate() {
            for i in 1..=k {
                *c += &self.coeffs[k-i] * &other.coeffs[i]
            }
        }
        self.coeffs = c;
    }
}

impl<Var, C: Coeff> MulAssign
    for Series<Var, C>
where for<'a> Series<Var, C>: MulAssign<&'a Series<Var, C>>
{
    /// Set s = s * t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s *= &s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul_assign(& mut self, other: Series<Var, C>) {
        *self *= &other
    }
}

impl<Var, C: Coeff> Mul for Series<Var, C>
where Series<Var, C>: MulAssign<Series<Var, C>> {
    type Output = Series<Var, C>;

    /// Multiply two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s * t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul(mut self, other: Series<Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<&'a Series<Var, C>> for Series<Var, C>
where Series<Var, C>: MulAssign<&'a Series<Var, C>> {
    type Output = Series<Var, C>;

    /// Multiply two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, s * &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul(mut self, other: &'a Series<Var, C>) -> Self::Output {
        self *= other;
        self
    }
}

impl<'a, Var, C: Coeff> Mul<Series<Var, C>> for &'a Series<Var, C>
where Series<Var, C>: MulAssign<&'a Series<Var, C>> {
    type Output = Series<Var, C>;

    /// Multiply two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, &s * t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul(self, other: Series<Var, C>) -> Self::Output {
        other * self
    }
}

impl<'a, 'b, Var: Clone + PartialEq + fmt::Debug,
     C: Coeff + Mul<Output=C> + AddAssign>
    Mul<&'b Series<Var, C>>
    for &'a Series<Var, C>
where
    for<'c, 'd> &'c C: Mul<&'d C,Output=C>
{
    type Output = Series<Var, C>;

    /// Multiply two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let t = s.clone();
    /// let res = Series::new("x", -6, vec!(1.,0.,-6.));
    /// assert_eq!(res, &s * &t);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn mul(self, other: &'b Series<Var, C>) -> Self::Output {
        assert_eq!(self.var, other.var);
        let product_min_pow = self.min_pow() + other.min_pow();
        let num_coeffs = min(self.coeffs.len(), other.coeffs.len());
        // compute Cauchy product
        let mut c: Vec<_> = self.coeffs.iter().take(num_coeffs)
            .map(|c| c * &other.coeffs[0])
            .collect();
        for (k,c) in c.iter_mut().enumerate() {
            for i in 1..=k {
                *c += &self.coeffs[k-i] * &other.coeffs[i]
            }
        }
        Series::new(self.var.clone(), product_min_pow, c)
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign>
    DivAssign<&'a Series<Var, C>>
    for Series<Var, C>
where
    Series<Var, C>: MulAssign,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>,
    for<'c> &'c Series<Var, C>: MulInverse<Output=Series<Var, C>>
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s /= &s.clone();
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(& mut self, other: &'a Series<Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<'a, Var: Clone, C: Coeff + SubAssign>
    DivAssign
    for Series<Var, C>
where
    Series<Var, C>: MulAssign + MulInverse<Output = Series<Var,C>>,
    for<'b> &'b C: Div<Output = C> + Mul<Output = C>
{
    /// Sets s = s / t for two series s,t
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// s /= s.clone();
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div_assign(& mut self, other: Series<Var, C>) {
        *self *= other.mul_inverse();
    }
}

impl<Var, C: Coeff> Div for Series<Var, C>
where Series<Var, C>: DivAssign
{
    type Output = Series<Var, C>;

    /// Divides two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s.clone() / s.clone());
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div(mut self, other: Series<Var, C>) -> Series<Var, C> {
        self /= other;
        self
    }
}

impl<'a, Var, C: Coeff> Div<&'a Series<Var, C>> for Series<Var, C>
where Series<Var, C>: DivAssign<&'a Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Divides two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s.clone() / &s.clone());
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div(mut self, other: &'a Series<Var, C>) -> Self::Output {
        self /= other;
        self
    }
}

impl<'a, Var, C: Coeff> Div<Series<Var, C>> for &'a Series<Var, C>
where Series<Var, C>: Clone + Div<Output = Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Divides two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, s.clone() / &s.clone());
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div(self, other: Series<Var, C>) -> Self::Output {
        self.clone() / other
    }
}

impl<'a, 'b, Var, C: Coeff> Div<&'b Series<Var, C>> for &'a Series<Var, C>
where
    for<'c> Series<Var, C>: Clone + Div<&'c Series<Var,C>, Output=Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Divides two series
    ///
    /// # Example
    ///
    /// ```rust
    /// use series::Series;
    /// let s = Series::new("x", -3, vec!(1.,0.,-3.));
    /// let res = Series::new("x", 0, vec!(1.,0.,0.));
    /// assert_eq!(res, &s / &s);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables.
    fn div(self, other: &'b Series<Var, C>) -> Self::Output {
        self.clone() / other
    }
}

// TODO: understand why there is a compiler error when removing this trait
//       (with rust 1.28, ok with 1.24)
//       and just implement the method
trait ExpCoeff {
    type Output;
    fn exp_coeff(&self) -> Self::Output;
}

impl<Var, C: Coeff> ExpCoeff for Series<Var, C>
where
    for<'a> &'a C: Mul<Output=C>,
    for<'a> C: MulAssign<&'a C>,
    C: Clone + From<i32> + Div<Output=C> + Mul<Output=C>
    + AddAssign + Exp<Output=C>
{
    type Output = Vec<C>;

    fn exp_coeff(&self) -> Vec<C> {
        assert!(self.min_pow() >= 0);
        let mut b = Vec::with_capacity(min(self.coeffs.len(), 1));
        b.push(C::from(1));
        debug_assert!(self.max_pow() >= 0);
        for n in 1..self.max_pow() as usize {
            let mut b_n = C::from(0);
            for i in 1..n+1 {
                let num_factor = C::from(i as i32)/C::from(n as i32);
                let a_i = self.coeff(i as isize).unwrap();
                b_n += num_factor*(a_i * &b[n-i]);
            }
            b.push(b_n);
        }
        if self.min_pow() == 0 {
            let exp_a_0 = self.coeff(0).unwrap().clone().exp();
            for mut b_n in & mut b {
                *b_n *= &exp_a_0;
            }
        }
        b
    }
}

impl<Var, C: Coeff> Exp for Series<Var, C>
where
    for<'a> &'a C: Mul<Output=C>,
    for<'a> C: MulAssign<&'a C>,
    C: Clone + From<i32> + Div<Output=C> + Mul<Output=C>
    + AddAssign + Exp<Output=C>
{
    type Output = Self;

    /// Computes the exponential of a series
    ///
    /// # Panics
    ///
    /// Panics if the series contains negative powers of the expansion
    /// variable
    fn exp(self) -> Self::Output {
        let coeff = self.exp_coeff();
        Series::new(self.var, 0, coeff)
    }
}

impl<'a, Var, C: Coeff> Exp for &'a Series<Var, C>
where
    for<'b> &'b C: Mul<Output=C>,
    for<'b> C: MulAssign<&'b C>,
    Var: Clone,
    C: Clone + From<i32> + Div<Output=C> + Mul<Output=C>
    + AddAssign + Exp<Output=C>
{
    type Output = Series<Var, C>;

    /// Computes the exponential of a series
    ///
    /// # Panics
    ///
    /// Panics if the series contains negative powers of the expansion
    /// variable
    fn exp(self) -> Self::Output {
        Series::new(self.var.clone(), 0,  self.exp_coeff())
    }
}

impl<Var, C: Coeff> Ln for Series<Var, C>
where
    for <'a> C: DivAssign<&'a C>,
    for <'a> &'a C: Mul<Output=C>,
    C: Clone + From<i32>
    + SubAssign
    + Add<Output=C>
    + Mul<Output=C>
    + Div<Output=C>
    + Ln<Output=C>
    + From<Var>,
    Var: Clone,
{
    type Output = Self;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(mut self) -> Self {
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = self.coeffs[0].clone();
        self.coeffs[0] = C::from(1);
        for i in 1..self.coeffs.len() {
            self.coeffs[i] /= &c_k0;
        }
        let a = self.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = if k0 != 0 {
           let var = self.var.clone();
            c_k0.ln() + C::from(k0 as i32)*C::from(var).ln()
        }
        else {
            c_k0.ln()
        };
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone());
            for i in 1..n {
                let num_factor = C::from(i as i32)/C::from(n as i32);
                b[n] -= num_factor * (&a[n-i] * &b[i]);
            }
        }
        Series::new(self.var, 0, b)
    }
}

impl<'a, Var, C: Coeff> Ln for &'a Series<Var, C>
where
    for <'b> C: Div<&'b C, Output=C>,
    for <'b> &'b C: Mul<Output=C> + Ln<Output=C>,
    C: Clone + From<i32>
    + SubAssign
    + Add<Output=C>
    + Mul<Output=C>
    + Div<Output=C>
    + From<Var>,
    Var: Clone,
{
    type Output = Series<Var, C>;

    /// Computes the logarithm of a series
    ///
    /// # Panics
    ///
    /// Panics if the series has no (non-zero) coefficients
    fn ln(self) -> Self::Output {
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = &self.coeffs[0];
        // self.coeffs[0] = C::from(1);
        // for i in 1..self.coeffs.len() {
        //     self.coeffs[i] /= &c_k0;
        // }
        let a = &self.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = if k0 != 0 {
           let var = self.var.clone();
            c_k0.ln() + C::from(k0 as i32)*C::from(var).ln()
        }
        else {
            c_k0.ln()
        };
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone()/c_k0);
            for i in 1..n {
                let num_factor = C::from(i as i32)/C::from(n as i32);
                b[n] -= num_factor * (&a[n-i] * &b[i])/c_k0;
            }
        }
        Series::new(self.var.clone(), 0, b)
    }
}


impl<Var, C: Coeff> Pow<Series<Var, C>> for Series<Var, C>
where Series<Var, C>: Ln<Output=Self> + Exp<Output=Self> + Mul<Output=Self>
{
    type Output = Self;

    /// Computes s^t for two series s,t
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables, s is
    /// zero, or t contains negative powers of the expansion variable
    fn pow(self, exponent: Series<Var, C>) -> Self {
        (exponent * self.ln()).exp()
    }
}

impl<'a, Var, C: Coeff> Pow<Series<Var, C>> for &'a Series<Var, C>
where
    for <'b> &'b Series<Var, C>: Ln<Output=Series<Var, C>>,
    Series<Var, C>: Exp<Output=Series<Var, C>> + Mul<Output=Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Computes s^t for two series s,t
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables, s is
    /// zero, or t contains negative powers of the expansion variable
    fn pow(self, exponent: Series<Var, C>) -> Self::Output {
        (exponent * self.ln()).exp()
    }
}

impl<'a, Var, C: Coeff> Pow<&'a Series<Var, C>> for Series<Var, C>
where
    for <'b> &'b Series<Var, C>: Mul<Series<Var, C>, Output=Series<Var, C>>,
    Series<Var, C>: Ln<Output=Self> + Exp<Output=Self>
{
    type Output = Self;

    /// Computes s^t for two series s,t
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables, s is
    /// zero, or t contains negative powers of the expansion variable
    fn pow(self, exponent: &'a Series<Var, C>) -> Self {
        (exponent * self.ln()).exp()
    }
}

impl<'a, 'b, Var, C: Coeff> Pow<&'a Series<Var, C>> for &'b Series<Var, C>
where
    for <'c> &'c Series<Var, C>:
    Ln<Output=Series<Var, C>>
    + Mul<Series<Var, C>, Output=Series<Var, C>>,
    Series<Var, C>: Exp<Output=Series<Var, C>>
{
    type Output = Series<Var, C>;

    /// Computes s^t for two series s,t
    ///
    /// # Panics
    ///
    /// Panics if the series have different expansion variables, s is
    /// zero, or t contains negative powers of the expansion variable
    fn pow(self, exponent: &'a Series<Var, C>) -> Self::Output {
        (exponent * self.ln()).exp()
    }
}

impl<'a, Var, C: Coeff>
    MulAssign<&'a C>
    for Series<Var, C>
    where C: MulAssign<&'a C>
{
    fn mul_assign(& mut self, rhs: &'a C) {
        for mut coeff in &mut self.coeffs {
            *coeff *= rhs
        }
    }
}

impl<Var, C: Coeff>
    MulAssign<C>
    for Series<Var, C>
    where for <'a> Series<Var, C>: MulAssign<&'a C>
{
    fn mul_assign(& mut self, rhs: C) {
        *self *= &rhs
    }
}

impl<'a, Var, C: Coeff>
    Mul<&'a C>
    for Series<Var, C>
    where C: MulAssign<&'a C>
{
    type Output = Self;

    fn mul(mut self, rhs: &'a C) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a, 'b, Var, C: Coeff>
    Mul<&'a C>
    for &'b Series<Var, C>
where
    C: MulAssign<&'a C>,
    Series<Var, C>: Clone
{
    type Output = Series<Var, C>;

    fn mul(self, rhs: &'a C) -> Self::Output {
        self.clone() * rhs
    }
}

impl<Var, C: Coeff>
    Mul<C>
    for Series<Var, C>
    where Series<Var, C>: MulAssign<C>
{
    type Output = Self;

    fn mul(mut self, rhs: C) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<'a, Var, C: Coeff>
    Mul<C>
    for &'a Series<Var, C>
where
    Series<Var, C>: Clone + MulAssign<C>
{
    type Output = Series<Var, C>;

    fn mul(self, rhs: C) -> Self::Output {
        self.clone() * rhs
    }
}

impl<'a, Var, C: Coeff>
    DivAssign<&'a C>
    for Series<Var, C>
    where C: DivAssign<&'a C>
{
    fn div_assign(& mut self, rhs: &'a C) {
        for mut coeff in &mut self.coeffs {
            *coeff /= rhs
        }
    }
}

impl<Var, C: Coeff>
    DivAssign<C>
    for Series<Var, C>
    where for <'a> Series<Var, C>: DivAssign<&'a C>
{
    fn div_assign(& mut self, rhs: C) {
        *self /= &rhs
    }
}

impl<'a, Var, C: Coeff>
    Div<&'a C>
    for Series<Var, C>
    where C: DivAssign<&'a C>
{
    type Output = Self;

    fn div(mut self, rhs: &'a C) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<'a, 'b, Var, C: Coeff>
    Div<&'a C>
    for &'b Series<Var, C>
where
    C: DivAssign<&'a C>,
    Series<Var, C>: Clone
{
    type Output = Series<Var, C>;

    fn div(self, rhs: &'a C) -> Self::Output {
        self.clone() / rhs
    }
}

impl<Var, C: Coeff>
    Div<C>
    for Series<Var, C>
    where Series<Var, C>: DivAssign<C>
{
    type Output = Self;

    fn div(mut self, rhs: C) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<'a, Var, C: Coeff>
    Div<C>
    for &'a Series<Var, C>
where
    Series<Var, C>: Clone + DivAssign<C>
{
    type Output = Series<Var, C>;

    fn div(self, rhs: C) -> Self::Output {
        self.clone() / rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tst_series() {
        let var = String::from("x");
        let min_pow = -10;
        let coeffs = vec!();
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-11), Some(&0));
        assert_eq!(s.coeff(-10), None);

        let min_pow = -3;
        let coeffs = vec!(1.,2.,3.);
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-4), Some(&0.));
        assert_eq!(s.coeff(-3), Some(&1.));
        assert_eq!(s.coeff(-2), Some(&2.));
        assert_eq!(s.coeff(-1), Some(&3.));
        assert_eq!(s.coeff(0), None);

        let min_pow = -2;
        let coeffs = vec!(0.,0.,3.);
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow + 2);
        assert_eq!(s.coeff(-2), Some(&0.));
        assert_eq!(s.coeff(-1), Some(&0.));
        assert_eq!(s.coeff(0), Some(&3.));
        assert_eq!(s.coeff(1), None);

        let s = Series::new(var.clone(), -2, vec!(0.,0.,1.));
        let t = Series::new(var.clone(), 0, vec!(1.));
        assert_eq!(s,t);

        let s = Series::new(var.clone(), -3, vec!(0.,0.,0.));
        let t = Series::new(var.clone(), 0, vec!());
        assert_eq!(s,t);
    }

    #[test]
    fn tst_display() {
        // let s = Series::new("x", -10, vec!());
        // assert_eq!(format!("{}", s), "O(x^-10)");
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        assert_eq!(format!("{}", s), "(1)*x^-3 + (-3)*x^-1 + O(x^0)");
        let s = Series::new("x", -1, vec!(1.,2.,-3.));
        assert_eq!(format!("{}", s), "(1)*x^-1 + (2) + (-3)*x + O(x^2)");
    }

    #[test]
    fn tst_neg() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", -3, vec!(-1.,0.,3.));
        assert_eq!(res, -&s);
        assert_eq!(res, -s);
    }

    #[test]
    fn tst_add() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", -3, vec!(2.,0.,-6.));
        assert_eq!(res, &s + &s);
        assert_eq!(res, &s + s.clone());
        assert_eq!(res, s.clone() + &s);
        assert_eq!(res, s.clone() + s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        assert_eq!(s, &s + &t);
        assert_eq!(s, &t + &s);
        assert_eq!(s, &s + t.clone());
        assert_eq!(s, &t + s.clone());
        assert_eq!(s, s.clone() + &t);
        assert_eq!(s, t.clone() + &s);
        assert_eq!(s, s.clone() + t.clone());
        assert_eq!(s, t.clone() + s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -3, vec!(-1., 0., 3.));
        let res = Series::new("x", 0, vec!());
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());
    }

    #[test]
    fn tst_add_assign() {
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", -3, vec!(2.,0.,-6.));
        s += s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s += t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -1, vec!(3.,4.,5.));
        let t = Series::new("x", -3, vec!(1.,0.,-3.));
        s += t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = s.clone();
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        s += &t;
        assert_eq!(s, res);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s += t;
        assert_eq!(s, res);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -3, vec!(-1., 0., 3.));
        let res = Series::new("x", 0, vec!());
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s += t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_sub() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", 0, vec!());
        assert_eq!(res, &s - &s);
        assert_eq!(res, &s - s.clone());
        assert_eq!(res, s.clone() - &s);
        assert_eq!(res, s.clone() - s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(-3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        assert_eq!(res, &s - &t);
        assert_eq!(res, &s - t.clone());
        assert_eq!(res, s.clone() - &t);
        assert_eq!(res, s - t);

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        assert_eq!(s, &s - &t);
        assert_eq!(s, &s - t.clone());
        assert_eq!(s, s.clone() - &t);
        assert_eq!(s, s.clone() - t);
    }

    #[test]
    fn tst_sub_assign() {
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", 0, vec!());
        s -= s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(-3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s -= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = s.clone();
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s -= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_mul() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", -6, vec!(1.,0.,-6.));
        assert_eq!(res, &s * &s);
        assert_eq!(res, &s * s.clone());
        assert_eq!(res, s.clone() * &s);
        assert_eq!(res, s.clone() * s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5., 7.));
        let res = Series::new("x", -4, vec!(3.,4.,-4.));
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);

        let s = Series::new("x", -3, vec!(1., 7.,-3.));
        let t = Series::new("x", 3, vec!(1., -7., 52.));
        let res = Series::new("x", 0, vec!(1.,0., 0.));
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, s.clone() * &t);
        assert_eq!(res, t.clone() * &s);
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);
    }

    #[test]
    fn tst_mul_assign() {
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s *= s.clone();
        let res = Series::new("x", -6, vec!(1.,0.,-6.));
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5., 7.));
        s *= &t;
        let res = Series::new("x", -4, vec!(3.,4.,-4.));
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s *= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1., 7.,-3.));
        let t = Series::new("x", 3, vec!(1., -7., 52.));
        s *= &t;
        let res = Series::new("x", 0, vec!(1.,0., 0.));
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1., 7.,-3.));
        s *= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_div() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let res = Series::new("x", 0, vec!(1.,0.,0.));
        assert_eq!(res, &s / &s);
        assert_eq!(res, &s / s.clone());
        assert_eq!(res, s.clone() / &s);
        assert_eq!(res, s.clone() / s.clone());

        // disabled for floating-point inaccuracy
        // let s = Series::new("x", -3, vec!(1.,0.,-3.));
        // let t = Series::new("x", -1, vec!(3., 4., 5., 7.));
        // let res = Series::new("x", -2, vec!(1./3.,-4./9.,-26./27.));
        // assert_eq!(res, &s / &t);
        // assert_eq!(res, &t / &s);
        // assert_eq!(res, s / t);

        let s = Series::new("x", -3, vec!(1., 7.,-3.));
        let t = Series::new("x", 3, vec!(1., -7., 52.));
        let res = Series::new("x", -6, vec!(1.,14., 43.));
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);

        let s = Series::new("x", 1, vec!(1., 7.,-3.));
        let t = Series::new("x", 5, vec!());
        let res = Series::new("x", -4, vec!());
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);
    }

    #[test]
    fn tst_div_assign() {
        let mut s = Series::new("x", -3, vec!(1.,0.,-3.));
        s /= s.clone();
        let res = Series::new("x", 0, vec!(1.,0.,0.));
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec!(1., 7.,-3.));
        let t = Series::new("x", 3, vec!(1., -7., 52.));
        s /= &t;
        let res = Series::new("x", -6, vec!(1.,14., 43.));
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec!(1., 7.,-3.));
        s /= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", 1, vec!(1., 7.,-3.));
        let t = Series::new("x", 5, vec!());
        s /= &t;
        let res = Series::new("x", -4, vec!());
        assert_eq!(res, s);
        let mut s = Series::new("x", 1, vec!(1., 7.,-3.));
        s /= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_var() {
        let _ = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let _ = Series::new('j', -3, vec!(1.,0.,-3.));
        let _ = Series::new(8, -3, vec!(1.,0.,-3.));
    }

    #[derive(Debug,Clone,PartialEq)]
    struct Mystr<'a>(&'a str);

    impl<'a> From<Mystr<'a>> for f64 {
        fn from(_s: Mystr<'a>) -> f64 {
            panic!("can't turn str to f64")
        }
    }

    #[test]
    fn tst_ln() {
        let s = Series::new(Mystr("x"), 0, vec!(1., 7.,-3.));
        let res = Series::new(Mystr("x"), 1, vec!(7.,-55./2.));
        assert_eq!(res, (&s).ln());
        //assert_eq!(res, s.ln());

        let s = Series::new(Mystr("x"), 0, vec!(4., 7.,-3.));
        let res = Series::new(Mystr("x"), 0, vec!(4_f64.ln(),7./4.,-73./32.));
        assert_eq!(res, (&s).ln());
        //assert_eq!(res, s.ln());
    }

    #[test]
    fn tst_exp() {
        let s = Series::new("x", 1, vec!(7.,-3.));
        let res = Series::new("x", 0, vec!(1., 7.,43./2.));
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 2, vec!(0.));
        let res = Series::new("x", 0, vec!(1.,0.,0.));
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 0, vec!(5., 11., -7.));
        let e5 = 5_f64.exp();
        let res = Series::new("x", 0, vec!(e5,e5*11.,e5*107./2.));
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());
    }

    #[test]
    fn tst_pow() {
        let base = Series::new(Mystr("x"), 0, vec!(1., 7., 0.));
        let exp = Series::new(Mystr("x"), -1, vec!(1., -5., 43.));
        let e7 = 7_f64.exp();
        let res = Series::new(Mystr("x"), 0, vec!(e7, -119./2.*e7));
        assert_eq!(res, (&base).pow(&exp));
        assert_eq!(res, (&base).pow(exp.clone()));
        assert_eq!(res, base.clone().pow(&exp));
        assert_eq!(res, base.pow(exp));

        let base = Series::new(Mystr("x"), 0, vec!(2., 7., 0.));
        let exp = Series::new(Mystr("x"), 0, vec!(3., -5., 11.));
        // rescale result so we can use round and still get decent precision
        let rescale = Series::new(Mystr("x"), 0, vec!(1e13, 0., 0., 0.));
        let test = &rescale * &base.pow(exp);
        let ln2 = 2_f64.ln();
        let res = Series::new(Mystr("x"), 0, vec!(8., 84.-40.*ln2, 154.+ln2*(-332.+100.*ln2)));
        let res = rescale * res;
        assert_eq!(res.min_pow(), test.min_pow());
        assert_eq!(res.max_pow(), test.max_pow());
        for i in res.min_pow()..res.max_pow() {
            assert_eq!(
                res.coeff(i).unwrap().round(),
                test.coeff(i).unwrap().round()
            );
        }
    }

    #[test]
    fn tst_scalar() {
        let s = Series::new("x", -3, vec!(1.,0.,-2.));
        let res = Series::new("x", -3, vec!(1./2.,0.,-1.));
        assert_eq!(res, &s / 2.);
        let mut s = s;
        s /= 2.;
        assert_eq!(res, s);

        let s = Series::new("x", -3, vec!(1./2.,0.,-1.));
        let res = Series::new("x", -3, vec!(1.,0.,-2.));
        assert_eq!(res, &s * 2.);
        let mut s = s;
        s *= 2.;
        assert_eq!(res, s);
    }

    #[test]
    #[should_panic]
    fn tst_bad_add() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        let _ = s + t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_sub() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        let _ = s - t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_mul() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        let _ = s * t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_div() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        let _ = s / t;
    }
}
