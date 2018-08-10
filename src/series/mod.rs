use std::fmt;
use std::ops::{Add,AddAssign,Sub,SubAssign,Mul,MulAssign,Div,DivAssign,Neg};
use std::convert::From;
use std::cmp::min;

pub mod ops;
use self::ops::{Ln,Exp};

type Pow = i32;

#[derive(PartialEq,Debug,Clone)]
struct Series<Var, Coeff> {
    var: Var,
    min_pow: Pow,
    coeffs: Vec<Coeff>,
    zero: Coeff // TODO: evil hack, learn how to do this better
}

impl<Var, Coeff> Series<Var, Coeff> {
    pub fn min_pow(&self) -> Pow {
        self.min_pow
    }

    pub fn max_pow(&self) -> Pow {
        // TODO: replace i32 by bigger type
        self.min_pow + Pow::from(self.coeffs.len() as i32)
    }
}

impl<Var, Coeff: From<i32> + PartialEq> Series<Var, Coeff> {
    pub fn new(var: Var, min_pow: Pow, coeffs: Vec<Coeff>) -> Series<Var, Coeff> {
        let mut res = Series{var, min_pow, coeffs, zero: Coeff::from(0)};
        res.trim();
        res
    }

    pub fn coeff(&self, pow: Pow) -> Option<&Coeff> {
        if pow < self.min_pow() {
            return Some(&self.zero) // TODO this is a bad hack
        }
        if pow >= self.max_pow() {
            return None
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(&self.coeffs[idx])
    }
}

impl<
        Var: Clone, Coeff: From<i32> + PartialEq + SubAssign
    > Series<Var, Coeff> where
    for<'a> &'a Coeff: Div<Output = Coeff> + Mul<Output = Coeff>
{
    pub fn inverse(&self) -> Series<Var, Coeff> {
        let inv_min_pow = -self.min_pow;
        if self.coeffs.is_empty() {
            return Series::new(self.var.clone(), inv_min_pow, vec!())
        }
        let a : Vec<_> = self.coeffs.iter().map(|c| c/&self.coeffs[0]).collect();
        let mut b = Vec::with_capacity(a.len());
        b.push(Coeff::from(1));
        for n in 1..a.len() {
            let mut b_n = Coeff::from(0);
            for i in 0..n {
                b_n -= &a[n-i] * &b[i];
            }
            b.push(b_n);
        }
        let inv_coeffs : Vec<_> = b.iter().map(|b| b/&self.coeffs[0]).collect();
        Series::new(self.var.clone(), inv_min_pow, inv_coeffs)
    }
}

impl<Var, Coeff: From<i32> + PartialEq> Series<Var, Coeff> {
    fn trim(& mut self) {
        let idx = {
            let non_zero = self.coeffs.iter().enumerate().
                find(|&(_, c)| *c != Coeff::from(0));
            match non_zero {
                Some((idx, _)) => Some(idx),
                _ => None
            }
        };
        match idx {
            Some(idx) => {
                if idx > 0 {
                    self.min_pow += Pow::from(idx as i32);
                    self.coeffs.drain(0..idx);
                }
            },
            None => {
                self.min_pow += Pow::from(self.coeffs.len() as i32);
                self.coeffs = vec!();
            }
        }
    }
}

impl<Var: fmt::Display, Coeff: From<i32> + PartialEq + fmt::Display> fmt::Display
    for Series<Var, Coeff>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if *coeff == Coeff::from(0) { continue; }
            let cur_pow = self.min_pow() + Pow::from(i as i32);
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

impl<
    Var, Coeff: From<i32> + PartialEq
>
    Neg for Series<Var, Coeff>
    where for<'c> &'c Coeff: Neg<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn neg(self) -> Series<Var, Coeff> {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Series::new(self.var, self.min_pow, neg_coeff)
    }
}

impl<'a, Var: Clone, Coeff: From<i32> + PartialEq>
    Neg for &'a Series<Var, Coeff>
    where for<'c> &'c Coeff: Neg<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn neg(self) -> Series<Var, Coeff> {
        let neg_coeff = self.coeffs.iter().map(|c| -c).collect();
        Series::new(self.var.clone(), self.min_pow, neg_coeff)
    }
}

impl<
    'a,
    Var: PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq + Clone
>
    AddAssign<&'a Series<Var, Coeff>>
    for Series<Var, Coeff>
    where for<'c> Coeff: AddAssign<&'c Coeff>
{
    fn add_assign(& mut self, other: &'a Series<Var, Coeff>) {
        assert_eq!(self.var, other.var);
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
        let offset = self.min_pow();
        for (i, c) in self.coeffs.iter_mut().enumerate() {
            let power = offset + i as i32;
            *c += other.coeff(power).unwrap();
        }
        if other.min_pow() < self.min_pow() {
            let num_leading = min(
                (self.min_pow() - other.min_pow()) as usize,
                other.coeffs.len()
            );
            let leading_coeff = other.coeffs[0..num_leading].iter().cloned();
            self.coeffs.splice(0..0, leading_coeff);
            self.min_pow = other.min_pow;
        }
        debug_assert!(other.max_pow() >= self.max_pow());
        self.trim();
    }
}

impl<Var, Coeff>
    AddAssign<Series<Var, Coeff>> for Series<Var, Coeff>
    where for<'c> Series<Var, Coeff>: AddAssign<&'c Series<Var, Coeff>>
{
    fn add_assign(& mut self, other: Series<Var, Coeff>) {
        self.add_assign(&other)
    }
}

impl<'a, 'b, Var: Clone, Coeff: Clone> Add<&'b Series<Var, Coeff>>
    for &'a Series<Var, Coeff>
    where for<'c> Series<Var, Coeff>: AddAssign<&'c Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn add(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        let mut res = self.clone();
        res += other;
        res
    }
}

impl<Var, Coeff> Add for Series<Var, Coeff>
    where Series<Var, Coeff>: AddAssign<Series<Var,Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn add(mut self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        self += other;
        self
    }
}

impl<'a, 'b, Var, Coeff>
    Sub<&'b Series<Var, Coeff>> for &'a Series<Var, Coeff>
    where for<'c> &'c Series<Var, Coeff>:
    Add<Output = Series<Var, Coeff>> + Neg<Output = Series<Var,Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn sub(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        self + &(-other)
    }
}

impl<Var, Coeff> Sub for Series<Var, Coeff>
   where Series<Var, Coeff>: SubAssign<Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn sub(mut self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        self -= other;
        self
    }
}

impl<'a, Var, Coeff>
    SubAssign<&'a Series<Var, Coeff>>
    for Series<Var, Coeff>
where
    for<'c> &'c Series<Var, Coeff>: Neg<Output=Series<Var, Coeff>>,
    Series<Var, Coeff>: AddAssign<Series<Var, Coeff>>
{
    fn sub_assign(& mut self, other: &'a Series<Var, Coeff>) {
        *self += -other;
    }
}

impl<Var, Coeff>
    SubAssign<Series<Var, Coeff>> for Series<Var, Coeff>
    where for<'c> Series<Var, Coeff>: SubAssign<&'c Series<Var, Coeff>>
{
    fn sub_assign(& mut self, other: Series<Var, Coeff>) {
        *self -= &other;
    }
}


impl<
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq + AddAssign
>
    Mul<&'b Series<Var, Coeff>>
    for &'a Series<Var, Coeff>
    where for<'c> &'c Coeff: Mul<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn mul(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        assert_eq!(self.var, other.var); // TODO: handle this better
        let res_min_pow = self.min_pow() + other.min_pow();
        let num_coeffs = min(self.coeffs.len(), other.coeffs.len());
        // compute Cauchy product
        let mut res_coeff = Vec::with_capacity(num_coeffs);
        for k in 0..num_coeffs {
            let mut c_k = &self.coeffs[0] * &other.coeffs[k];
            for i in 1..(k+1) {
                c_k += &self.coeffs[i] * &other.coeffs[k-i];
            }
            res_coeff.push(c_k);
        }
        Series::new(self.var.clone(), res_min_pow, res_coeff)
    }
}

impl<Var, Coeff> Mul for Series<Var, Coeff>
where for<'a> &'a Series<Var, Coeff>: Mul<Output = Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn mul(self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        &self * &other
    }
}

impl<'a, Var, Coeff>
    MulAssign<&'a Series<Var, Coeff>>
    for Series<Var, Coeff>
where for<'c> &'c Series<Var, Coeff>: Mul<Output=Series<Var, Coeff>>
{
    fn mul_assign(& mut self, other: &'a Series<Var, Coeff>) {
        let res = &*self * other;
        *self = res
    }
}

impl<Var, Coeff: fmt::Display>
    MulAssign<Series<Var, Coeff>>
    for Series<Var, Coeff>
    where for <'c> Series<Var, Coeff>: MulAssign<&'c Series<Var, Coeff>>
{
    fn mul_assign(& mut self, other: Series<Var, Coeff>) {
        *self *= &other;
    }
}

impl<
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq + SubAssign
>
    Div<&'b Series<Var, Coeff>>
    for &'a Series<Var, Coeff>
    where for<'c> &'c Series<Var, Coeff>: Mul<Output = Series<Var, Coeff>>,
          for<'d> &'d Coeff: Div<Output = Coeff> + Mul<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn div(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        assert_eq!(self.var, other.var); // TODO: handle this better
        let inv = other.inverse();
        self * &inv
    }
}

impl<Var, Coeff> Div for Series<Var, Coeff>
where for<'a> &'a Series<Var, Coeff>: Div<Output = Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn div(self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        &self / &other
    }
}

impl<'a, Var, Coeff>
    DivAssign<&'a Series<Var, Coeff>>
    for Series<Var, Coeff>
where for<'c> &'c Series<Var, Coeff>: Div<Output=Series<Var, Coeff>>
{
    fn div_assign(& mut self, other: &'a Series<Var, Coeff>) {
        let res = &*self / other;
        *self = res
    }
}

impl<Var, Coeff: fmt::Display>
    DivAssign<Series<Var, Coeff>>
    for Series<Var, Coeff>
    where for <'c> Series<Var, Coeff>: DivAssign<&'c Series<Var, Coeff>>
{
    fn div_assign(& mut self, other: Series<Var, Coeff>) {
        *self /= &other;
    }
}

impl<
    Var: Clone,
    Coeff: Clone + PartialEq + From<i32>
        + Add<Output=Coeff> + AddAssign<Coeff>
        + Mul<Output=Coeff> + Div<Output=Coeff> + Exp
     > Exp for Series<Var, Coeff>
where
    for <'c> &'c Coeff: Mul<Output=Coeff>,
    for <'c> Coeff: MulAssign<&'c Coeff>,
{
    fn exp(self) -> Series<Var, Coeff>{
        assert!(self.min_pow() >= 0);
        let mut b = Vec::with_capacity(min(self.coeffs.len(), 1));
        b.push(Coeff::from(1));
        debug_assert!(self.max_pow() >= 0);
        for n in 1..self.max_pow() as usize {
            let mut b_n = Coeff::from(0);
            for i in 1..n+1 {
                let num_factor = Coeff::from(i as i32)/Coeff::from(n as i32);
                let a_i = self.coeff(i as i32).unwrap();
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
        Series::new(self.var, 0, b)
    }
}

impl<
    Var: Clone,
    Coeff: Clone + PartialEq + From<i32> + From<Var>
        + Add<Output=Coeff> + SubAssign<Coeff>
        + Mul<Output=Coeff> + Div<Output=Coeff> + Ln
    > Ln for Series<Var, Coeff>
where
    for <'c> Coeff: DivAssign<&'c Coeff>,
    for <'c> &'c Coeff: Mul<Output=Coeff>,
{
    fn ln(mut self) -> Series<Var, Coeff> {
        assert!(!self.coeffs.is_empty());
        let k0 = self.min_pow();
        let c_k0 = self.coeffs[0].clone();
        self.coeffs[0] = Coeff::from(1);
        for i in 1..self.coeffs.len() {
            self.coeffs[i] /= &c_k0;
        }
        let a = self.coeffs;
        let mut b = Vec::with_capacity(a.len());
        let b_0 = if k0 != 0 {
           let var = self.var.clone();
            c_k0.ln() + Coeff::from(k0)*Coeff::from(var).ln()
        }
        else {
            c_k0.ln()
        };
        b.push(b_0);
        for n in 1..a.len() {
            b.push(a[n].clone());
            for i in 1..n {
                let num_factor = Coeff::from(i as i32)/Coeff::from(n as i32);
                b[n] -= num_factor * (&a[n-i] * &b[i]);
            }
        }
        Series::new(self.var, 0, b)
    }
}

impl<Var, Coeff> Series<Var, Coeff>
where Series<Var, Coeff>: Ln + Exp + Mul<Output=Series<Var, Coeff>>
{
    pub fn pow(self, exponent: Series<Var, Coeff>) -> Series<Var, Coeff> {
        (exponent * self.ln()).exp()
    }
}

impl<Var, Coeff> Series<Var, Coeff>
where
    for<'c> &'c Series<Var, Coeff>: Ln + Mul<Output=Series<Var, Coeff>>,
    Series<Var, Coeff>: Exp
{
    pub fn pow(self, exponent: & Series<Var, Coeff>) -> Series<Var, Coeff> {
        (exponent * self.ln()).exp()
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
        assert_eq!(res, s.clone() + s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, s + t);

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        assert_eq!(s, &s + &t);
        assert_eq!(s, &t + &s);
        assert_eq!(s, s.clone() + t);

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -3, vec!(-1., 0., 3.));
        let res = Series::new("x", 0, vec!());
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, s + t);
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
        assert_eq!(res, s.clone() - s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(-3., 4., 5.));
        let res = Series::new("x", -3, vec!(1.,0.,0.));
        assert_eq!(res, &s - &t);
        assert_eq!(res, s - t);

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", 1, vec!(3., 4., 5.));
        assert_eq!(s, &s - &t);
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
        assert_eq!(res, s.clone() * s.clone());

        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("x", -1, vec!(3., 4., 5., 7.));
        let res = Series::new("x", -4, vec!(3.,4.,-4.));
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, s * t);

        let s = Series::new("x", -3, vec!(1., 7.,-3.));
        let t = Series::new("x", 3, vec!(1., -7., 52.));
        let res = Series::new("x", 0, vec!(1.,0., 0.));
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
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
        assert_eq!(res.inverse(), &t / &s);
        assert_eq!(res, s / t);

        let s = Series::new("x", 1, vec!(1., 7.,-3.));
        let t = Series::new("x", 5, vec!());
        let res = Series::new("x", -4, vec!());
        assert_eq!(res, &s / &t);
        assert_eq!(res.inverse(), &t / &s);
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
        assert_eq!(res, s.ln());

        let s = Series::new(Mystr("x"), 0, vec!(4., 7.,-3.));
        let res = Series::new(Mystr("x"), 0, vec!(4_f64.ln(),7./4.,-73./32.));
        assert_eq!(res, s.ln());
    }

    #[test]
    fn tst_exp() {
        let s = Series::new("x", 1, vec!(7.,-3.));
        let res = Series::new("x", 0, vec!(1., 7.,43./2.));
        assert_eq!(res, s.exp());

        let s = Series::new("x", 2, vec!(0.));
        let res = Series::new("x", 0, vec!(1.,0.,0.));
        assert_eq!(res, s.exp());

        let s = Series::new("x", 0, vec!(5., 11., -7.));
        let e5 = 5_f64.exp();
        let res = Series::new("x", 0, vec!(e5,e5*11.,e5*107./2.));
        assert_eq!(res, s.exp());
    }

    #[test]
    fn tst_pow() {
        let base = Series::new(Mystr("x"), 0, vec!(1., 7., 0.));
        let exp = Series::new(Mystr("x"), -1, vec!(1., -5., 43.));
        let e7 = 7_f64.exp();
        let res = Series::new(Mystr("x"), 0, vec!(e7, -119./2.*e7));
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
    #[should_panic]
    fn tst_bad_add() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        s + t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_sub() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        s - t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_mul() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        s * t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_div() {
        let s = Series::new("x", -3, vec!(1.,0.,-3.));
        let t = Series::new("y", -3, vec!(1.,0.,-3.));
        s / t;
    }
}
