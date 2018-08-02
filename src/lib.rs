use std::fmt;
use std::ops::{Add,Sub,Mul,Div};
use std::convert::From;
use std::cmp::min;

type Pow = i32;

#[derive(PartialEq,Debug,Clone)]
struct Series<Var, Coeff> {
    var: Var,
    min_pow: Pow,
    coeffs: Vec<Coeff>,
    zero: Coeff // TODO: evil hack, learn how to do this better
}

impl<Var, Coeff: From<i32> + PartialEq> Series<Var, Coeff> {
    pub fn new(var: Var, min_pow: Pow, coeffs: Vec<Coeff>) -> Series<Var, Coeff> {
        let mut res = Series{var, min_pow, coeffs, zero: Coeff::from(0)};
        res.trim();
        res
    }

    pub fn min_pow(&self) -> Pow {
        self.min_pow
    }

    pub fn max_pow(&self) -> Pow {
        // TODO: replace i32 by bigger type
        self.min_pow + Pow::from(self.coeffs.len() as i32)
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
        Var: Clone, Coeff: From<i32> + PartialEq + Sub<Output = Coeff>
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
                b_n = b_n - &a[n-i] * &b[i];
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
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq
>
    Add<&'b Series<Var, Coeff>>
    for &'a Series<Var, Coeff>
    where for<'c> &'c Coeff: Add<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn add(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        assert_eq!(self.var, other.var); // TODO: handle this better
        let res_min_pow = min(self.min_pow(), other.min_pow());
        let res_max_pow = min(self.max_pow(), other.max_pow());
        assert!(res_max_pow >= res_min_pow);
        let num_coeffs = (res_max_pow - res_min_pow) as usize;
        let mut res_coeff = Vec::with_capacity(num_coeffs);
        for idx in res_min_pow..res_max_pow {
            res_coeff.push(
                self.coeff(idx).unwrap() +  other.coeff(idx).unwrap()
            );
        }
        Series::new(self.var.clone(), res_min_pow, res_coeff)
    }
}

impl<Var, Coeff> Add for Series<Var, Coeff>
where for<'a> &'a Series<Var, Coeff>: Add<Output = Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn add(self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        &self + &other
    }
}

impl<
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq
>
    Sub<&'b Series<Var, Coeff>>
    for &'a Series<Var, Coeff>
    where for<'c> &'c Coeff: Sub<Output = Coeff>
{
    type Output = Series<Var, Coeff>;

    fn sub(self, other: &'b Series<Var, Coeff>) -> Series<Var, Coeff> {
        assert_eq!(self.var, other.var); // TODO: handle this better
        let res_min_pow = min(self.min_pow(), other.min_pow());
        let res_max_pow = min(self.max_pow(), other.max_pow());
        assert!(res_max_pow >= res_min_pow);
        let num_coeffs = (res_max_pow - res_min_pow) as usize;
        let mut res_coeff = Vec::with_capacity(num_coeffs);
        for idx in res_min_pow..res_max_pow {
            res_coeff.push(
                self.coeff(idx).unwrap() -  other.coeff(idx).unwrap()
            );
        }
        Series::new(self.var.clone(), res_min_pow, res_coeff)
    }
}

impl<Var, Coeff> Sub for Series<Var, Coeff>
where for<'a> &'a Series<Var, Coeff>: Sub<Output = Series<Var, Coeff>>
{
    type Output = Series<Var, Coeff>;

    fn sub(self, other: Series<Var, Coeff>) -> Series<Var, Coeff> {
        &self - &other
    }
}

impl<
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq + Add<Output = Coeff>
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
                c_k = c_k + &self.coeffs[i] * &other.coeffs[k-i];
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

impl<
    'a, 'b,
    Var: Clone + PartialEq + fmt::Debug,
    Coeff: From<i32> + PartialEq + Sub<Output = Coeff>
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
    fn tst_var() {
        let _ = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let _ = Series::new('j', -3, vec!(1.,0.,-3.));
        let _ = Series::new(8, -3, vec!(1.,0.,-3.));
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
