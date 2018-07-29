use std::fmt;
use std::ops::Add;
use std::cmp::min;

type Var = String;
type Pow = i32;
type Coeff = f64;

#[derive(PartialEq,Debug,Clone)]
struct Series {
    var: Var,
    min_pow: Pow,
    coeffs: Vec<Coeff>
}

impl Series {
    pub fn new(var: Var, min_pow: Pow, coeffs: Vec<Coeff>) -> Series {
        let mut res = Series{var, min_pow, coeffs};
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

    pub fn coeff(&self, pow: Pow) -> Option<Coeff> {
        if pow < self.min_pow() {
            return Some(Coeff::from(0))
        }
        if pow >= self.max_pow() {
            return None
        }
        let idx = (pow - self.min_pow()) as usize;
        Some(self.coeffs[idx])
    }

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

impl fmt::Display for Series {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if *coeff == Coeff::from(0) { continue; }
            let cur_pow = self.min_pow() + Pow::from(i as i32);
            if i > 0 { write!(f, " + ")?; }
            write!(f, "({})*{}^{}", coeff, self.var, cur_pow)?;
        }
        if ! self.coeffs.is_empty() { write!(f, " + ")?; }
        write!(f, "O({}^{})", self.var, self.max_pow())
    }
}

impl<'a, 'b> Add<&'b Series> for &'a Series {
    type Output = Series;

    fn add(self, other: &'b Series) -> Series {
        // TODO: handle this better
        assert_eq!(self.var, other.var);
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

impl Add for Series {
    type Output = Series;

    fn add(self, other: Series) -> Series {
        &self + &other
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
        assert_eq!(s.coeff(-11), Some(Coeff::from(0)));
        assert_eq!(s.coeff(-10), None);

        let min_pow = -3;
        let coeffs = vec!(1.,2.,3.);
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-4), Some(Coeff::from(0)));
        assert_eq!(s.coeff(-3), Some(Coeff::from(1)));
        assert_eq!(s.coeff(-2), Some(Coeff::from(2)));
        assert_eq!(s.coeff(-1), Some(Coeff::from(3)));
        assert_eq!(s.coeff(0), None);

        let min_pow = -2;
        let coeffs = vec!(0.,0.,3.);
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow + 2);
        assert_eq!(s.coeff(-2), Some(Coeff::from(0)));
        assert_eq!(s.coeff(-1), Some(Coeff::from(0)));
        assert_eq!(s.coeff(0), Some(Coeff::from(3)));
        assert_eq!(s.coeff(1), None);
    }

    #[test]
    fn tst_display() {
        let s = Series::new(String::from("x"), -10, vec!());
        assert_eq!(format!("{}", s), "O(x^-10)");
        let s = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        assert_eq!(format!("{}", s), "(1)*x^-3 + (-3)*x^-1 + O(x^0)");
    }

    #[test]
    fn tst_add() {
        let s = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let res = Series::new(String::from("x"), -3, vec!(2.,0.,-6.));
        assert_eq!(res, &s + &s);
        assert_eq!(res, s.clone() + s.clone());

        let s = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let t = Series::new(String::from("x"), -1, vec!(3., 4., 5.));
        let res = Series::new(String::from("x"), -3, vec!(1.,0.,0.));
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, s + t);

        let s = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let t = Series::new(String::from("x"), 1, vec!(3., 4., 5.));
        assert_eq!(s, &s + &t);
        assert_eq!(s, &t + &s);
        assert_eq!(s, s.clone() + t);

        let s = Series::new(String::from("x"), -3, vec!(1.,0.,-3.));
        let t = Series::new(String::from("x"), -3, vec!(-1., 0., 3.));
        let res = Series::new(String::from("x"), 0, vec!());
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, s + t);
    }

}
