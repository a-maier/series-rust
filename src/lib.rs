use std::fmt;

type Var = String;
type Pow = i32;
type Coeff = f64;

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
        self.min_pow + Pow::from(self.coeffs.len() as i32) - 1
    }

    pub fn coeff(&self, pow: Pow) -> Option<Coeff> {
        if pow < self.min_pow() {
            return Some(Coeff::from(0))
        }
        if pow > self.max_pow() {
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
        if let Some(idx) = idx {
            if idx > 0 {
                self.min_pow += Pow::from(idx as i32);
                self.coeffs.drain(0..idx);
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
        write!(f, "O({}^{})", self.var, self.max_pow() + 1)
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
}
