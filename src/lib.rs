#![allow(clippy::suspicious_op_assign_impl)]
#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
extern crate num_traits;
#[macro_use]
extern crate log;

pub mod ops;
pub mod poly;
pub mod polyslice;
pub mod series;
pub mod slice;
pub use self::ops::{Exp, Ln, Pow};
pub use self::poly::{Polynomial, PolynomialParts};
pub use self::polyslice::PolynomialSlice;
pub use self::series::{Series, SeriesParts};
pub use self::slice::SeriesSlice;
mod traits;
pub use self::traits::{AsSlice, KaratsubaMul, MulInverse};
mod util;

use num_traits::{One, Zero};
use std::iter::Zip;
use std::ops::RangeFrom;

/// Minimum requirements on series coefficients
pub trait Coeff: From<i32> + Zero + One + PartialEq {}
impl<T: From<i32> + Zero + One + PartialEq> Coeff for T {}

/// Immutable `Series` iterator.
///
/// This `struct` is created by the `iter` method on `Series`
pub type Iter<'a, C> = Zip<RangeFrom<isize>, std::slice::Iter<'a, C>>;
/// An iterator that moves out of a series.
///
/// This `struct` is created by the `into_iter` method on `Series`
pub type IntoIter<C> = Zip<RangeFrom<isize>, std::vec::IntoIter<C>>;

#[cfg(test)]
mod tests {
    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn tst_series() {
        log_init();

        let var = String::from("x");
        let min_pow = -10;
        let coeffs = vec![];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-11), Some(&0));
        assert_eq!(s.coeff(-10), None);

        let min_pow = -3;
        let coeffs = vec![1., 2., 3.];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow);
        assert_eq!(s.coeff(-4), Some(&0.));
        assert_eq!(s.coeff(-3), Some(&1.));
        assert_eq!(s.coeff(-2), Some(&2.));
        assert_eq!(s.coeff(-1), Some(&3.));
        assert_eq!(s.coeff(0), None);

        let min_pow = -2;
        let coeffs = vec![0., 0., 3.];
        let s = Series::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), min_pow + 2);
        assert_eq!(s.coeff(-2), Some(&0.));
        assert_eq!(s.coeff(-1), Some(&0.));
        assert_eq!(s.coeff(0), Some(&3.));
        assert_eq!(s.coeff(1), None);

        let s = Series::new(var.clone(), -2, vec![0., 0., 1.]);
        let t = Series::new(var.clone(), 0, vec![1.]);
        assert_eq!(s, t);

        let s = Series::new(var.clone(), -3, vec![0., 0., 0.]);
        let t = Series::new(var.clone(), 0, vec![]);
        assert_eq!(s, t);
    }

    #[test]
    fn tst_series_with_cutoff() {
        log_init();

        let s = Series::with_cutoff("x", -10, 1, Vec::<i32>::new());
        let t = Series::new("x", 1, vec![]);
        assert_eq!(s, t);

        let s = Series::with_cutoff("x", 0, 5, vec![1, 2, 3]);
        let t = Series::new("x", 0, vec![1, 2, 3, 0, 0]);
        assert_eq!(s, t);

        let s = Series::with_cutoff("x", 0, 2, vec![1, 2, 3]);
        let t = Series::new("x", 0, vec![1, 2]);
        assert_eq!(s, t);
    }
    #[test]
    #[should_panic]
    fn tst_bad_cutoff() {
        log_init();

        let _ = Series::with_cutoff("x", 0, -2, vec![1, 2, 3]);
    }

    #[test]
    fn tst_display() {
        log_init();

        // let s = Series::new("x", -10, vec!());
        // assert_eq!(format!("{}", s), "O(x^-10)");
        let s = Series::new("x", -3, vec![1., 0., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-3 + (-3)*x^-1 + O(x^0)");
        let s = Series::new("x", -1, vec![1., 2., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-1 + (2) + (-3)*x + O(x^2)");
    }

    #[test]
    fn tst_neg() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![-1., 0., 3.]);
        assert_eq!(res, -&s);
        assert_eq!(res, -s);
    }

    #[test]
    fn tst_add() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![2., 0., -6.]);
        assert_eq!(res, &s + &s);
        assert_eq!(res, &s + s.clone());
        assert_eq!(res, s.clone() + &s);
        assert_eq!(res, s.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        assert_eq!(s, &s + &t);
        assert_eq!(s, &t + &s);
        assert_eq!(s, &s + t.clone());
        assert_eq!(s, &t + s.clone());
        assert_eq!(s, s.clone() + &t);
        assert_eq!(s, t.clone() + &s);
        assert_eq!(s, s.clone() + t.clone());
        assert_eq!(s, t.clone() + s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -3, vec![-1., 0., 3.]);
        let res = Series::new("x", 0, vec![]);
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
        log_init();

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -3, vec![2., 0., -6.]);
        s += s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -1, vec![3., 4., 5.]);
        let t = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = s.clone();
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        s += &t;
        assert_eq!(s, res);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(s, res);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -3, vec![-1., 0., 3.]);
        let res = Series::new("x", 0, vec![]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_sub() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![]);
        assert_eq!(res, &s - &s);
        assert_eq!(res, &s - s.clone());
        assert_eq!(res, s.clone() - &s);
        assert_eq!(res, s.clone() - s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![-3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        assert_eq!(res, &s - &t);
        assert_eq!(res, &s - t.clone());
        assert_eq!(res, s.clone() - &t);
        assert_eq!(res, s - t);

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        assert_eq!(s, &s - &t);
        assert_eq!(s, &s - t.clone());
        assert_eq!(s, s.clone() - &t);
        assert_eq!(s, s.clone() - t);
    }

    #[test]
    fn tst_sub_assign() {
        log_init();

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![]);
        s -= s.clone();
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![-3., 4., 5.]);
        let res = Series::new("x", -3, vec![1., 0., 0.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = s.clone();
        let t = Series::new("x", 1, vec![3., 4., 5.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_mul() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", -6, vec![1., 0., -6.]);
        assert_eq!(res, &s * &s);
        assert_eq!(res, &s * s.clone());
        assert_eq!(res, s.clone() * &s);
        assert_eq!(res, s.clone() * s.clone());

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5., 7.]);
        let res = Series::new("x", -4, vec![3., 4., -4.]);
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);

        let s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
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
        log_init();

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s *= s.clone();
        let res = Series::new("x", -6, vec![1., 0., -6.]);
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("x", -1, vec![3., 4., 5., 7.]);
        s *= &t;
        let res = Series::new("x", -4, vec![3., 4., -4.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s *= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        s *= &t;
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        s *= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_div() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
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

        let s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        let res = Series::new("x", -6, vec![1., 14., 43.]);
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);

        let s = Series::new("x", 1, vec![1., 7., -3.]);
        let t = Series::new("x", 5, vec![]);
        let res = Series::new("x", -4, vec![]);
        assert_eq!(res, &s / &t);
        assert_eq!(res, s.clone() / &t);
        assert_eq!(res, &s / t.clone());
        assert_eq!((&res).mul_inverse(), &t / &s);
        assert_eq!(res, s / t);
    }

    #[test]
    fn tst_div_assign() {
        log_init();

        let mut s = Series::new("x", -3, vec![1., 0., -3.]);
        s /= s.clone();
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, s);

        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        let t = Series::new("x", 3, vec![1., -7., 52.]);
        s /= &t;
        let res = Series::new("x", -6, vec![1., 14., 43.]);
        assert_eq!(res, s);
        let mut s = Series::new("x", -3, vec![1., 7., -3.]);
        s /= t;
        assert_eq!(res, s);

        let mut s = Series::new("x", 1, vec![1., 7., -3.]);
        let t = Series::new("x", 5, vec![]);
        s /= &t;
        let res = Series::new("x", -4, vec![]);
        assert_eq!(res, s);
        let mut s = Series::new("x", 1, vec![1., 7., -3.]);
        s /= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_var() {
        log_init();

        let _ = Series::new(String::from("x"), -3, vec![1., 0., -3.]);
        let _ = Series::new('j', -3, vec![1., 0., -3.]);
        let _ = Series::new(8, -3, vec![1., 0., -3.]);
    }

    #[derive(Debug, Clone, PartialEq)]
    struct Mystr<'a>(&'a str);

    impl<'a> From<Mystr<'a>> for f64 {
        fn from(_s: Mystr<'a>) -> f64 {
            panic!("can't turn str to f64")
        }
    }

    #[test]
    fn tst_ln() {
        log_init();

        let s = Series::new(Mystr("x"), 0, vec![1., 7., -3.]);
        let res = Series::new(Mystr("x"), 1, vec![7., -55. / 2.]);
        assert_eq!(res, (&s).ln());
        assert_eq!(res, s.ln());

        let s = Series::new(Mystr("x"), 0, vec![4., 7., -3.]);
        let res =
            Series::new(Mystr("x"), 0, vec![4_f64.ln(), 7. / 4., -73. / 32.]);
        assert_eq!(res, (&s).ln());
        assert_eq!(res, s.ln());
    }

    #[test]
    fn tst_exp() {
        log_init();

        let s = Series::new("x", 1, vec![7., -3.]);
        let res = Series::new("x", 0, vec![1., 7., 43. / 2.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 2, vec![0.]);
        let res = Series::new("x", 0, vec![1., 0., 0.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());

        let s = Series::new("x", 0, vec![5., 11., -7.]);
        let e5 = 5_f64.exp();
        let res = Series::new("x", 0, vec![e5, e5 * 11., e5 * 107. / 2.]);
        assert_eq!(res, (&s).exp());
        assert_eq!(res, s.exp());
    }

    #[test]
    fn tst_pow() {
        log_init();

        let base = Series::new(Mystr("x"), 0, vec![1., 7., 0.]);
        let exp = Series::new(Mystr("x"), -1, vec![1., -5., 43.]);
        let e7 = 7_f64.exp();
        let res = Series::new(Mystr("x"), 0, vec![e7, -119. / 2. * e7]);
        assert_eq!(res, (&base).pow(&exp));
        assert_eq!(res, (&base).pow(exp.clone()));
        assert_eq!(res, base.clone().pow(&exp));
        assert_eq!(res, base.pow(exp));

        let base = Series::new(Mystr("x"), 0, vec![2., 7., 0.]);
        let exp = Series::new(Mystr("x"), 0, vec![3., -5., 11.]);
        // rescale result so we can use round and still get decent precision
        let rescale = Series::new(Mystr("x"), 0, vec![1e13, 0., 0., 0.]);
        let test = &rescale * &base.pow(exp);
        let ln2 = 2_f64.ln();
        let res = Series::new(
            Mystr("x"),
            0,
            vec![8., 84. - 40. * ln2, 154. + ln2 * (-332. + 100. * ln2)],
        );
        let res = rescale * res;
        assert_eq!(res.min_pow(), test.min_pow());
        assert_eq!(res.cutoff_pow(), test.cutoff_pow());
        for i in res.min_pow()..res.cutoff_pow() {
            assert_eq!(
                res.coeff(i).unwrap().round(),
                test.coeff(i).unwrap().round()
            );
        }
    }

    #[test]
    fn tst_scalar() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -2.]);
        let res = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(res, &s / 2.);
        let mut s = s;
        s /= 2.;
        assert_eq!(res, s);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -3, vec![1., 0., -2.]);
        assert_eq!(res, &s * 2.);
        let mut s = s;
        s *= 2.;
        assert_eq!(res, s);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(s, &s + 2.);
        let s = Series::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -2, vec![1. / 2., 0., 1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Series::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", 0, vec![2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Series::new("x", 0, vec![-2., 0., -1.]);
        let res = Series::new("x", 2, vec![-1.]);
        assert_eq!(res, s + 2.);

        let s = Series::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(s, &s - 2.);
        let s = Series::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", -2, vec![1. / 2., 0., -3.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Series::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Series::new("x", 0, vec![-2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Series::new("x", 0, vec![2., 0., -1.]);
        let res = Series::new("x", 2, vec![-1.]);
        assert_eq!(res, s - 2.);

        let base = Series::new(Mystr("x"), 0, vec![1., 7., 0.]);
        assert_eq!(base, (&base).pow(1.));
        assert_eq!(base, (&base).pow(&1.));
        let res = Series::new(Mystr("x"), 0, vec![1., 21., 147.]);
        assert_eq!(res, (&base).pow(3.));
        assert_eq!(res, base.pow(&3.));
    }

    #[test]
    #[should_panic]
    fn tst_bad_add() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s + t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_sub() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s - t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_mul() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s * t;
    }

    #[test]
    #[should_panic]
    fn tst_bad_div() {
        log_init();

        let s = Series::new("x", -3, vec![1., 0., -3.]);
        let t = Series::new("y", -3, vec![1., 0., -3.]);
        let _ = s / t;
    }

    #[test]
    fn tst_poly() {
        log_init();

        let var = String::from("x");
        let min_pow = -10;
        let coeffs = vec![0];
        let s = Polynomial::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), None);
        assert_eq!(s.coeff(-11), (&0));
        assert_eq!(s.coeff(-10), (&0));

        let min_pow = -3;
        let coeffs = vec![1., 2., 3.];
        let s = Polynomial::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), Some(min_pow));
        assert_eq!(s.coeff(-4), (&0.));
        assert_eq!(s.coeff(-3), (&1.));
        assert_eq!(s.coeff(-2), (&2.));
        assert_eq!(s.coeff(-1), (&3.));
        assert_eq!(s.coeff(0), (&0.));

        let min_pow = -2;
        let coeffs = vec![0., 0., 3.];
        let s = Polynomial::new(var.clone(), min_pow, coeffs);
        assert_eq!(s.min_pow(), Some(min_pow + 2));
        assert_eq!(s.coeff(-2), (&0.));
        assert_eq!(s.coeff(-1), (&0.));
        assert_eq!(s.coeff(0), (&3.));
        assert_eq!(s.coeff(1), (&0.));

        let s = Polynomial::new(var.clone(), -2, vec![0., 0., 1.]);
        let t = Polynomial::new(var.clone(), 0, vec![1.]);
        assert_eq!(s, t);

        let s = Polynomial::new(var.clone(), -3, vec![0., 0., 0.]);
        let t = Polynomial::new(var.clone(), 0, vec![]);
        assert_eq!(s, t);
    }

    #[test]
    fn tst_poly_display() {
        log_init();

        let s = Polynomial::new("x", -10, vec![0]);
        assert_eq!(format!("{}", s), "");
        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-3 + (-3)*x^-1");
        let s = Polynomial::new("x", -1, vec![1., 2., -3.]);
        assert_eq!(format!("{}", s), "(1)*x^-1 + (2) + (-3)*x");
    }

    #[test]
    fn tst_poly_neg() {
        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", -3, vec![-1., 0., 3.]);
        assert_eq!(res, -&s);
        assert_eq!(res, -s);
    }

    #[test]
    fn tst_poly_add() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", -3, vec![2., 0., -6.]);
        assert_eq!(res, &s + &s);
        assert_eq!(res, &s + s.clone());
        assert_eq!(res, s.clone() + &s);
        assert_eq!(res, s.clone() + s.clone());

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., 0., 4., 5.]);
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", 1, vec![3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., -3., 0., 3., 4., 5.]);
        assert_eq!(res, &s + &t);
        assert_eq!(res, &t + &s);
        assert_eq!(res, &s + t.clone());
        assert_eq!(res, &t + s.clone());
        assert_eq!(res, s.clone() + &t);
        assert_eq!(res, t.clone() + &s);
        assert_eq!(res, s.clone() + t.clone());
        assert_eq!(res, t.clone() + s.clone());

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -3, vec![-1., 0., 3.]);
        let res = Polynomial::new("x", 0, vec![]);
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
    fn tst_poly_add_assign() {
        log_init();

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", -3, vec![2., 0., -6.]);
        s += s.clone();
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., 0., 4., 5.]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -1, vec![3., 4., 5.]);
        let t = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", -3, vec![1., 0., 0., 4., 5.]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -1, vec![3., 4., 5.]);
        s += t;
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", 1, vec![3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., -3., 0., 3., 4., 5.]);
        s += &t;
        assert_eq!(s, res);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(s, res);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -3, vec![-1., 0., 3.]);
        let res = Polynomial::new("x", 0, vec![]);
        s += &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s += t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_poly_sub() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", 0, vec![]);
        assert_eq!(res, &s - &s);
        assert_eq!(res, &s - s.clone());
        assert_eq!(res, s.clone() - &s);
        assert_eq!(res, s.clone() - s.clone());

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![-3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., 0., -4., -5.]);
        assert_eq!(res, &s - &t);
        assert_eq!(res, &s - t.clone());
        assert_eq!(res, s.clone() - &t);
        assert_eq!(res, s - t);

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", 1, vec![3., 4., 5.]);
        let res =
            Polynomial::new("x", -3, vec![1., 0., -3., 0., -3., -4., -5.]);
        assert_eq!(res, &s - &t);
        assert_eq!(res, &s - t.clone());
        assert_eq!(res, s.clone() - &t);
        assert_eq!(res, s.clone() - t);
    }

    #[test]
    fn tst_poly_sub_assign() {
        log_init();

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", 0, vec![]);
        s -= s.clone();
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![-3., 4., 5.]);
        let res = Polynomial::new("x", -3, vec![1., 0., 0., -4., -5.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", 1, vec![3., 4., 5.]);
        let res =
            Polynomial::new("x", -3, vec![1., 0., -3., 0., -3., -4., -5.]);
        s -= &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s -= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_poly_mul() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let res = Polynomial::new("x", -6, vec![1., 0., -6., 0., 9.]);
        assert_eq!(res, &s * &s);
        assert_eq!(res, &s * s.clone());
        assert_eq!(res, s.clone() * &s);
        assert_eq!(res, s.clone() * s.clone());

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![3., 4., 5., 7.]);
        let res = Polynomial::new("x", -4, vec![3., 4., -4., -5., -15., -21.]);
        assert_eq!(res, &s * &t);
        assert_eq!(res, &t * &s);
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, &s * t.clone());
        assert_eq!(res, &t * s.clone());
        assert_eq!(res, t.clone() * s.clone());
        assert_eq!(res, s * t);

        let s = Polynomial::new("x", -3, vec![1., 7., -3.]);
        let t = Polynomial::new("x", 3, vec![1., -7., 52.]);
        let res = Polynomial::new("x", 0, vec![1., 0., 0., 385., -156.]);
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
    fn tst_poly_mul_karatsuba() {
        log_init();

        let s = Polynomial::new("x", -3, (0..=100).collect());
        let res = s.karatsuba_mul(&s, 200);
        assert_eq!(res.min_pow(), Some(-4));
        let res2 = s.karatsuba_mul(&s, 8);
        assert_eq!(res, res2);
    }

    #[test]
    fn tst_poly_mul_assign() {
        log_init();

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s *= s.clone();
        let res = Polynomial::new("x", -6, vec![1., 0., -6., 0., 9.]);
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("x", -1, vec![3., 4., 5., 7.]);
        let res = Polynomial::new("x", -4, vec![3., 4., -4., -5., -15., -21.]);
        s *= &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        s *= t;
        assert_eq!(res, s);

        let mut s = Polynomial::new("x", -3, vec![1., 7., -3.]);
        let t = Polynomial::new("x", 3, vec![1., -7., 52.]);
        let res = Polynomial::new("x", 0, vec![1., 0., 0., 385., -156.]);
        s *= &t;
        assert_eq!(res, s);
        let mut s = Polynomial::new("x", -3, vec![1., 7., -3.]);
        s *= t;
        assert_eq!(res, s);
    }

    #[test]
    fn tst_poly_var() {
        log_init();

        let _ = Polynomial::new(String::from("x"), -3, vec![1., 0., -3.]);
        let _ = Polynomial::new('j', -3, vec![1., 0., -3.]);
        let _ = Polynomial::new(8, -3, vec![1., 0., -3.]);
    }

    #[test]
    fn tst_poly_scalar() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -2.]);
        let res = Polynomial::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(res, std::ops::Div::div(&s, 2.));
        assert_eq!(res, &s / 2.);
        let mut s = s;
        s /= 2.;
        assert_eq!(res, s);

        let s = Polynomial::new("x", -3, vec![1. / 2., 0., -1.]);
        let res = Polynomial::new("x", -3, vec![1., 0., -2.]);
        assert_eq!(res, &s * 2.);
        let mut s = s;
        s *= 2.;
        assert_eq!(res, s);

        let s = Polynomial::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        let res = Polynomial::new("x", -3, vec![1. / 2., 0., -1., 2.]);
        assert_eq!(res, &s + 2.);
        let s = Polynomial::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Polynomial::new("x", -2, vec![1. / 2., 0., 1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Polynomial::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Polynomial::new("x", 0, vec![2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s + 0.);
        assert_eq!(res, s + 2.);
        let s = Polynomial::new("x", 0, vec![-2., 0., -1.]);
        let res = Polynomial::new("x", 2, vec![-1.]);
        assert_eq!(res, s + 2.);

        let s = Polynomial::new("x", -3, vec![1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        let res = Polynomial::new("x", -3, vec![1. / 2., 0., -1., -2.]);
        assert_eq!(res, &s - 2.);
        let s = Polynomial::new("x", -2, vec![1. / 2., 0., -1.]);
        let res = Polynomial::new("x", -2, vec![1. / 2., 0., -3.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Polynomial::new("x", 2, vec![1. / 2., 0., -1.]);
        let res = Polynomial::new("x", 0, vec![-2., 0., 1. / 2., 0., -1.]);
        assert_eq!(s, &s - 0.);
        assert_eq!(res, s - 2.);
        let s = Polynomial::new("x", 0, vec![2., 0., -1.]);
        let res = Polynomial::new("x", 2, vec![-1.]);
        assert_eq!(res, s - 2.);
    }

    #[test]
    #[should_panic]
    fn tst_poly_bad_add() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("y", -3, vec![1., 0., -3.]);
        let _ = s + t;
    }

    #[test]
    #[should_panic]
    fn tst_poly_bad_sub() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("y", -3, vec![1., 0., -3.]);
        let _ = s - t;
    }

    #[test]
    #[should_panic]
    fn tst_poly_bad_mul() {
        log_init();

        let s = Polynomial::new("x", -3, vec![1., 0., -3.]);
        let t = Polynomial::new("y", -3, vec![1., 0., -3.]);
        let _ = s * t;
    }
}
