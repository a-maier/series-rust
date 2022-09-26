#[macro_use]
extern crate criterion;
#[macro_use]
extern crate lazy_static;

use criterion::Criterion;

use rand::prelude::*;
use rand::SeedableRng;

use series::{Series, Polynomial, KaratsubaMul};

const MAX_ELEMENTS: usize = 2000;
const MAX_DIGITS: usize = 20;

#[derive(Clone,Debug,PartialOrd,PartialEq,Ord,Eq)]
struct Integer(rug::Integer);

impl Integer {
    pub fn from_digits<T>(digits: &[T], order: rug::integer::Order) -> Self
    where T: rug::integer::UnsignedPrimitive
    {
        Integer(rug::Integer::from_digits(digits, order))
    }
}

impl std::ops::Neg for Integer {
    type Output = Self;
    fn neg(self) -> Self {
        Integer(-self.0)
    }
}

impl<'a> std::ops::Neg for &'a Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        Integer((-&self.0).into())
    }
}

impl std::ops::AddAssign<Integer> for Integer {
    fn add_assign(&mut self, rhs: Integer) {
        self.0 += rhs.0;
    }
}

impl<'a> std::ops::AddAssign<&'a Integer> for Integer {
    fn add_assign(&mut self, rhs: &'a Integer) {
        self.0 += &rhs.0;
    }
}

impl std::ops::Add<Integer> for Integer {
    type Output = Integer;

    fn add(self, rhs: Integer) -> Integer {
        Integer(rug::Integer::from(self.0 + rhs.0))
    }
}

impl<'a> std::ops::MulAssign<&'a Integer> for Integer {
    fn mul_assign(&mut self, rhs: &'a Integer) {
        self.0 *= &rhs.0;
    }
}

impl<'a, 'b> std::ops::Mul<&'b Integer> for &'a Integer {
    type Output = Integer;

    fn mul(self, rhs: &'b Integer) -> Integer {
        Integer(rug::Integer::from(&self.0 * &rhs.0))
    }
}

impl std::ops::Mul<Integer> for Integer {
    type Output = Integer;

    fn mul(self, rhs: Integer) -> Integer {
        Integer(rug::Integer::from(self.0 * rhs.0))
    }
}

impl num_traits::Zero for Integer {
    fn zero() -> Self {
        Integer(rug::Integer::new())
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl num_traits::One for Integer {
    fn one() -> Self {
        Integer(rug::Integer::from(1))
    }

    fn is_one(&self) -> bool {
        self.0 == 1
    }
}

impl std::convert::From<i32> for Integer {
    fn from(i: i32) -> Integer {
        Integer(rug::Integer::from(i))
    }
}

impl std::fmt::Display for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}


lazy_static! {
    static ref RAN_F64: [f64; MAX_ELEMENTS] = {
        let mut rng = rand_pcg::Pcg64::seed_from_u64(0);
        let mut array = [0.; MAX_ELEMENTS];
        for e in array.iter_mut() {
            *e = rng.gen_range(-100.0..=100.0);
        }
        array
    };

    static ref RAN_INT: Vec<Integer> = {
        let mut rng = rand_pcg::Pcg64::seed_from_u64(0);
        let mut array  = Vec::new();
        let mut digits = [0u8; MAX_DIGITS];
        for _ in 0..MAX_ELEMENTS {
            let ndigits = rng.gen_range(1..=MAX_DIGITS);
            rng.fill(&mut digits[..ndigits]);
            array.push(
                Integer::from_digits(
                    &digits[..ndigits],
                    rug::integer::Order::Lsf
                )
            );
        }
        array
    };

}

fn mul_f64_1(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..1].to_owned());
    c.bench_function(
        "multiply series with 1 f64 coefficient",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_10(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..10].to_owned());
    c.bench_function(
        "multiply series with 10 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_100(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..100].to_owned());
    c.bench_function(
        "multiply series with 100 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_1000(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "multiply series with 1000 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_1(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..1].to_owned());
    c.bench_function(
        "multiply series with 1 integer coefficient",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_10(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..10].to_owned());
    c.bench_function(
        "multiply series with 10 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_100(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..100].to_owned());
    c.bench_function(
        "multiply series with 100 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_1000(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..1000].to_owned());
    let mut group = c.benchmark_group("dummy name");
    group.sample_size(20);
    group.bench_function(
        "multiply series with 1000 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_f64_1(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..1].to_owned());
    c.bench_function(
        "multiply polynomials with 1 f64 coefficient",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_f64_10(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..10].to_owned());
    c.bench_function(
        "multiply polynomials with 10 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_f64_100(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..100].to_owned());
    c.bench_function(
        "multiply polynomials with 100 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_f64_1000(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "multiply polynomials with 1000 f64 coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_int_1(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..1].to_owned());
    c.bench_function(
        "multiply polynomials with 1 integer coefficient",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_int_10(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..10].to_owned());
    c.bench_function(
        "multiply polynomials with 10 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_int_100(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..100].to_owned());
    c.bench_function(
        "multiply polynomials with 100 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_int_1000(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..1000].to_owned());
    let mut group = c.benchmark_group("dummy name");
    group.sample_size(20);
    group.bench_function(
        "multiply polynomials with 1000 integer coefficients",
        move |b| {
            b.iter(|| &s * &s)
        }
    );
}

fn mul_poly_int_karatsuba_4(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..1000].to_owned());
    let mut group = c.benchmark_group("dummy name");
    group.sample_size(20);
    group.bench_function(
        "multiply int polynomials with Karatsuba threshold 4",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 4))
        }
    );
}

fn mul_poly_int_karatsuba_8(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..1000].to_owned());
    let mut group = c.benchmark_group("dummy name");
    group.sample_size(20);
    group.bench_function(
        "multiply int polynomials with Karatsuba threshold 8",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 8))
        }
    );
}

fn mul_poly_int_karatsuba_16(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_INT[..1000].to_owned());
    let mut group = c.benchmark_group("dummy name");
    group.sample_size(20);
    group.bench_function(
        "multiply int polynomials with Karatsuba threshold 16",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 16))
        }
    );
}

fn mul_poly_f64_karatsuba_4(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "multiply f64 polynomials with Karatsuba threshold 4",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 4))
        }
    );
}

fn mul_poly_f64_karatsuba_8(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "multiply f64 polynomials with Karatsuba threshold 8",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 8))
        }
    );
}

fn mul_poly_f64_karatsuba_16(c: &mut Criterion) {
    let s = Polynomial::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "multiply f64 polynomials with Karatsuba threshold 16",
        move |b| {
            b.iter(|| s.karatsuba_mul(&s, 16))
        }
    );
}

criterion_group!(
    benches,
    mul_f64_1, mul_f64_10, mul_f64_100, mul_f64_1000,
    mul_int_1, mul_int_10, mul_int_100, mul_int_1000,
    mul_poly_f64_1, mul_poly_f64_10, mul_poly_f64_100, mul_poly_f64_1000,
    mul_poly_int_1, mul_poly_int_10, mul_poly_int_100, mul_poly_int_1000,
    mul_poly_f64_karatsuba_4, mul_poly_f64_karatsuba_8, mul_poly_f64_karatsuba_16,
    mul_poly_int_karatsuba_4, mul_poly_int_karatsuba_8, mul_poly_int_karatsuba_16,
);
criterion_main!(benches);
