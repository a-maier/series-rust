#[macro_use]
extern crate criterion;
#[macro_use]
extern crate lazy_static;

use criterion::Criterion;

use rand::prelude::*;
use rand::SeedableRng;

use series::Series;

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

impl<'a> std::ops::AddAssign<Integer> for Integer {
    fn add_assign(&mut self, rhs: Integer) {
        self.0 += rhs.0;
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
            *e = rng.gen_range(-100.0, 100.0);
        }
        array
    };

    static ref RAN_INT: Vec<Integer> = {
        let mut rng = rand_pcg::Pcg64::seed_from_u64(0);
        let mut array  = Vec::new();
        let mut digits = [0u8; MAX_DIGITS];
        for _ in 0..MAX_ELEMENTS {
            let ndigits = rng.gen_range(1, MAX_DIGITS);
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
        "1 f64 coefficient",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_10(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..10].to_owned());
    c.bench_function(
        "10 f64 coefficients",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_100(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..100].to_owned());
    c.bench_function(
        "100 f64 coefficients",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_f64_1000(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_F64[..1000].to_owned());
    c.bench_function(
        "1000 f64 coefficients",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_1(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..1].to_owned());
    c.bench_function(
        "1 integer coefficient",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_10(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..10].to_owned());
    c.bench_function(
        "10 integer coefficients",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_100(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..100].to_owned());
    c.bench_function(
        "100 integer coefficients",
        move |b| {
            let s = s.clone();
            b.iter(|| &s * &s)
        }
    );
}

fn mul_int_1000(c: &mut Criterion) {
    let s = Series::new("x", -2, RAN_INT[..1000].to_owned());
    c.bench(
        "",
        criterion::Benchmark::new(
            "1000 integer coefficients",
            move |b| {
                let s = s.clone();
                b.iter(|| &s * &s)
            }
        ).sample_size(20)
    );
}


criterion_group!(
    benches,
    mul_f64_1, mul_f64_10, mul_f64_100, mul_f64_1000,
    mul_int_1, mul_int_10, mul_int_100, mul_int_1000,
);
criterion_main!(benches);
