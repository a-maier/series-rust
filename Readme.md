# Series

This is a crate for handling truncated Laurent series and Laurent
polynomials in a single variable about zero, i.e. expressions of the
form

```no_test
s = a_n0*x^n0 + ... + a_N*x^N + O(x^{N+1})
p = a_n0*x^n0 + ... + a_N*x^N
```

where `n0` and `N` are integers and `^` denotes exponentiation. Such
expressions can be added, subtracted, and multiplied. For Laurent series,
Some additional simple functions like division, powers of series,
natural logarithms and exponentials are also implemented.

The kinds of operations that can be performed depends on the data
type of the variable and the coefficients. For example, we usually
have to calculate the logarithm of both the leading coefficient and
the expansion variable if we take the logarithm of a Laurent
series. This crate is therefore most useful in combination with a
library providing at least basic symbolic math.

# Usage

Add this to your Cargo.toml:

```toml
[dependencies]
series = "0.13"
```

# Examples

```rust
use series::{MulInverse, SeriesIn, PolynomialIn};
use series::ops::{Ln,Exp,Pow};

// Create a new series in x, starting at order x^2 with coefficients 1, 2, 3,
// i.e. s = 1*x^2 + 2*x^3 + 3*x^4 + O(x^5).
let s = SeriesIn::new("x", 2, vec![1, 2, 3]);
println!("s = {}", s);

// The corresponding polynomial
// p = 1*x^2 + 2*x^3 + 3*x^4.
let p = PolynomialIn::new("x", 2, vec![1, 2, 3]);
assert_eq!(p, PolynomialIn::from(s));

// series with a cutoff power of 7
// s = 1*x^2 + 2*x^3 + 3*x^4 + O(x^7).
let s = SeriesIn::with_cutoff("x", 2, 7, vec![1, 2, 3]);

// To show various kinds of operations we now switch to floating-point
// coefficients

// Now s = 1 - x + O(x^5).
let s = SeriesIn::with_cutoff("x", 0, 5, vec![1., -1.]);
// Expand 1/(1-x) up to x^4.
let t = (&s).mul_inverse();
println!("1/(1-x) = {}", t);

// Series and polynomials can be added, subtracted, multiplied.
// Series can also be divided by other series.
// We can either move the arguments or use references
println!("s+t = {}", &s + &t);
println!("s-t = {}", &s - &t);
println!("s*t = {}", &s * &t);
println!("s/t = {}", &s/t);

// We can also multiply or divide each coefficient by a number
println!("s*3 = {}", &s * 3.);
println!("s/3 = {}", &s / 3.);

// More advanced operations on Laurent series in general require the
// variable type to be convertible to the coefficient type by
// implementing the From trait.
// In the examples shown here, this conversion is actually never used,
// so we can get away with a dummy implementation.
#[derive(Debug,Clone,PartialEq)]
struct Variable<'a>(&'a str);

impl<'a> From<Variable<'a>> for f64 {
    fn from(_s: Variable<'a>) -> f64 {
        panic!("Can't convert variable to f64")
    }
}

impl<'a> std::fmt::Display for Variable<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

// Now we can calculate logarithms, exponentials, and powers:
let s = SeriesIn::new(Variable("x"), 0, vec![1., -3., 5.]);
println!("exp(s) = {}", s.clone().exp());
println!("ln(s) = {}", s.clone().ln());
let t = s.clone();
println!("s^s = {}", (&s).pow(&t));
println!("s^4 = {}", s.powi(4));
```

# Multivariate Series and Polynomials

There is limited support for Laurent series and polynomials in more
than one variable via nested structures. Beyond the outermost level,
any polynomials have to be anonymous:
```
use num_traits::{One, Zero};
use series::{Polynomial, PolynomialIn};

// The polynomial 1 + x*y, where ? is an anonymous variable
let inner = Polynomial::new(1, vec![1i32]); // anonymous!
let p = PolynomialIn::new("x", 0, vec![One::one(), inner]);
assert!((&p - &p).iter().all(|(_, c)| c.is_zero()));
```
For nested series, we must use [InnerSeries],
which is a sum type of an anonymous series and its coefficient type:
```
use num_traits::{One, Zero};
use series::{InnerSeries, Series, SeriesIn};

// The nested series 1 + x*(1 + ? + O(?^3)) + O(x^2)
// where ? is an anonymous variable
let inner = InnerSeries::from(Series::with_cutoff(1..3, vec![1i32]));
let s = SeriesIn::new("x", 0, vec![One::one(), inner]);
// Note that s - s = x * O(?^3) + (x^2).
// The coefficient of x^1 is non-zero!
assert!(!(&s - &s).coeff(1).unwrap().is_zero());
// However, the coefficient of x^0 is zero:
assert!((&s - &s).coeff(0).unwrap().is_zero());
```
