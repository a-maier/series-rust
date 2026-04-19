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
use series::{Laurent, MulInverse, Series, Polynomial, var};

fn main() {

   // Define a variable `X`, rendered as "x"
   //
   // This defines a struct `X` with name `X::name() == "x"`.
   var!(X);

   // Create a new series in x
   let s: Series<X, i32> = "x^2 + 2*x^3 + 3*x^4 + O(x^5)".parse().unwrap();
   // We can also use the more efficient direct construction
   let t = Series::new(X, 2, vec![1, 2, 3]);
   assert_eq!(s, t);

   // There is also a constructor with an explicit cutoff power
   let s = Series::with_cutoff(X, 2..7, vec![1, 2, 3]);
   assert_eq!(s, "x^2 + 2*x^3 + 3*x^4 + O(x^7)".parse().unwrap());

   // The corresponding polynomial
   let p: Polynomial<X, i32> = "x^2 + 2*x^3 + 3*x^4".parse().unwrap();
   assert_eq!(p, Polynomial::from(s));

   // Series and Polynomials are generic in their variable and coefficient types
   let p: Polynomial<String, i32> = "x^2 + 2*x^3 + 3*x^4".parse().unwrap();

   // To show various kinds of operations we now switch to floating-point
   // coefficients

   let s: Series<X, f64> = "1. - x + O(x^5)".parse().unwrap();
   // Expand 1 / (1 - x) up to x^4.
   let t = (&s).mul_inverse();
   println!("1/(1-x) = {t}");

   // Series and polynomials can be added, subtracted, multiplied.
   // Series can also be divided by other series.
   // We can either move the arguments or use references
   println!("s + t = {}", &s + &t);
   println!("s - t = {}", &s - &t);
   println!("s * t = {}", &s * &t);
   println!("s / t = {}", &s / t);

   // We can also multiply or divide each coefficient by a number
   println!("s * 3 = {}", &s * 3.);
   println!("s / 3 = {}", &s / 3.);
   // Polynomials can be multiplied by 0
   assert_eq!(&p * 0, Polynomial::zero());
   // But for series, multiplication by 0 is not defined,
   // as the result is not a series anymore
   assert!(std::panic::catch_unwind(|| &s * 0.).is_err());
   // We can instead use the Laurent sum type, which can be either
   // a series or a polynomial
   let l = Laurent::from(s);
   assert_eq!(&l * 0., Laurent::zero());
}
```

# Exponentials, logarithms, and powers

More advanced operations on Laurent series in general require the
variable type to be convertible to the coefficient type by
implementing the
[From](https://doc.rust-lang.org/std/convert/trait.From.html) trait:

```
use series::{Series, var};
use series::ops::{Ln, Exp, Pow};

var!(X);

fn main() {
   // In the examples shown here, this conversion is actually never used,
   // so we can get away with a dummy implementation.
   impl From<X> for f64 {
       fn from(_: X) -> f64 {
           panic!("Can't convert variable `x` to f64")
       }
   }

   // Now we can calculate logarithms, exponentials, and powers:
   let s = Series::new(X, 0, vec![1., -3., 5.]);
   println!("exp(s) = {}", s.clone().exp());
   println!("ln(s) = {}", s.clone().ln());
   let t = s.clone();
   println!("s^s = {}", (&s).pow(&t));
   println!("s^4 = {}", s.powi(4));
}
```

# Multivariate series and polynomials

Multivariate polynomials can be created via nesting:

```
use series::{Polynomial, var};

var!(X);
var!(Y);

fn main() {
   let p: Polynomial<X, Polynomial<Y, i32>> =
      "+ 1 + 2*y + 3*y^2
       + (4 + 2*y + 3*y^2)*x
       + (2 + 4*y + 6*y^2)*x^2".parse().unwrap();
   let q = Polynomial::new(Y, 0, vec![1, 2, 3]);
   let q = Polynomial::new(X, 0, vec![q.clone(), &q + 3,  &q * 2]);
   assert_eq!(p, q);
}
```

Series are not supported as coefficients, since they can never
be identically zero. A way around this limitation is to use the
[Laurent] struct:

```
use series::{Laurent, Polynomial, Series, var};

var!(X);
var!(Y);

fn main() {
   let s = Series::new(Y, 0, vec![1, 2, 3]);
   let s = Laurent::from(s);
   let p = Polynomial::new(X, 0, vec![s.clone(), &s + 3,  &s * 2]);
}
```
