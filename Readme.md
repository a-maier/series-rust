This is a crate for handling truncated Laurent series in a single
variable about zero, i.e. expressions of the form

`s = a_n0*x^n0 + ... + a_N*x^N + O(x^{N+1}),`

where `n0` and `N` are integers and `^` denotes exponentiation. Such
series can be added, subtracted, multiplied, and divided. Some
simple functions like powers of series, natural logarithms and
exponentials are also implemented.

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
series = "0.1.0"
```

and this to your crate root:
```rust
extern crate series;
```

# Examples

```rust
use series::Series;
use series::ops::{Ln,Exp};

// Create a new series in x, starting at order x^2 with coefficients 1, 2, 3,
// i.e. s = 1*x^2 + 2*x^3 + 3*x^4 + O(x^5).
let s = Series::new("x", 2, vec!(1, 2, 3));
println!("s = {}", s);

// To show various kinds of operations we now switch to floating-point
// coefficients

// Now s = 1 - x + O(x^5).
let s = Series::new("x", 0, vec!(1., -1., 0., 0.));
// Expand 1/(1-x) up to x^4.
let t = s.inverse();
println!("1/(1-x) = {}", t);

// Series can be added, subtracted, multiplied, divided.
// We can either move the arguments or use references
println!("s+t = {}", &s + &t);
println!("s-t = {}", &s - &t);
println!("s*t = {}", &s * &t);
println!("s/t = {}", s/t);

// More advanced operations in general require the variable type to be
// convertible to the coefficient type by implementing the From trait.
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
let s = Series::new(Variable("x"), 0, vec!(1., -3., 5.));
println!("exp(s) = {}", s.clone().exp());
println!("ln(s) = {}", s.clone().ln());
let t = s.clone();
println!("s^s = {}", s.pow(t));
```
