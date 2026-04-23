# Version 0.14.0

A major rewrite with many backwards-incompatible changes.

- Added macros `var!` and `once_var!` for defining light-weight
  variables with name set at compile time (`var!`) or once at runtime
  (`once_var!`).
- Added `Laurent` as a sum type of series and polynomials. `Laurent`
  replaces `InnerSeries`, which was removed.
- Redesigned `Polynomial` as a sum type of constants and non-constant
  polynomials. This type replaces the old `PolynomialIn`. Polynomials
  with anonymous variables (`Polynomial` in versions 0.12 and 0.13)
  are no longer supported.
- Renamed
  `SeriesIn` -> `Series`
  `Series` -> `AnonSeries`
- Multiplying series by zero now leads to a panic instead of an
  incorrect result.
- Rewrote `Display` implementations to support series and polynomials
  whose coefficients implement the new `SplitSign` and
  `NeedsCoeffBracket` traits. This includes integers and floats, as
  well as nested series and polynomials.
- Added `O!` macro for representing series without coefficients.
- Extended support for arithmetic operations, e.g. combining series
  and polynomials.
- Added basic `FromStr` implementations for series and polynomials via
  the `parse` feature.
- Added optional support for arbitrary-precision integer and rational
  coefficients via [rug](https://crates.io/crates/rug).
- Added `map` method to replace series and polynomial coefficients.
- Added `truncate_at` method to truncate series at a lower expansion order.
- Added `shift_pow` method to shift all powers in a series.
- Added `replace_var` method to replace the variable of a series or
  polynomial.
- Updated the crate to the 2024 Rust edition.

# Version 0.13.0

- Added `InnerSeries` type to support nested series.
- Added `for_each` and `apply_at` methods for modifying coefficients.
- Implemented `Default` for `Polynomial`.
- Fixed lifetimes of `var` methods in `PolynomialSliceIn`, `SeriesSliceIn`.
- Implemented `AddAssign` for references to the coefficient type.

# Version 0.12.0

- Added Conversions
  `SeriesIn` -> `Series`
  `PolynomialIn` -> `Polynomial`
  `SeriesSliceIn` -> `SeriesSlice`
  `PolynomialSliceIn` -> `PolynomialSlice`

# Version 0.11.0

- Fixed the behaviour of `powi` for series without coefficients,
  i.e. series of the form O(x^n).
- Replaced the custom `Pow` trait by the
  [`num_traits`](https://crates.io/crates/num-traits) equivalent.

# Version 0.10.0

- `Series`, `SeriesSlice`, `Polynomial`, `PolynomialSlice` now
  represent Laurent series and polynomials in an anonymous expansion
  variable.
- Relaxed trait bounds on coefficient types.

# Version 0.9.0

- Renamed all structs:
  `Series` -> `SeriesIn`
  `Polynomial` -> `PolynomialIn`
  `SeriesSlice` -> `SeriesSliceIn`
  `PolynomialSlice` -> `PolynomialSliceIn`
  `SeriesParts` -> `SeriesInParts`
  `PolynomialParts` -> `PolynomialInParts`
- Fixed lifetime in `coeff` method return values.
- Fixed off-by-one error in the `PolynomialSlice` `coeff` method.

# Version 0.8.0

- Update to Rust 2021.
- Fixed broken `serde` support.
- Updated dependencies.

# Version 0.7.1

- Fixed typo in Readme.

# Version 0.7.0

- Implemented Laurent polynomials.
- Added `SeriesSlice` and `PolynomialSlice` as immutable views into
  series and polynomials.
- Implemented `Iter` and `IntoIter` for iteration over series and
  polynomial coefficients.
- Added multiplication benchmarks.

# Version 0.6.0

- Update to Rust 2018 (1.35.0)

# Version 0.5.0

- Implemented `powi` as a better method to compute a series to an
  integer power.
- Added support for destructuring series via the new `SeriesData` struct.
- Fixed regression: `pow` method once again works on series references.

# Version 0.4.0

- Implemented exponentiation with an exponent of the same type as the
  series coefficients.
- Renamed `max_pow` method to `cutoff_pow`.
- Added a constructor with explicit series cutoff power.

# Version 0.3.0

- Implemented multiplication and division by objects with the same type
  as the series coefficients.
- Implemented `Ord` and `PartialOrd` traits for `Series`.

# Version 0.2.0

- Fixed compilation issues with rustc 1.28.0.
- Implemented all arithmetic operations for all combinations of
  references and values.
- Updated traits for logarithm and exponential function and added trait
  for powers.
