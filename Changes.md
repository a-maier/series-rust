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
