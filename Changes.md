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
