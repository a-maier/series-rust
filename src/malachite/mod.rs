//! Suppport for [malachite](https://crates.io/crates/malachite) types as coefficients
#[cfg(feature = "malachite-integer")]
pub mod integer;
#[cfg(all(feature = "malachite-integer", feature = "parse"))]
pub mod parse_integer;
#[cfg(all(feature = "malachite-rational", feature = "parse"))]
pub mod parse_rational;
#[cfg(feature = "malachite-rational")]
pub mod rational;

#[cfg(feature = "malachite-integer")]
pub use integer::Integer;

#[cfg(feature = "malachite-rational")]
pub use rational::Rational;
