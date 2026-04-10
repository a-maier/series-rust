//! Suppport for [rug](https://crates.io/crates/rug) types as coefficients
#[cfg(feature = "rug-integer")]
pub mod integer;
#[cfg(feature = "rug-rational")]
pub mod rational;

#[cfg(feature = "rug-integer")]
pub use integer::Integer;

#[cfg(feature = "rug-rational")]
pub use rational::Rational;
