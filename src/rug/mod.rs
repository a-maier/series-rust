//! Suppport for [rug](https://crates.io/crates/rug) types as coefficients
#[cfg(feature = "rug-integer")]
pub mod integer;
#[cfg(feature = "rug-rational")]
pub mod rational;
