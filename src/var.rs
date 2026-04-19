use std::{fmt::Display, marker::PhantomData};

/// Error parsing a variable
#[derive(Debug, Default)]
pub struct VarParseError<Var>(PhantomData<Var>);

impl<Var> VarParseError<Var> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

impl<Var: Default + Display> Display for VarParseError<Var> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "String is not {}", Var::default())
    }
}

/// Define a compile-time variable
///
/// This is a convenience macro for defining variables for
/// [Series](crate::Series) or [Polynomial](crate::Polynomial).
/// Effectively, `var!(VarName)` defines a unit-like struct VarName`
/// with an associated string constant `var_name`. [Display] and
/// [FromStr](std::str::FromStr) are implemented matching the value of
/// this string constant.
///
/// # Example
/// ```
/// # use series::{Polynomial, var};
/// var!(X);
/// let p: Polynomial<X, i32> = "x".parse().unwrap();
/// assert_eq!(p.to_string(), "x");
///
/// // with custom visibility
/// var!{pub(crate) Y};
/// let p: Polynomial<Y, i32> = "y".parse().unwrap();
/// assert_eq!(p.to_string(), "y");
/// ```
#[macro_export]
macro_rules! var {
    ($var:ident) => {
        var!(pub(self) $var);
    };
    ($v:vis $var:ident) => {
        $crate::paste::paste! {
            #[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
            $v struct [< $var:camel >];

            impl [< $var:camel >] {
                $v const fn name() -> &'static str {
                    stringify!([< $var:snake:lower >])
                }
            }


            impl std::fmt::Display for [< $var:camel >] {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    std::fmt::Display::fmt(Self::name(), f)
                }
            }

            impl std::str::FromStr for [< $var:camel >] {
                type Err = $crate::var::VarParseError<Self>;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    if s == Self::name() {
                        Ok([< $var:camel >])
                    } else {
                        Err($crate::var::VarParseError::new())
                    }
                }
            }
        }
    };

}

/// Define a variable with name known at run time.
///
/// This is an alternative to [var] for cases where the name of the
/// variable is only known at run time. It can and must be set exactly
/// once.
///
/// # Examples
/// ```
/// # use series::{Polynomial, once_var};
/// once_var!(X);
/// // We need to set a name before parsing or printing `X` for the first time
/// X::set_name("x").unwrap();
/// let p: Polynomial<X, i32> = "x".parse().unwrap();
/// assert_eq!(p.to_string(), "x");
/// ```
///
/// Since the variable name is only known at run time it is quite likely
/// that its lifetime is not `'static`. To extend it we can use a
/// [Box] or [String] and leak the allocation:
///
/// ```
/// # use series::{Polynomial, once_var};
/// once_var!(Z);
/// let name = "z".to_string();
/// Z::set_name(name.leak()).unwrap();
/// ```
///
/// Variables can be defined with custom visibility:
/// ```
/// # use series::{Polynomial, once_var};
/// once_var!{pub(crate) Y};
/// Y::set_name("y").unwrap();
/// let p: Polynomial<Y, i32> = "y".parse().unwrap();
/// assert_eq!(p.to_string(), "y");
/// ```
#[macro_export]
macro_rules! once_var {
    ($var:ident) => {
        once_var!(pub(self) $var);
    };
    ($v:vis $var:ident) => {
        $crate::paste::paste! {
            #[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
            $v struct [< $var:camel >];

            impl [< $var:camel >] {
                $v fn set_name(name: &'static str) -> Result<(), &'static str> {
                    [< $var:snake:upper _STR >].set(name)
                }

                $v fn name() -> Option<&'static str> {
                    [< $var:snake:upper _STR >].get().map(|v| *v)
                }
            }

            static [< $var:snake:upper _STR >]: std::sync::OnceLock<&str> = std::sync::OnceLock::new();

            impl std::fmt::Display for [< $var:camel >] {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    let name = Self::name()
                        .expect("Variable name not set. Use `set_name` before using the variable.");
                    std::fmt::Display::fmt(name, f)
                }
            }

            impl std::str::FromStr for [< $var:camel >] {
                type Err = $crate::var::VarParseError<Self>;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    let name = [< $var:camel >]::name()
                        .expect("Variable name not set. Use `set_name` before using the variable.");
                    if s == name {
                        Ok([< $var:camel >])
                    } else {
                        Err($crate::var::VarParseError::new())
                    }
                }
            }
        }
    };
}
