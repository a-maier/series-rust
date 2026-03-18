use std::{fmt::Display, marker::PhantomData};

/// Error parsing a variable
#[derive(Debug)]
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

#[macro_export]
macro_rules! var {
    ($var:ident) => {
        $crate::paste::paste! {
            #[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
            struct [< $var:camel >] { }

            const [< $var:snake:upper >]: [< $var:camel >] = [< $var:camel >]{};

            const [< $var:snake:upper _STR >]: &str = stringify!([< $var:snake:lower >]);

            impl std::fmt::Display for [< $var:camel >] {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    std::fmt::Display::fmt([< $var:snake:upper _STR >], f)
                }
            }

            impl std::str::FromStr for [< $var:camel >] {
                type Err = $crate::var::VarParseError<Self>;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    if s == [< $var:snake:upper _STR >] {
                        Ok([< $var:snake:upper >])
                    } else {
                        Err($crate::var::VarParseError::new())
                    }
                }
            }
        }
    };
    ($v:vis $var:ident) => {
        $crate::paste::paste! {
            #[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
            $v struct [< $var:camel >] { }

            $v const [< $var:snake:upper >]: [< $var:camel >] = [< $var:camel >]{};

            $v const [< $var:snake:upper _STR >]: &str = stringify!([< $var:snake:lower >]);

            impl std::fmt::Display for [< $var:camel >] {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    std::fmt::Display::fmt([< $var:snake:upper _STR >], f)
                }
            }

            impl std::str::FromStr for [< $var:camel >] {
                type Err = $crate::var::VarParseError<Self>;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    if s == [< $var:snake:upper _STR >] {
                        Ok([< $var:snake:upper >])
                    } else {
                        Err($crate::var::VarParseError::new())
                    }
                }
            }
        }
    };

}
