use std::{
    fmt::Debug,
    ops::{AddAssign, Mul, Neg},
    str::FromStr,
};

use derive_more::Display;
use winnow::{
    ModalResult, Parser,
    ascii::{dec_int, dec_uint, float, multispace0},
    combinator::{alt, delimited, opt, preceded, repeat, separated_pair, trace},
    error::ContextError,
    token::{any, take_while},
};

use crate::{Coeff, Laurent, Polynomial, Series, Sign};

/// Error parsing a polynomial
#[derive(Debug, Display)]
pub struct ParseError(String); // TODO: make this more useful

impl std::error::Error for ParseError {}

impl<Var, C> FromStr for Polynomial<Var, C>
where
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
    Var: FromStr + PartialEq,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        poly.parse(s).map_err(|e| ParseError(e.to_string()))
    }
}

impl<Var, C> FromStr for Series<Var, C>
where
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
    Var: Clone + Debug + FromStr + PartialEq,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        series.parse(s).map_err(|e| ParseError(e.to_string()))
    }
}

impl<Var, C> FromStr for Laurent<Var, C>
where
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
    Var: Clone + Debug + FromStr + PartialEq,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        laurent.parse(s).map_err(|e| ParseError(e.to_string()))
    }
}

fn poly<Var, C>(input: &mut &str) -> ModalResult<Polynomial<Var, C>>
where
    Var: FromStr + PartialEq,
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
{
    (opt_signed_monomial, repeat(0.., signed_monomial))
        .verify_map(|(first, mut rest): (_, Vec<_>)| {
            rest.push(first);
            poly_from_monomials(rest).ok()
        })
        .parse_next(input)
}

fn series<Var, C>(input: &mut &str) -> ModalResult<Series<Var, C>>
where
    Var: Debug + Clone + FromStr + PartialEq,
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
{
    alt((
        separated_pair(poly, plus, cutoff)
            .map(|(p, (var, pow))| p.cutoff_at(&var, pow)),
        preceded(opt(plus), cutoff)
            .map(|(var, pow)| Series::new(var, pow, vec![])),
    ))
        .parse_next(input)
}

fn laurent<Var, C>(input: &mut &str) -> ModalResult<Laurent<Var, C>>
where
    Var: Debug + Clone + FromStr + PartialEq,
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
{
    alt((
        (poly, opt(preceded(plus, cutoff)))
            .map(|(p, c)| if let Some((var, pow)) = c {
                p.cutoff_at(&var, pow).into()
            } else {
                p.into()
            }),
        preceded(opt(plus), cutoff)
            .map(|(var, pow)| Series::new(var, pow, vec![]).into()),
    ))
        .parse_next(input)
}

#[derive(Debug)]
struct VarMismatch {}

fn poly_from_monomials<Var: PartialEq, C: AddAssign + Coeff>(
    monomials: Vec<Monomial<Var, C>>,
) -> Result<Polynomial<Var, C>, VarMismatch> {
    assert!(!monomials.is_empty());
    let res = if is_const(&monomials)? {
        debug_assert!(monomials.iter().all(|m| m.pow == 0));
        let c = monomials
            .into_iter()
            .map(|m| m.c)
            .reduce(std::ops::Add::add)
            .unwrap();
        Polynomial::from_const(c)
    } else {
        let pows = monomials.iter().map(|m| m.pow);
        let min_pow = pows.clone().min().unwrap();
        let max_pow = pows.max().unwrap();
        let mut c = Vec::from_iter((min_pow..=max_pow).map(|_| C::zero()));
        let mut var = None;
        for m in monomials {
            let pos = (m.pow - min_pow) as usize;
            c[pos] += m.c;
            if let Some(v) = m.var {
                var = Some(v);
            }
        }
        Polynomial::new(var.unwrap(), min_pow, c)
    };
    Ok(res)
}

fn is_const<Var: PartialEq, C>(
    monomials: &[Monomial<Var, C>],
) -> Result<bool, VarMismatch> {
    let mut vars = monomials.iter().flat_map(|m| m.var.as_ref());
    let Some(var) = vars.next() else {
        return Ok(true);
    };
    if vars.all(|v| v == var) {
        Ok(false)
    } else {
        Err(VarMismatch {})
    }
}

fn opt_signed_monomial<
    Var: FromStr,
    C: Coeff + Neg<Output = C> + ParseCoeff,
>(
    input: &mut &str,
) -> ModalResult<Monomial<Var, C>> {
    (opt(sign).map(|s| s.unwrap_or(Sign::Plus)), monomial)
        .map(|(s, m)| match s {
            Sign::Plus => m,
            Sign::Minus => -m,
        })
        .parse_next(input)
}

fn signed_monomial<
    Var: FromStr,
    C: Coeff + Neg<Output = C> + ParseCoeff,
>(
    input: &mut &str,
) -> ModalResult<Monomial<Var, C>> {
    (sign, monomial)
        .map(|(s, m)| match s {
            Sign::Plus => m,
            Sign::Minus => -m,
        })
        .parse_next(input)
}

fn monomial<Var: FromStr, C: Coeff + ParseCoeff>(
    input: &mut &str,
) -> ModalResult<Monomial<Var, C>> {
    trace("monomial", alt((coeff_times_var_pow, var_pow_as_monomial))).parse_next(input)
}

fn coeff_times_var_pow<Var: FromStr, C: Coeff + ParseCoeff>(
    input: &mut &str,
) -> ModalResult<Monomial<Var, C>> {
    (coeff, opt((times_or_div, var_pow_as_monomial)))
        .map(|(c, op_var_pow)| match op_var_pow {
            None => Monomial {
                c,
                var: None,
                pow: 0,
            },
            Some((MulOp::Times, m)) => m * c,
            Some((MulOp::Div, mut m)) => {
                m.pow *= -1;
                m * c
            }
        })
        .parse_next(input)
}

fn coeff<C: ParseCoeff>(input: &mut &str) -> ModalResult<C> {
    alt((
        delimited(open_bracket, coeff_bracket, closing_bracket),
        coeff_no_bracket,
    ))
    .parse_next(input)
}

fn coeff_bracket<C: ParseCoeff>(input: &mut &str) -> ModalResult<C> {
    let res = delimited(open_bracket, coeff_bracket, closing_bracket)
        .parse_next(input);
    if res.is_ok() {
        return res
    }
    C::parse_coeff(input, true)
        // TODO: how to create a proper error?
        .map_err(|_| winnow::error::ErrMode::Backtrack(ContextError::new()))
}

fn coeff_no_bracket<C: ParseCoeff>(input: &mut &str) -> ModalResult<C> {
    C::parse_coeff(input, false)
        // TODO: how to create a proper error?
        .map_err(|_| winnow::error::ErrMode::Backtrack(ContextError::new()))
}

fn var_pow_as_monomial<Var: FromStr, C: Coeff>(
    input: &mut &str,
) -> ModalResult<Monomial<Var, C>> {
    var_pow
        .map(|(var, pow)| Monomial {
                c: C::one(),
                var: Some(var),
                pow,
        })
        .parse_next(input)
}

fn var<Var: FromStr>(input: &mut &str) -> ModalResult<Var> {
    (
        any.verify(|c: &char| c.is_alphabetic()),
        take_while(0.., |c: char| c.is_alphanumeric() || c == '_'),
    )
        .take()
        .parse_to()
        .parse_next(input)
}

fn num_pow(input: &mut &str) -> ModalResult<isize> {
    preceded(pow, dec_int).parse_next(input)
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Monomial<Var, C> {
    c: C,
    var: Option<Var>,
    pow: isize,
}

impl<C: Neg<Output = C>, Var> Neg for Monomial<Var, C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let Self { c, var, pow } = self;
        Self { c: -c, var, pow }
    }
}

impl<C: Mul<Output = C>, Var> Mul<C> for Monomial<Var, C> {
    type Output = Self;

    fn mul(self, rhs: C) -> Self::Output {
        let Self { c, var, pow } = self;
        Self {
            c: c * rhs,
            var,
            pow,
        }
    }
}

pub trait ParseCoeff
where
    Self: Sized,
{
    type Error;

    fn parse_coeff(
        input: &mut &str,
        inside_bracket: bool
    ) -> Result<Self, Self::Error>;
}

fn open_bracket(input: &mut &str) -> ModalResult<()> {
    delimited(multispace0, '(', multispace0)
        .void()
        .parse_next(input)
}

fn closing_bracket(input: &mut &str) -> ModalResult<()> {
    delimited(multispace0, ')', multispace0)
        .void()
        .parse_next(input)
}

fn sign(input: &mut &str) -> ModalResult<Sign> {
    alt((plus, minus)).parse_next(input)
}

fn plus(input: &mut &str) -> ModalResult<Sign> {
    delimited(multispace0, '+', multispace0)
        .value(Sign::Plus)
        .parse_next(input)
}

fn minus(input: &mut &str) -> ModalResult<Sign> {
    delimited(multispace0, '-', multispace0)
        .value(Sign::Minus)
        .parse_next(input)
}

fn times_or_div(input: &mut &str) -> ModalResult<MulOp> {
    alt((times, div)).parse_next(input)
}

fn times(input: &mut &str) -> ModalResult<MulOp> {
    delimited(multispace0, '*', multispace0)
        .value(MulOp::Times)
        .parse_next(input)
}

fn div(input: &mut &str) -> ModalResult<MulOp> {
    delimited(multispace0, '/', multispace0)
        .value(MulOp::Div)
        .parse_next(input)
}

fn pow(input: &mut &str) -> ModalResult<()> {
    delimited(multispace0, alt(('^'.void(), "**".void())), multispace0)
        .parse_next(input)
}

fn var_pow<Var: FromStr>(
    input: &mut &str,
) -> ModalResult<(Var, isize)> {
    (var, opt(num_pow))
        .map(|(v, p)| (v, p.unwrap_or(1)))
        .parse_next(input)
}

fn cutoff<Var: FromStr>(input: &mut &str) -> ModalResult<(Var, isize)> {
    delimited(("O(", multispace0), var_pow, closing_bracket)
        .parse_next(input)
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum MulOp {
    Times,
    Div,
}

macro_rules! impl_parse_coeff_int {
    ($($t:ty), *) => {
        $(
            impl ParseCoeff for $t {
                type Error = winnow::error::ErrMode<ContextError>;

                fn parse_coeff(input: &mut &str, _inside_bracket: bool) -> Result<Self, Self::Error> {
                    dec_int.parse_next(input)
                }
            }
        )*
    };
}

impl_parse_coeff_int!(i8, i16, i32, i64, i128);

macro_rules! impl_parse_coeff_uint {
    ($($t:ty), *) => {
        $(
            impl ParseCoeff for $t {
                type Error = winnow::error::ErrMode<ContextError>;

                fn parse_coeff(input: &mut &str, _inside_bracket: bool) -> Result<Self, Self::Error> {
                    dec_uint.parse_next(input)
                }
            }
        )*
    };
}

impl_parse_coeff_uint!(u8, u16, u32, u64, u128);

macro_rules! impl_parse_coeff_float {
    ($($t:ty), *) => {
        $(
            impl ParseCoeff for $t {
                type Error = winnow::error::ErrMode<ContextError>;

                fn parse_coeff(input: &mut &str, _inside_bracket: bool) -> Result<Self, Self::Error> {
                    float.parse_next(input)
                }
            }
        )*
    };
}

impl_parse_coeff_float!(f32, f64);

impl<Var, C> ParseCoeff for Polynomial<Var, C>
where
    Var: FromStr + PartialEq,
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
{
    type Error = winnow::error::ErrMode<ContextError>;

    fn parse_coeff(input: &mut &str, inside_bracket: bool) -> Result<Self, Self::Error> {
        if inside_bracket {
            poly.parse_next(input)
        } else {
            monomial.parse_next(input)
                .map(
                    |Monomial { c, var, pow }| if let Some(var) = var {
                        Polynomial::new(var, pow, vec![c])
                    } else {
                        debug_assert_eq!(pow, 0);
                        Polynomial::Const(c)
                    })
        }
    }
}

impl<Var, C> ParseCoeff for Laurent<Var, C>
where
    Var: Clone + Debug + FromStr + PartialEq,
    C: AddAssign + Coeff + Neg<Output = C> + ParseCoeff,
{
    type Error = winnow::error::ErrMode<ContextError>;

    fn parse_coeff(input: &mut &str, inside_bracket: bool) -> Result<Self, Self::Error> {
        if inside_bracket {
            laurent.parse_next(input)
        } else {
            alt((
                monomial
                    .map(
                        |Monomial { c, var, pow }| if let Some(var) = var {
                            Polynomial::new(var, pow, vec![c]).into()
                        } else {
                            debug_assert_eq!(pow, 0);
                            Polynomial::Const(c).into()
                        }),
                cutoff
                    .map(|(var, pow)| Series::new(var, pow, vec![]).into())
            )).parse_next(input)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{var, O};

    #[test]
    fn constants() {
        let p: Polynomial<String, i32> = "0".parse().unwrap();
        assert!(p.is_zero());

        let p: Polynomial<String, i32> = "1".parse().unwrap();
        assert!(p.is_one());

        let p: Polynomial<String, i32> = "2".parse().unwrap();
        assert_eq!(p, Polynomial::from_const(2));

        let p: Polynomial<String, i32> = "2 + 2".parse().unwrap();
        assert_eq!(p, Polynomial::from_const(4));
    }

    #[test]
    fn linear() {
        let p: Polynomial<String, i32> = "x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![1]));

        let p: Polynomial<String, i32> = "-x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![-1]));

        let p: Polynomial<String, i32> = "x - x".parse().unwrap();
        assert!(p.is_zero());

        let p: Polynomial<String, i32> = "2*x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![2]));

        let p: Polynomial<String, i32> = "+2*x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![2]));

        let p: Polynomial<String, i32> = "-2*x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![-2]));

        let p: Polynomial<String, i32> = "x+x".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![2]));
    }

    #[test]
    fn monomial() {
        let p: Polynomial<String, i32> = "x^0".parse().unwrap();
        assert!(p.is_one());

        let p: Polynomial<String, i32> = "x^-1".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), -1, vec![1]));

        let p: Polynomial<String, i32> = "x^1".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 1, vec![1]));

        let p: Polynomial<String, i32> = "2*x^2".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 2, vec![2]));

        let p: Polynomial<String, i32> = "1/x^2".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), -2, vec![1]));

        let p: Polynomial<String, i32> = "1/x^0".parse().unwrap();
        assert!(p.is_one());
    }

    #[test]
    fn poly() {
        let p: Polynomial<String, i32> = "1 + 0*x^0".parse().unwrap();
        assert!(p.is_one());

        let p: Polynomial<String, i32> = "2 + 3*x^0".parse().unwrap();
        assert_eq!(p, Polynomial::new("x".to_owned(), 0, vec![5]));

        let p: Polynomial<String, i32> = "2/x - 3*x^3".parse().unwrap();
        assert_eq!(
            p,
            Polynomial::new("x".to_owned(), -1, vec![2, 0, 0, 0, -3])
        );

        let p: Polynomial<String, i32> = "(0)*x^10".parse().unwrap();
        assert!(p.is_zero());
    }

    #[test]
    fn nested_poly() {
        var!(X);
        var!(Y);

        let p: Polynomial<X, Polynomial<Y, i32>> = "x + y".parse().unwrap();
        let res = Polynomial::new(
            X,
            0,
            vec![Polynomial::new(Y, 1, vec![1]), Polynomial::one()]
        );
        assert_eq!(p, res);

        let p: Polynomial<X, Polynomial<Y, i32>> = "(1 + y)*x".parse().unwrap();
        let res = Polynomial::new(
            X,
            1,
            vec![Polynomial::new(Y, 0, vec![1, 1])]
        );
        assert_eq!(p, res);
    }

    #[test]
    fn series() {
        var!(X);
        let p: Series<X, i32> = "1 + 0*x^0 + O(x)".parse().unwrap();
        let res = Series::new(X, 0, vec![1]);
        assert_eq!(p, res);

        let p: Series<X, i32> = "2 + 3*x^0 + O(x^2)".parse().unwrap();
        assert_eq!(p, Series::new(X, 0, vec![5, 0]));

        let p: Series<X, i32> = "2/x - 3*x^3 + O(x^4)".parse().unwrap();
        assert_eq!(
            p,
            Series::new(X, -1, vec![2, 0, 0, 0, -3])
        );

        let p: Series<X, i32> = "O(x^10)".parse().unwrap();
        assert_eq!(p, O!(X^10));
    }

    #[test]
    fn laurent() {
        var!(X);
        var!(Y);

        let l: Laurent<X, i32> = "1 + 0*x^0 + O(x)".parse().unwrap();
        let res = Laurent::from(Series::new(X, 0, vec![1]));
        assert_eq!(l, res);

        let l: Laurent<Y, Laurent<X, i32>> = "y + (1/x + 3*x^2 + O(x^3))*y^2".parse().unwrap();
        let res = Laurent::from(
            Polynomial::new(
                Y,
                1,
                vec![
                    Laurent::one(),
                    Laurent::from(Series::new(X, -1, vec![1, 0, 0, 3]))
                ]
            )
        );
        assert_eq!(l, res);
    }

    #[test]
    fn errs() {
        let p: Result<Polynomial<String, i32>, _> = "x + y".parse();
        assert!(p.is_err());

        let p: Result<Polynomial<String, i32>, _> = "x + ".parse();
        assert!(p.is_err());

        let p: Result<Polynomial<String, i32>, _> = "x@".parse();
        assert!(p.is_err());

        let p: Result<Polynomial<String, i32>, _> = "".parse();
        assert!(p.is_err());
    }
}
