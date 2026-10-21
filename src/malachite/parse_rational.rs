use winnow::{
    Parser,
    ascii::digit1,
    combinator::{delimited, opt, preceded},
    error::ContextError,
};

use crate::parse::ParseCoeff;

use super::Rational;

impl ParseCoeff for Rational {
    type Error = winnow::error::ErrMode<ContextError>;

    fn parse_coeff(
        input: &mut &str,
        _inside_bracket: bool,
    ) -> Result<Self, Self::Error> {
        delimited(opt('-'), digit1, preceded('/', digit1))
            .take()
            .parse_to()
            .parse_next(input)
    }
}
