use winnow::{
    Parser,
    ascii::digit1,
    combinator::{opt, preceded},
    error::ContextError,
};

use crate::parse::ParseCoeff;

use super::Integer;

impl ParseCoeff for Integer {
    type Error = winnow::error::ErrMode<ContextError>;

    fn parse_coeff(
        input: &mut &str,
        _inside_bracket: bool,
    ) -> Result<Self, Self::Error> {
        preceded(opt('-'), digit1)
            .take()
            .parse_to()
            .parse_next(input)
    }
}
