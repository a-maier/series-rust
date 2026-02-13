use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

use num_traits::{One, Zero};

use crate::{zero_ref::zero_ref, Coeff, MulInverse, anon_series::AnonSeries};

/// A sum type of a Laurent series and its coefficient ("inner") type.
///
/// The main reason for this struct is to allow nested Laurent series.
///
/// [Series] cannot be nested, because Laurent series have neither an
/// additive nor a multiplicative identity. That means that [Series]
/// cannot implement [Coeff] and is not a valid coefficient type.
///
/// However, [Coeff] itself has identity elements. They are also
/// identities for Laurent series, so the sum type of a Laurent series
/// and its coefficient type can be nested.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq, Debug, Clone, Hash, Ord, PartialOrd)]
pub enum InnerSeries<C: Coeff> {
    Series(AnonSeries<C>),
    Inner(C),
}

impl<C: Coeff + Send + Sync + 'static> InnerSeries<C> {
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        use InnerSeries::*;
        match self {
            Series(s) => s.coeff(pow),
            Inner(i) => Some(if pow == 0 { i } else { zero_ref() }),
        }
    }
}

impl<C: Coeff + Default> Default for InnerSeries<C> {
    fn default() -> Self {
        Self::Inner(C::default())
    }
}

impl<C: Coeff> From<AnonSeries<C>> for InnerSeries<C> {
    fn from(s: AnonSeries<C>) -> Self {
        Self::Series(s)
    }
}

impl<C: Coeff> From<C> for InnerSeries<C> {
    fn from(c: C) -> Self {
        Self::Inner(c)
    }
}

impl<C: Coeff> Neg for InnerSeries<C>
where
    C: Neg<Output = C>,
    AnonSeries<C>: Neg<Output = AnonSeries<C>>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.neg().into(),
            Inner(c) => c.neg().into(),
        }
    }
}

impl<'a, C: Coeff> Neg for &'a InnerSeries<C>
where
    &'a C: Neg<Output = C>,
    &'a AnonSeries<C>: Neg<Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn neg(self) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.neg().into(),
            Inner(c) => c.neg().into(),
        }
    }
}

impl<C: Coeff> MulInverse for InnerSeries<C>
where
    C: MulInverse<Output = C>,
    AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    type Output = Self;

    fn mul_inverse(self) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_inverse().into(),
            Inner(c) => c.mul_inverse().into(),
        }
    }
}

impl<'a, C: Coeff> MulInverse for &'a InnerSeries<C>
where
    &'a C: MulInverse<Output = C>,
    &'a AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul_inverse(self) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_inverse().into(),
            Inner(c) => c.mul_inverse().into(),
        }
    }
}

impl<C: Coeff> Add<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Add<Output = AnonSeries<C>> + Add<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<C: Coeff> Add<InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Add<Output = AnonSeries<C>> + Add<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<C: Coeff> Add<C> for InnerSeries<C>
where
    AnonSeries<C>: Add<C, Output = AnonSeries<C>>,
    C: Add<Output = C>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<C: Coeff> Add for InnerSeries<C>
where
    InnerSeries<C>: Add<C, Output = Self> + Add<AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add(s),
            Inner(c) => self.add(c),
        }
    }
}

impl<'a, C: Coeff> Add<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Add<&'a AnonSeries<C>, Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Add<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Add<InnerSeries<C>> for &'a AnonSeries<C>
where
    AnonSeries<C>: Add<&'a AnonSeries<C>, Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Add<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, C: Coeff> Add<&'a C> for InnerSeries<C>
where
    C: Add<&'a C, Output = C>,
    AnonSeries<C>: Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => c.add(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Add<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>:
        Add<&'a C, Output = Self> + Add<&'a AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add(s),
            Inner(c) => self.add(c),
        }
    }
}

impl<'a, C: Coeff> Add<AnonSeries<C>> for &'a InnerSeries<C>
where
    AnonSeries<C>:
        Add<&'a AnonSeries<C>, Output = AnonSeries<C>> + Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.add(s).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Add<&'a InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>:
        Add<&'a AnonSeries<C>, Output = AnonSeries<C>> + Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, C: Coeff> Add<C> for &'a InnerSeries<C>
where
    C: Add<&'a C, Output = C>,
    &'a AnonSeries<C>: Add<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Add<InnerSeries<C>> for &'a InnerSeries<C>
where
    InnerSeries<C>: Add<&'a C, Output = InnerSeries<C>>
        + Add<&'a AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, 'b, C: Coeff> Add<&'b AnonSeries<C>> for &'a InnerSeries<C>
where
    &'a AnonSeries<C>: Add<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'b AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Add<&'a InnerSeries<C>> for &'b AnonSeries<C>
where
    &'a AnonSeries<C>: Add<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, 'b, C: Coeff> Add<&'b C> for &'a InnerSeries<C>
where
    &'a C: Add<&'b C, Output = C>,
    &'a AnonSeries<C>: Add<&'b C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'b C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => c.add(rhs).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Add<&'b InnerSeries<C>> for &'a InnerSeries<C>
where
    &'a InnerSeries<C>: Add<&'b C, Output = InnerSeries<C>>
        + Add<&'b AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'b InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add(s),
            Inner(c) => self.add(c),
        }
    }
}

impl<C: Coeff + Default> AddAssign<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: AddAssign + AddAssign<C>,
{
    fn add_assign(&mut self, mut rhs: AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => {
                rhs.add_assign(std::mem::take(c));
                *self = rhs.into();
            }
        }
    }
}

impl<C: Coeff + Default> AddAssign<C> for InnerSeries<C>
where
    AnonSeries<C>: AddAssign<C>,
    C: AddAssign,
{
    fn add_assign(&mut self, rhs: C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => c.add_assign(rhs),
        }
    }
}

impl<C: Coeff + Default> AddAssign for InnerSeries<C>
where
    InnerSeries<C>: AddAssign<AnonSeries<C>> + AddAssign<C>,
{
    fn add_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add_assign(s),
            Inner(c) => self.add_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> AddAssign<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: AddAssign<&'a AnonSeries<C>>,
    &'a AnonSeries<C>: Add<C, Output = AnonSeries<C>>,
{
    fn add_assign(&mut self, rhs: &'a AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => {
                let res = rhs.add(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<'a, C: Coeff + Default> AddAssign<&'a C> for InnerSeries<C>
where
    AnonSeries<C>: AddAssign<&'a C>,
    C: AddAssign<&'a C>,
{
    fn add_assign(&mut self, rhs: &'a C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => c.add_assign(rhs),
        }
    }
}

impl<'a, C: Coeff + Default> AddAssign<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>: AddAssign<&'a AnonSeries<C>> + AddAssign<&'a C>,
{
    fn add_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add_assign(s),
            Inner(c) => self.add_assign(c),
        }
    }
}

impl<C: Coeff> Sub<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Sub<Output = AnonSeries<C>>
        + Add<C, Output = AnonSeries<C>>
        + Neg<Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<C: Coeff> Sub<InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Sub<Output = AnonSeries<C>> + Sub<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub(s).into(),
            Inner(c) => self.sub(c).into(),
        }
    }
}

impl<C: Coeff> Sub<C> for InnerSeries<C>
where
    AnonSeries<C>: Sub<C, Output = AnonSeries<C>>,
    C: Sub<Output = C>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => c.sub(rhs).into(),
        }
    }
}

impl<C: Coeff> Sub for InnerSeries<C>
where
    InnerSeries<C>: Sub<C, Output = Self> + Sub<AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub(s),
            Inner(c) => self.sub(c),
        }
    }
}

impl<'a, C: Coeff> Sub<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Sub<&'a AnonSeries<C>, Output = AnonSeries<C>> + Neg<Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Sub<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.sub(c).neg().into(),
        }
    }
}

impl<'a, C: Coeff> Sub<InnerSeries<C>> for &'a AnonSeries<C>
where
    InnerSeries<C>: Neg<Output = InnerSeries<C>>
        + Add<&'a AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.neg().add(self)
    }
}

impl<'a, C: Coeff> Sub<&'a C> for InnerSeries<C>
where
    C: Sub<&'a C, Output = C>,
    AnonSeries<C>: Sub<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => c.sub(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Sub<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>:
        Sub<&'a C, Output = Self> + Sub<&'a AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub(s),
            Inner(c) => self.sub(c),
        }
    }
}

impl<'a, C: Coeff> Sub<AnonSeries<C>> for &'a InnerSeries<C>
where
    AnonSeries<C>: Add<&'a AnonSeries<C>, Output = AnonSeries<C>>
        + Add<&'a C, Output = AnonSeries<C>>
        + Neg<Output = AnonSeries<C>>,
    C: Neg<Output = C>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.neg().add(s).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Sub<&'a InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Add<InnerSeries<C>, Output = InnerSeries<C>>,
    &'a InnerSeries<C>: Neg<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        self.add(rhs.neg())
    }
}

impl<'a, C: Coeff> Sub<C> for &'a InnerSeries<C>
where
    &'a C: Sub<C, Output = C>,
    &'a AnonSeries<C>: Sub<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => c.sub(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Sub<InnerSeries<C>> for &'a InnerSeries<C>
where
    InnerSeries<C>: Add<&'a C, Output = InnerSeries<C>>
        + Add<&'a AnonSeries<C>, Output = InnerSeries<C>>
        + Neg<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.neg().add(self)
    }
}

impl<'a, 'b, C: Coeff> Sub<&'b AnonSeries<C>> for &'a InnerSeries<C>
where
    &'a AnonSeries<C>: Sub<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: Neg<Output = AnonSeries<C>>,
    AnonSeries<C>: Add<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'b AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Sub<&'a InnerSeries<C>> for &'b AnonSeries<C>
where
    &'b AnonSeries<C>:
        Sub<&'a C, Output = AnonSeries<C>> + Sub<&'a AnonSeries<C>, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub(s).into(),
            Inner(c) => self.sub(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Sub<&'b C> for &'a InnerSeries<C>
where
    &'a C: Sub<&'b C, Output = C>,
    &'a AnonSeries<C>: Sub<&'b C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'b C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => c.sub(rhs).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Sub<&'b InnerSeries<C>> for &'a InnerSeries<C>
where
    &'a InnerSeries<C>: Sub<&'b C, Output = InnerSeries<C>>
        + Sub<&'b AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'b InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub(s),
            Inner(c) => self.sub(c),
        }
    }
}

impl<C: Coeff + Default> SubAssign<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: SubAssign + AddAssign<C> + Neg<Output = AnonSeries<C>>,
{
    fn sub_assign(&mut self, rhs: AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => {
                let mut res = rhs.neg();
                res.add_assign(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<C: Coeff + Default> SubAssign<C> for InnerSeries<C>
where
    AnonSeries<C>: SubAssign<C>,
    C: SubAssign,
{
    fn sub_assign(&mut self, rhs: C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => c.sub_assign(rhs),
        }
    }
}

impl<C: Coeff + Default> SubAssign for InnerSeries<C>
where
    InnerSeries<C>: SubAssign<AnonSeries<C>> + SubAssign<C>,
{
    fn sub_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub_assign(s),
            Inner(c) => self.sub_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> SubAssign<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: AddAssign<C> + SubAssign<&'a AnonSeries<C>>,
    &'a AnonSeries<C>: Neg<Output = AnonSeries<C>>,
{
    fn sub_assign(&mut self, rhs: &'a AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => {
                let mut res = rhs.neg();
                res.add_assign(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<'a, C: Coeff + Default> SubAssign<&'a C> for InnerSeries<C>
where
    AnonSeries<C>: SubAssign<&'a C>,
    C: SubAssign<&'a C>,
{
    fn sub_assign(&mut self, rhs: &'a C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => c.sub_assign(rhs),
        }
    }
}

impl<'a, C: Coeff + Default> SubAssign<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>: SubAssign<&'a AnonSeries<C>> + SubAssign<&'a C>,
{
    fn sub_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub_assign(s),
            Inner(c) => self.sub_assign(c),
        }
    }
}

impl<C: Coeff> Mul<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Mul<Output = AnonSeries<C>> + Mul<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<C: Coeff> Mul<InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Mul<Output = AnonSeries<C>> + Mul<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<C: Coeff> Mul<C> for InnerSeries<C>
where
    AnonSeries<C>: Mul<C, Output = AnonSeries<C>>,
    C: Mul<Output = C>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<C: Coeff> Mul for InnerSeries<C>
where
    InnerSeries<C>: Mul<C, Output = Self> + Mul<AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul(s),
            Inner(c) => self.mul(c),
        }
    }
}

impl<'a, C: Coeff> Mul<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Mul<&'a AnonSeries<C>, Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Mul<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<InnerSeries<C>> for &'a AnonSeries<C>
where
    AnonSeries<C>: Mul<&'a AnonSeries<C>, Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Mul<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, C: Coeff> Mul<&'a C> for InnerSeries<C>
where
    C: Mul<&'a C, Output = C>,
    AnonSeries<C>: Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => c.mul(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>:
        Mul<&'a C, Output = Self> + Mul<&'a AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul(s),
            Inner(c) => self.mul(c),
        }
    }
}

impl<'a, C: Coeff> Mul<AnonSeries<C>> for &'a InnerSeries<C>
where
    AnonSeries<C>:
        Mul<&'a AnonSeries<C>, Output = AnonSeries<C>> + Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.mul(s).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<&'a InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>:
        Mul<&'a AnonSeries<C>, Output = AnonSeries<C>> + Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, C: Coeff> Mul<C> for &'a InnerSeries<C>
where
    C: Mul<&'a C, Output = C>,
    &'a AnonSeries<C>: Mul<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<InnerSeries<C>> for &'a InnerSeries<C>
where
    InnerSeries<C>: Mul<&'a C, Output = InnerSeries<C>>
        + Mul<&'a AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b AnonSeries<C>> for &'a InnerSeries<C>
where
    &'a AnonSeries<C>: Mul<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'b AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Mul<&'a InnerSeries<C>> for &'b AnonSeries<C>
where
    &'a AnonSeries<C>: Mul<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b C> for &'a InnerSeries<C>
where
    &'a C: Mul<&'b C, Output = C>,
    &'a AnonSeries<C>: Mul<&'b C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'b C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => c.mul(rhs).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b InnerSeries<C>> for &'a InnerSeries<C>
where
    &'a InnerSeries<C>: Mul<&'b C, Output = InnerSeries<C>>
        + Mul<&'b AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'b InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul(s),
            Inner(c) => self.mul(c),
        }
    }
}

impl<C: Coeff + Default> MulAssign<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: MulAssign + MulAssign<C>,
{
    fn mul_assign(&mut self, mut rhs: AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => {
                rhs.mul_assign(std::mem::take(c));
                *self = rhs.into();
            }
        }
    }
}

impl<C: Coeff + Default> MulAssign<C> for InnerSeries<C>
where
    AnonSeries<C>: MulAssign<C>,
    C: MulAssign,
{
    fn mul_assign(&mut self, rhs: C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => c.mul_assign(rhs),
        }
    }
}

impl<C: Coeff + Default> MulAssign for InnerSeries<C>
where
    InnerSeries<C>: MulAssign<AnonSeries<C>> + MulAssign<C>,
{
    fn mul_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul_assign(s),
            Inner(c) => self.mul_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> MulAssign<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: MulAssign<&'a AnonSeries<C>>,
    &'a AnonSeries<C>: Mul<C, Output = AnonSeries<C>>,
{
    fn mul_assign(&mut self, rhs: &'a AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => {
                let res = rhs.mul(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<'a, C: Coeff + Default> MulAssign<&'a C> for InnerSeries<C>
where
    AnonSeries<C>: MulAssign<&'a C>,
    C: MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => c.mul_assign(rhs),
        }
    }
}

impl<'a, C: Coeff + Default> MulAssign<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>: MulAssign<&'a AnonSeries<C>> + MulAssign<&'a C>,
{
    fn mul_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul_assign(s),
            Inner(c) => self.mul_assign(c),
        }
    }
}

impl<C: Coeff> Div<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: Div<Output = AnonSeries<C>>
        + Mul<C, Output = AnonSeries<C>>
        + MulInverse<Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<C: Coeff> Div<InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Div<Output = AnonSeries<C>> + Div<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div(s).into(),
            Inner(c) => self.div(c).into(),
        }
    }
}

impl<C: Coeff> Div<C> for InnerSeries<C>
where
    AnonSeries<C>: Div<C, Output = AnonSeries<C>>,
    C: Div<Output = C>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => c.div(rhs).into(),
        }
    }
}

impl<C: Coeff> Div for InnerSeries<C>
where
    InnerSeries<C>: Div<C, Output = Self> + Div<AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div(s),
            Inner(c) => self.div(c),
        }
    }
}

impl<'a, C: Coeff> Div<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>:
        Div<&'a AnonSeries<C>, Output = AnonSeries<C>> + MulInverse<Output = AnonSeries<C>>,
    &'a AnonSeries<C>: Div<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.div(c).mul_inverse().into(),
        }
    }
}

impl<'a, C: Coeff> Div<InnerSeries<C>> for &'a AnonSeries<C>
where
    InnerSeries<C>: MulInverse<Output = InnerSeries<C>>
        + Mul<&'a AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul_inverse().mul(self)
    }
}

impl<'a, C: Coeff> Div<&'a C> for InnerSeries<C>
where
    C: Div<&'a C, Output = C>,
    AnonSeries<C>: Div<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => c.div(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Div<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>:
        Div<&'a C, Output = Self> + Div<&'a AnonSeries<C>, Output = Self>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div(s),
            Inner(c) => self.div(c),
        }
    }
}

impl<'a, C: Coeff> Div<AnonSeries<C>> for &'a InnerSeries<C>
where
    AnonSeries<C>: Mul<&'a AnonSeries<C>, Output = AnonSeries<C>>
        + Mul<&'a C, Output = AnonSeries<C>>
        + MulInverse<Output = AnonSeries<C>>,
    C: MulInverse<Output = C>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.mul_inverse().mul(s).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Div<&'a InnerSeries<C>> for AnonSeries<C>
where
    AnonSeries<C>: Mul<InnerSeries<C>, Output = InnerSeries<C>>,
    &'a InnerSeries<C>: MulInverse<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        self.mul(rhs.mul_inverse())
    }
}

impl<'a, C: Coeff> Div<C> for &'a InnerSeries<C>
where
    &'a C: Div<C, Output = C>,
    &'a AnonSeries<C>: Div<C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => c.div(rhs).into(),
        }
    }
}

impl<'a, C: Coeff> Div<InnerSeries<C>> for &'a InnerSeries<C>
where
    InnerSeries<C>: Mul<&'a C, Output = InnerSeries<C>>
        + Mul<&'a AnonSeries<C>, Output = InnerSeries<C>>
        + MulInverse<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul_inverse().mul(self)
    }
}

impl<'a, 'b, C: Coeff> Div<&'b AnonSeries<C>> for &'a InnerSeries<C>
where
    &'a AnonSeries<C>: Div<&'b AnonSeries<C>, Output = AnonSeries<C>>,
    &'b AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
    AnonSeries<C>: Mul<&'a C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'b AnonSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Div<&'a InnerSeries<C>> for &'b AnonSeries<C>
where
    &'b AnonSeries<C>:
        Div<&'a C, Output = AnonSeries<C>> + Div<&'a AnonSeries<C>, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div(s).into(),
            Inner(c) => self.div(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Div<&'b C> for &'a InnerSeries<C>
where
    &'a C: Div<&'b C, Output = C>,
    &'a AnonSeries<C>: Div<&'b C, Output = AnonSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'b C) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => c.div(rhs).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Div<&'b InnerSeries<C>> for &'a InnerSeries<C>
where
    &'a InnerSeries<C>: Div<&'b C, Output = InnerSeries<C>>
        + Div<&'b AnonSeries<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'b InnerSeries<C>) -> Self::Output {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div(s),
            Inner(c) => self.div(c),
        }
    }
}

impl<C: Coeff + Default> DivAssign<AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: DivAssign + MulAssign<C> + MulInverse<Output = AnonSeries<C>>,
{
    fn div_assign(&mut self, rhs: AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => {
                let mut res = rhs.mul_inverse();
                res.mul_assign(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<C: Coeff + Default> DivAssign<C> for InnerSeries<C>
where
    AnonSeries<C>: DivAssign<C>,
    C: DivAssign,
{
    fn div_assign(&mut self, rhs: C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => c.div_assign(rhs),
        }
    }
}

impl<C: Coeff + Default> DivAssign for InnerSeries<C>
where
    InnerSeries<C>: DivAssign<AnonSeries<C>> + DivAssign<C>,
{
    fn div_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div_assign(s),
            Inner(c) => self.div_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> DivAssign<&'a AnonSeries<C>> for InnerSeries<C>
where
    AnonSeries<C>: MulAssign<C> + DivAssign<&'a AnonSeries<C>>,
    &'a AnonSeries<C>: MulInverse<Output = AnonSeries<C>>,
{
    fn div_assign(&mut self, rhs: &'a AnonSeries<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => {
                let mut res = rhs.mul_inverse();
                res.mul_assign(std::mem::take(c));
                *self = res.into();
            }
        }
    }
}

impl<'a, C: Coeff + Default> DivAssign<&'a C> for InnerSeries<C>
where
    AnonSeries<C>: DivAssign<&'a C>,
    C: DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a C) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => c.div_assign(rhs),
        }
    }
}

impl<'a, C: Coeff + Default> DivAssign<&'a InnerSeries<C>> for InnerSeries<C>
where
    InnerSeries<C>: DivAssign<&'a AnonSeries<C>> + DivAssign<&'a C>,
{
    fn div_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div_assign(s),
            Inner(c) => self.div_assign(c),
        }
    }
}

// TODO: pow, exp, ln

impl<C: Coeff> Zero for InnerSeries<C>
where
    InnerSeries<C>: Add<Output = Self>,
{
    fn zero() -> Self {
        C::zero().into()
    }

    fn is_zero(&self) -> bool {
        use InnerSeries::*;
        match self {
            Series(_) => false,
            Inner(c) => c.is_zero(),
        }
    }
}

impl<C: Coeff> One for InnerSeries<C>
where
    InnerSeries<C>: Mul<Output = Self>,
{
    fn one() -> Self {
        C::one().into()
    }

    fn is_one(&self) -> bool {
        use InnerSeries::*;
        match self {
            Series(_) => false,
            Inner(c) => c.is_one(),
        }
    }
}
