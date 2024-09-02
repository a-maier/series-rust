use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{One, Zero};

use crate::{zero_ref::zero_ref, Coeff, MulInverse, Series};

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
    Series(Series<C>),
    Inner(C),
}

impl<C: Coeff + Send + Sync + 'static> InnerSeries<C> {
    pub fn coeff(&self, pow: isize) -> Option<&C> {
        use InnerSeries::*;
        match self {
            Series(s) => s.coeff(pow),
            Inner(i) => Some(
                if pow == 0 {
                    i
                } else {
                    zero_ref()
                }
            ),
        }
    }
}

impl<C: Coeff + Default> Default for InnerSeries<C> {
    fn default() -> Self {
        Self::Inner(C::default())
    }
}

impl<C: Coeff> From<Series<C>> for InnerSeries<C> {
    fn from(s: Series<C>) -> Self {
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
    Series<C>: Neg<Output = Series<C>>,
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
    &'a Series<C>: Neg<Output = Series<C>>,
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
    Series<C>: MulInverse<Output = Series<C>>,
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
    &'a Series<C>: MulInverse<Output = Series<C>>,
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

impl<C: Coeff> Add<Series<C>> for InnerSeries<C>
where
    Series<C>: Add<Output = Series<C>> + Add<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<C: Coeff> Add<InnerSeries<C>> for Series<C>
where
    Series<C>: Add<Output = Series<C>> + Add<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<C: Coeff> Add<C> for InnerSeries<C>
where
    Series<C>: Add<C, Output = Series<C>>,
    C: Add<Output = C>
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
    InnerSeries<C>: Add<C, Output = Self> + Add<Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Add<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: Add<&'a Series<C>, Output = Series<C>>,
    &'a Series<C>: Add<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Add<InnerSeries<C>> for &'a Series<C>
where
    Series<C>: Add<&'a Series<C>, Output = Series<C>>,
    &'a Series<C>: Add<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, C: Coeff> Add<&'a C> for InnerSeries<C>
where
    C: Add<&'a C, Output = C>,
    Series<C>: Add<&'a C, Output = Series<C>>,
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
    InnerSeries<C>: Add<&'a C, Output = Self> + Add<&'a Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Add<Series<C>> for &'a InnerSeries<C>
where
    Series<C>: Add<&'a Series<C>, Output = Series<C>> + Add<&'a C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.add(s).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Add<&'a InnerSeries<C>> for Series<C>
where
    Series<C>: Add<&'a Series<C>, Output = Series<C>> + Add<&'a C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, C: Coeff> Add<C> for &'a InnerSeries<C>
where
    C: Add<&'a C, Output = C>,
    &'a Series<C>: Add<C, Output = Series<C>>,
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
    InnerSeries<C>: Add<&'a C, Output = InnerSeries<C>> + Add<&'a Series<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, 'b, C: Coeff> Add<&'b Series<C>> for &'a InnerSeries<C>
where
    &'a Series<C>: Add<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>: Add<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'b Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.add(rhs).into(),
            Inner(c) => rhs.add(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Add<&'a InnerSeries<C>> for &'b Series<C>
where
    &'a Series<C>: Add<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>: Add<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn add(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.add(self)
    }
}

impl<'a, 'b, C: Coeff> Add<&'b C> for &'a InnerSeries<C>
where
    &'a C: Add<&'b C, Output = C>,
    &'a Series<C>: Add<&'b C, Output = Series<C>>,
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
    &'a InnerSeries<C>: Add<&'b C, Output = InnerSeries<C>> + Add<&'b Series<C>, Output = InnerSeries<C>>,
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

impl<C: Coeff + Default> AddAssign<Series<C>> for InnerSeries<C>
where
    Series<C>: AddAssign + AddAssign<C>
{
    fn add_assign(&mut self, mut rhs: Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => {
                rhs.add_assign(std::mem::take(c));
                *self = rhs.into();
            },
        }
    }
}

impl<C: Coeff + Default> AddAssign<C> for InnerSeries<C>
where
    Series<C>: AddAssign<C>,
    C: AddAssign
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
    InnerSeries<C>: AddAssign<Series<C>> + AddAssign<C>
{
    fn add_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add_assign(s),
            Inner(c) => self.add_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> AddAssign<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: AddAssign<&'a Series<C>>,
    &'a Series<C>: Add<C, Output = Series<C>>,
{
    fn add_assign(&mut self, rhs: &'a Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.add_assign(rhs),
            Inner(c) => {
                let res = rhs.add(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<'a, C: Coeff + Default> AddAssign<&'a C> for InnerSeries<C>
where
    Series<C>: AddAssign<&'a C>,
    C: AddAssign<&'a C>
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
    InnerSeries<C>: AddAssign<&'a Series<C>> + AddAssign<&'a C>
{
    fn add_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.add_assign(s),
            Inner(c) => self.add_assign(c),
        }
    }
}

impl<C: Coeff> Sub<Series<C>> for InnerSeries<C>
where
    Series<C>: Sub<Output = Series<C>> + Add<C, Output = Series<C>> + Neg<Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<C: Coeff> Sub<InnerSeries<C>> for Series<C>
where
    Series<C>: Sub<Output = Series<C>> + Sub<C, Output = Series<C>>,
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
    Series<C>: Sub<C, Output = Series<C>>,
    C: Sub<Output = C>
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
    InnerSeries<C>: Sub<C, Output = Self> + Sub<Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Sub<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: Sub<&'a Series<C>, Output = Series<C>> + Neg<Output = Self>,
    &'a Series<C>: Sub<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'a Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.sub(c).neg().into(),
        }
    }
}

impl<'a, C: Coeff> Sub<InnerSeries<C>> for &'a Series<C>
where
    InnerSeries<C>: Neg<Output = InnerSeries<C>> + Add<&'a Series<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.neg().add(self)
    }
}

impl<'a, C: Coeff> Sub<&'a C> for InnerSeries<C>
where
    C: Sub<&'a C, Output = C>,
    Series<C>: Sub<&'a C, Output = Series<C>>,
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
    InnerSeries<C>: Sub<&'a C, Output = Self> + Sub<&'a Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Sub<Series<C>> for &'a InnerSeries<C>
where
    Series<C>: Add<&'a Series<C>, Output = Series<C>> + Add<&'a C, Output = Series<C>> + Neg<Output = Series<C>>,
C: Neg<Output = C>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.neg().add(s).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<'a, C: Coeff> Sub<&'a InnerSeries<C>> for Series<C>
where
    Series<C>: Add<InnerSeries<C>, Output = InnerSeries<C>>,
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
    &'a Series<C>: Sub<C, Output = Series<C>>,
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
    InnerSeries<C>: Add<&'a C, Output = InnerSeries<C>> + Add<&'a Series<C>, Output = InnerSeries<C>> + Neg<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.neg().add(self)
    }
}

impl<'a, 'b, C: Coeff> Sub<&'b Series<C>> for &'a InnerSeries<C>
where
    &'a Series<C>: Sub<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>:  Neg<Output = Series<C>>,
    Series<C>: Add<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn sub(self, rhs: &'b Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub(rhs).into(),
            Inner(c) => rhs.neg().add(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Sub<&'a InnerSeries<C>> for &'b Series<C>
where
    &'b Series<C>: Sub<&'a C, Output = Series<C>> + Sub<&'a Series<C>, Output = Series<C>>,
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
    &'a Series<C>: Sub<&'b C, Output = Series<C>>,
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
    &'a InnerSeries<C>: Sub<&'b C, Output = InnerSeries<C>> + Sub<&'b Series<C>, Output = InnerSeries<C>>,
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

impl<C: Coeff + Default> SubAssign<Series<C>> for InnerSeries<C>
where
    Series<C>: SubAssign + AddAssign<C> + Neg<Output = Series<C>>
{
    fn sub_assign(&mut self, rhs: Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => {
                let mut res = rhs.neg();
                res.add_assign(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<C: Coeff + Default> SubAssign<C> for InnerSeries<C>
where
    Series<C>: SubAssign<C>,
    C: SubAssign
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
    InnerSeries<C>: SubAssign<Series<C>> + SubAssign<C>
{
    fn sub_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub_assign(s),
            Inner(c) => self.sub_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> SubAssign<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: AddAssign<C> + SubAssign<&'a Series<C>>,
    &'a Series<C>: Neg<Output = Series<C>>,
{
    fn sub_assign(&mut self, rhs: &'a Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.sub_assign(rhs),
            Inner(c) => {
                let mut res = rhs.neg();
                res.add_assign(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<'a, C: Coeff + Default> SubAssign<&'a C> for InnerSeries<C>
where
    Series<C>: SubAssign<&'a C>,
    C: SubAssign<&'a C>
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
    InnerSeries<C>: SubAssign<&'a Series<C>> + SubAssign<&'a C>
{
    fn sub_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.sub_assign(s),
            Inner(c) => self.sub_assign(c),
        }
    }
}

impl<C: Coeff> Mul<Series<C>> for InnerSeries<C>
where
    Series<C>: Mul<Output = Series<C>> + Mul<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<C: Coeff> Mul<InnerSeries<C>> for Series<C>
where
    Series<C>: Mul<Output = Series<C>> + Mul<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<C: Coeff> Mul<C> for InnerSeries<C>
where
    Series<C>: Mul<C, Output = Series<C>>,
    C: Mul<Output = C>
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
    InnerSeries<C>: Mul<C, Output = Self> + Mul<Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Mul<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: Mul<&'a Series<C>, Output = Series<C>>,
    &'a Series<C>: Mul<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<InnerSeries<C>> for &'a Series<C>
where
    Series<C>: Mul<&'a Series<C>, Output = Series<C>>,
    &'a Series<C>: Mul<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, C: Coeff> Mul<&'a C> for InnerSeries<C>
where
    C: Mul<&'a C, Output = C>,
    Series<C>: Mul<&'a C, Output = Series<C>>,
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
    InnerSeries<C>: Mul<&'a C, Output = Self> + Mul<&'a Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Mul<Series<C>> for &'a InnerSeries<C>
where
    Series<C>: Mul<&'a Series<C>, Output = Series<C>> + Mul<&'a C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.mul(s).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Mul<&'a InnerSeries<C>> for Series<C>
where
    Series<C>: Mul<&'a Series<C>, Output = Series<C>> + Mul<&'a C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, C: Coeff> Mul<C> for &'a InnerSeries<C>
where
    C: Mul<&'a C, Output = C>,
    &'a Series<C>: Mul<C, Output = Series<C>>,
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
    InnerSeries<C>: Mul<&'a C, Output = InnerSeries<C>> + Mul<&'a Series<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b Series<C>> for &'a InnerSeries<C>
where
    &'a Series<C>: Mul<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>: Mul<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'b Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul(rhs).into(),
            Inner(c) => rhs.mul(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Mul<&'a InnerSeries<C>> for &'b Series<C>
where
    &'a Series<C>: Mul<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>: Mul<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn mul(self, rhs: &'a InnerSeries<C>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<'a, 'b, C: Coeff> Mul<&'b C> for &'a InnerSeries<C>
where
    &'a C: Mul<&'b C, Output = C>,
    &'a Series<C>: Mul<&'b C, Output = Series<C>>,
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
    &'a InnerSeries<C>: Mul<&'b C, Output = InnerSeries<C>> + Mul<&'b Series<C>, Output = InnerSeries<C>>,
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

impl<C: Coeff + Default> MulAssign<Series<C>> for InnerSeries<C>
where
    Series<C>: MulAssign + MulAssign<C>
{
    fn mul_assign(&mut self, mut rhs: Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => {
                rhs.mul_assign(std::mem::take(c));
                *self = rhs.into();
            },
        }
    }
}

impl<C: Coeff + Default> MulAssign<C> for InnerSeries<C>
where
    Series<C>: MulAssign<C>,
    C: MulAssign
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
    InnerSeries<C>: MulAssign<Series<C>> + MulAssign<C>
{
    fn mul_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul_assign(s),
            Inner(c) => self.mul_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> MulAssign<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: MulAssign<&'a Series<C>>,
    &'a Series<C>: Mul<C, Output = Series<C>>,
{
    fn mul_assign(&mut self, rhs: &'a Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.mul_assign(rhs),
            Inner(c) => {
                let res = rhs.mul(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<'a, C: Coeff + Default> MulAssign<&'a C> for InnerSeries<C>
where
    Series<C>: MulAssign<&'a C>,
    C: MulAssign<&'a C>
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
    InnerSeries<C>: MulAssign<&'a Series<C>> + MulAssign<&'a C>
{
    fn mul_assign(&mut self, rhs: &'a InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.mul_assign(s),
            Inner(c) => self.mul_assign(c),
        }
    }
}

impl<C: Coeff> Div<Series<C>> for InnerSeries<C>
where
    Series<C>: Div<Output = Series<C>> + Mul<C, Output = Series<C>> + MulInverse<Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<C: Coeff> Div<InnerSeries<C>> for Series<C>
where
    Series<C>: Div<Output = Series<C>> + Div<C, Output = Series<C>>,
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
    Series<C>: Div<C, Output = Series<C>>,
    C: Div<Output = C>
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
    InnerSeries<C>: Div<C, Output = Self> + Div<Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Div<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: Div<&'a Series<C>, Output = Series<C>> + MulInverse<Output = Self>,
    &'a Series<C>: Div<C, Output = Series<C>>
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'a Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.div(c).mul_inverse().into(),
        }
    }
}

impl<'a, C: Coeff> Div<InnerSeries<C>> for &'a Series<C>
where
    InnerSeries<C>: MulInverse<Output = InnerSeries<C>> + Mul<&'a Series<C>, Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul_inverse().mul(self)
    }
}

impl<'a, C: Coeff> Div<&'a C> for InnerSeries<C>
where
    C: Div<&'a C, Output = C>,
    Series<C>: Div<&'a C, Output = Series<C>>,
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
    InnerSeries<C>: Div<&'a C, Output = Self> + Div<&'a Series<C>, Output = Self>,
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

impl<'a, C: Coeff> Div<Series<C>> for &'a InnerSeries<C>
where
    Series<C>: Mul<&'a Series<C>, Output = Series<C>> + Mul<&'a C, Output = Series<C>> + MulInverse<Output = Series<C>>,
C: MulInverse<Output = C>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => rhs.mul_inverse().mul(s).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<'a, C: Coeff> Div<&'a InnerSeries<C>> for Series<C>
where
    Series<C>: Mul<InnerSeries<C>, Output = InnerSeries<C>>,
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
    &'a Series<C>: Div<C, Output = Series<C>>,
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
    InnerSeries<C>: Mul<&'a C, Output = InnerSeries<C>> + Mul<&'a Series<C>, Output = InnerSeries<C>> + MulInverse<Output = InnerSeries<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: InnerSeries<C>) -> Self::Output {
        rhs.mul_inverse().mul(self)
    }
}

impl<'a, 'b, C: Coeff> Div<&'b Series<C>> for &'a InnerSeries<C>
where
    &'a Series<C>: Div<&'b Series<C>, Output = Series<C>>,
    &'b Series<C>:  MulInverse<Output = Series<C>>,
    Series<C>: Mul<&'a C, Output = Series<C>>,
{
    type Output = InnerSeries<C>;

    fn div(self, rhs: &'b Series<C>) -> Self::Output {
        use InnerSeries::*;
        match self {
            Series(s) => s.div(rhs).into(),
            Inner(c) => rhs.mul_inverse().mul(c).into(),
        }
    }
}

impl<'a, 'b, C: Coeff> Div<&'a InnerSeries<C>> for &'b Series<C>
where
    &'b Series<C>: Div<&'a C, Output = Series<C>> + Div<&'a Series<C>, Output = Series<C>>,
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
    &'a Series<C>: Div<&'b C, Output = Series<C>>,
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
    &'a InnerSeries<C>: Div<&'b C, Output = InnerSeries<C>> + Div<&'b Series<C>, Output = InnerSeries<C>>,
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

impl<C: Coeff + Default> DivAssign<Series<C>> for InnerSeries<C>
where
    Series<C>: DivAssign + MulAssign<C> + MulInverse<Output = Series<C>>
{
    fn div_assign(&mut self, rhs: Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => {
                let mut res = rhs.mul_inverse();
                res.mul_assign(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<C: Coeff + Default> DivAssign<C> for InnerSeries<C>
where
    Series<C>: DivAssign<C>,
    C: DivAssign
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
    InnerSeries<C>: DivAssign<Series<C>> + DivAssign<C>
{
    fn div_assign(&mut self, rhs: InnerSeries<C>) {
        use InnerSeries::*;
        match rhs {
            Series(s) => self.div_assign(s),
            Inner(c) => self.div_assign(c),
        }
    }
}

impl<'a, C: Coeff + Default> DivAssign<&'a Series<C>> for InnerSeries<C>
where
    Series<C>: MulAssign<C> + DivAssign<&'a Series<C>>,
    &'a Series<C>: MulInverse<Output = Series<C>>,
{
    fn div_assign(&mut self, rhs: &'a Series<C>) {
        use InnerSeries::*;
        match self {
            Series(s) => s.div_assign(rhs),
            Inner(c) => {
                let mut res = rhs.mul_inverse();
                res.mul_assign(std::mem::take(c));
                *self = res.into();
            },
        }
    }
}

impl<'a, C: Coeff + Default> DivAssign<&'a C> for InnerSeries<C>
where
    Series<C>: DivAssign<&'a C>,
    C: DivAssign<&'a C>
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
    InnerSeries<C>: DivAssign<&'a Series<C>> + DivAssign<&'a C>
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
    InnerSeries<C>: Add<Output = Self>
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
    InnerSeries<C>: Mul<Output = Self>
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
