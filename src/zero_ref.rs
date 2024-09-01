use std::{any::{Any, TypeId}, sync::LazyLock};

use num_traits::Zero;

pub(crate) fn zero_ref<C: Sync + Send + Zero>() -> &'static C
{
    type TypeMap = elsa::sync::FrozenMap<TypeId, &'static (dyn Any + Sync + Send)>;
    static ZERO_MAP: LazyLock<TypeMap> = LazyLock::new(|| TypeMap::new());

    let res = ZERO_MAP.get(&TypeId::of::<C>());

    if let Some(res) = res {
        return res.downcast_ref().unwrap();
    }

    let res = Box::leak(Box::new(C::zero()));
    ZERO_MAP.insert(TypeId::of::<C>(), res);
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        assert!(zero_ref::<i32>().is_zero());
        assert!(zero_ref::<f64>().is_zero());
    }
}
