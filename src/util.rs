// removes `elem` from the start of `v` and returns the number of elements removed
pub(crate) fn trim_start<T>(v: &mut Vec<T>, elem: &T) -> usize
where
    for<'a> &'a T: PartialEq,
{
    let first_ok = v.iter().enumerate().find(|&(_, c)| c != elem);
    match first_ok {
        Some((idx, _)) => {
            v.drain(0..idx);
            idx
        }
        None => {
            let orig_len = v.len();
            v.clear();
            orig_len
        }
    }
}

// removes `elem` from the end of `v` and returns the number of elements removed
pub(crate) fn trim_end<T>(v: &mut Vec<T>, elem: &T) -> usize
where
    for<'a> &'a T: PartialEq,
{
    let last_ok = v.iter().rev().enumerate().find(|&(_, c)| c != elem);
    match last_ok {
        Some((idx, _)) => {
            let new_len = v.len() - idx;
            v.truncate(new_len);
            idx
        }
        None => {
            let orig_len = v.len();
            v.clear();
            orig_len
        }
    }
}

// returns a slice with leading occurences of `elem` removed and the
// number of removed elements
pub(crate) fn trim_slice_start<'a, 'b, T>(
    s: &'a [T],
    elem: &'b T,
) -> (&'a [T], usize)
where
    for<'c> &'c T: PartialEq,
{
    let remove = s
        .iter()
        .enumerate()
        .find(|&(_, c)| c != elem)
        .map(|(idx, _)| idx)
        .unwrap_or_else(|| s.len());
    (&s[remove..], remove)
}

// returns a slice with trailing occurences of `elem` removed and the
// number of removed elements
pub(crate) fn trim_slice_end<'a, 'b, T>(
    s: &'a [T],
    elem: &'b T,
) -> (&'a [T], usize)
where
    for<'c> &'c T: PartialEq,
{
    let remove = s
        .iter()
        .rev()
        .enumerate()
        .find(|&(_, c)| c != elem)
        .map(|(idx, _)| idx)
        .unwrap_or_else(|| s.len());
    (&s[..(s.len() - remove)], remove)
}
