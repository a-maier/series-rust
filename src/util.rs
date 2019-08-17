pub(crate) fn trim_start<T>(v: &mut Vec<T>, elem: &T) -> usize
where for<'a> &'a T: PartialEq
{
    let first_ok = v.iter()
        .enumerate()
        .find(|&(_, c)| c != elem);
    match first_ok {
        Some((idx, _)) => {
            if idx > 0 {
                v.drain(0..idx);
                idx
            } else {
                0
            }
        },
        None => {
            let orig_len = v.len();
            v.clear();
            orig_len
        }
    }
}
