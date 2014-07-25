use super::{stride};

#[test]
fn test_every_other()
{
    let count = vec![1i, 2i, 3i, 4i, 5i];
    let odds = stride(count.as_slice(), 2);
    let odds_should_be = vec![1i, 3i, 5i];
    let evens = stride(count.slice_from(1), 2);
    let evens_should_be = vec![2i, 4i];
    for (a, b) in odds_should_be.iter().zip(odds)
    {
        assert_eq!(a, b);
    }
    for (a, b) in evens_should_be.iter().zip(evens)
    {
        assert_eq!(a, b);
    }
}
