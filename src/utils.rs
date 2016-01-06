use VecLike;

use std::marker::PhantomData;

pub struct CombiIter<'a, V>
    where V: VecLike
{
    cur: usize,
    choices: &'a [f64],
    _marker: PhantomData<V>,
}

impl<'a, V> Iterator for CombiIter<'a, V> where V: VecLike
{
    type Item = V;
    fn next(&mut self) -> Option<Self::Item> {
        let dim = V::dim(None);
        let len = self.choices.len();
        if self.cur >= len.pow(dim as u32) {
            None
        } else {
            let mut result = V::zero();
            let mut div = self.cur;
            self.cur += 1;
            for n in 0..dim {
                let rem = div % len;
                div /= len;
                result[n] = self.choices[rem as usize];
            }
            Some(result)
        }
    }
}

pub fn each_combination<'a, V>(choices: &[f64]) -> CombiIter<V>
    where V: VecLike
{
    CombiIter {
        cur: 0,
        choices: choices,
        _marker: PhantomData,
    }
}

pub trait Inplace<T> {
    fn flat_map_inplace<F, I>(&mut self, f: F)
        where I: IntoIterator<Item = T>,
              F: FnMut(T) -> I;
}

impl<T> Inplace<T> for Vec<T> {
    fn flat_map_inplace<F, I>(&mut self, mut f: F)
        where I: IntoIterator<Item = T>,
              F: FnMut(T) -> I
    {
        for i in (0..self.len()).rev() {
            for t in f(self.swap_remove(i)) {
                self.push(t);
            }
        }
    }
}

#[test]
fn mapping_inplace_works() {
    let vec = vec![1, 2, 3, 4, 5, 6];
    let mut result = vec.clone();
    let func = |t| {
        match t % 3 {
            0 => (0..0),
            1 => (0..1),
            _ => (0..2),
        }
        .map(move |n| t + n)
    };
    result.flat_map_inplace(&func);
    let mut expected = vec.into_iter().flat_map(func).collect::<Vec<_>>();
    assert_eq!(expected.sort(), result.sort());
}
