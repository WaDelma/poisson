use VecLike;

use modulo::Mod;

use std::marker::PhantomData;
use std::mem::replace;

#[derive(Clone)]
pub struct Grid<V>
    where V: VecLike
{
    pub data: Vec<Vec<V>>,
    pub side: usize,
    pub cell: f64,
    periodicity: bool,
}

impl<V> Grid<V> where V: VecLike
{
    pub fn new(radius: f64, periodicity: bool) -> Grid<V> {
        let dim = V::dim(None);
        let cell = (2. * radius) / (dim as f64).sqrt();
        let side = (1. / cell) as usize;
        Grid {
            cell: cell,
            side: side,
            data: vec![vec![]; side.pow(dim as u32)],
            periodicity: periodicity,
        }
    }

    pub fn get(&self, index: V) -> Option<&Vec<V>> {
        encode(&index, self.side, self.periodicity).map(|t| &self.data[t])
    }

    pub fn get_mut(&mut self, index: V) -> Option<&mut Vec<V>> {
        encode(&index, self.side, self.periodicity).map(move |t| &mut self.data[t])
    }

    pub fn cells(&self) -> usize {
        self.data.len()
    }
}

pub fn encode<V>(v: &V, side: usize, periodicity: bool) -> Option<usize>
    where V: VecLike
{
    let mut index = 0;
    for &n in v.iter() {
        let mut cur = n as usize;
        if periodicity {
            cur = (n as isize).modulo(side as isize) as usize;
        } else if n < 0. || n >= side as f64 {
            return None;
        }
        index = (index + cur) * side;
    }
    Some(index / side)
}

#[cfg(test)]
pub fn decode<V>(index: usize, side: usize) -> Option<V>
    where V: VecLike
{
    use num::Zero;
    let dim = V::dim(None);
    if index >= side.pow(dim as u32) {
        return None;
    }
    let mut result = V::zero();
    let mut last = index;
    for n in result.iter_mut().rev() {
        let cur = last / side;
        let value = (last - cur * side) as f64;
        replace(n, value);
        last = cur;
    }
    Some(result)
}

#[test]
fn encoding_decoding_works() {
    let n = ::na::Vec2::new(10., 7.);
    assert_eq!(n, decode(encode(&n, 15, false).unwrap(), 15).unwrap());
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = ::na::Vec2::new(14., 14.);
    assert_eq!(n, decode(encode(&n, 15, false).unwrap(), 15).unwrap());
}

#[test]
fn encoding_outside_of_area_fails() {
    let n = ::na::Vec2::new(9., 7.);
    assert_eq!(None, encode(&n, 9, false));
    let n = ::na::Vec2::new(7., 9.);
    assert_eq!(None, encode(&n, 9, false));
}

#[test]
fn decoding_outside_of_area_fails() {
    assert_eq!(None, decode::<::na::Vec2<f64>>(100, 10));
}

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
            for n in result.iter_mut() {
                let rem = div % len;
                div /= len;
                replace(n, self.choices[rem as usize]);
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

pub fn calculate<T, F>(value: &mut T, fun: F)
    where F: FnOnce(&mut T) -> T
{
    let v = (fun)(value);
    ::std::mem::replace(value, v);
}
