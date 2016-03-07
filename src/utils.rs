use ::{PoissonType, VecLike};

use rand::{Rand, Rng};

use modulo::Mod;

use std::marker::PhantomData;
use std::f64;

#[derive(Clone)]
pub struct Grid<V>
    where V: VecLike
{
    pub data: Vec<Vec<V>>,
    pub side: usize,
    pub cell: f64,
    poisson_type: PoissonType,
}

impl<V> Grid<V> where V: VecLike
{
    pub fn new(radius: f64, poisson_type: PoissonType) -> Grid<V> {
        let dim = V::dim(None);
        let cell = (2. * radius) / (dim as f64).sqrt();
        let side = (1. / cell) as usize;
        Grid {
            cell: cell,
            side: side,
            data: vec![vec![]; side.pow(dim as u32)],
            poisson_type: poisson_type,
        }
    }

    pub fn get(&self, index: V) -> Option<&Vec<V>> {
        encode(&index, self.side, self.poisson_type).map(|t| &self.data[t])
    }

    pub fn get_mut(&mut self, index: V) -> Option<&mut Vec<V>> {
        encode(&index, self.side, self.poisson_type).map(move |t| &mut self.data[t])
    }

    pub fn cells(&self) -> usize {
        self.data.len()
    }
}

pub fn encode<V>(v: &V, side: usize, poisson_type: PoissonType) -> Option<usize>
    where V: VecLike
{
    use PoissonType::*;
    let mut index = 0;
    for &n in v.iter() {
        let mut cur = n as usize;
        match poisson_type {
            Perioditic => cur = (n as isize).modulo(side as isize) as usize,
            Normal => if n < 0. || n >= side as f64 {
                        return None;
                    },
        }
        index = (index + cur) * side;
    }
    Some(index / side)
}

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
        *n = (last - cur * side) as f64;
        last = cur;
    }
    Some(result)
}

#[test]
fn encoding_decoding_works() {
    let n = ::na::Vec2::new(10., 7.);
    assert_eq!(n, decode(encode(&n, 15, PoissonType::Normal).unwrap(), 15).unwrap());
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = ::na::Vec2::new(14., 14.);
    assert_eq!(n, decode(encode(&n, 15, PoissonType::Normal).unwrap(), 15).unwrap());
}

#[test]
fn encoding_outside_of_area_fails() {
    let n = ::na::Vec2::new(9., 7.);
    assert_eq!(None, encode(&n, 9, PoissonType::Normal));
    let n = ::na::Vec2::new(7., 9.);
    assert_eq!(None, encode(&n, 9, PoissonType::Normal));
}

#[test]
fn decoding_outside_of_area_fails() {
    assert_eq!(None, decode::<::na::Vec2<f64>>(100, 10));
}

pub fn choose_random_sample<V, R>(rand: &mut R, grid: &Grid<V>, index: V, level: usize) -> V
    where V: VecLike,
          R: Rng
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / side as f64;
    let mut result = index * spacing;
    for n in result.iter_mut() {
        let place = f64::rand(rand);
        *n += place * spacing;
    }
    result
}

#[test]
fn random_point_is_between_right_values_top_lvl() {
    use num::Zero;
    use rand::{SeedableRng, XorShiftRng};
    use ::na::Vec2;
    let mut rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let radius = 0.2;
    let grid = Grid::<Vec2<f64>>::new(radius, PoissonType::Normal);
    for _ in 0..1000 {
        let result = choose_random_sample(&mut rand, &grid, Vec2::<f64>::zero(), 0);
        assert!(result.x >= 0.);
        assert!(result.x < grid.cell);
        assert!(result.y >= 0.);
        assert!(result.y < grid.cell);
    }
}

pub fn sample_to_index<V>(value: &V, side: usize) -> V
    where V: VecLike
{
    let mut cur = value.clone();
    for c in cur.iter_mut() {
        *c = (*c * side as f64).floor();
    }
    cur
}

pub fn is_disk_free<V>(grid: &Grid<V>,
                      radius: f64,
                      poisson_type: PoissonType,
                      index: V,
                      level: usize,
                      c: V)
                      -> bool
    where V: VecLike
{
    let parent = get_parent(index, level);
    let sqradius = (2. * radius).powi(2);
    // TODO: This does unnessary checks at corners...
    each_combination(&[-2., -1., 0., 1., 2.])
        .filter_map(|t| grid.get(parent.clone() + t))
        .flat_map(|t| t)
        .all(|v| sqdist(v.clone(), c.clone(), poisson_type) >= sqradius)
}

pub fn sqdist<V>(v1: V, v2: V, poisson_type: PoissonType) -> f64
    where V: VecLike
{
    use PoissonType::*;
    let diff = v2 - v1;
    match poisson_type {
        Perioditic =>
            each_combination(&[-1., 0., 1.])
                .map(|v| (diff.clone() + v).sqnorm())
                .fold(f64::MAX, |a, b| a.min(b)),
        Normal => diff.sqnorm(),
    }
}

pub fn get_parent<V>(mut index: V, level: usize) -> V
    where V: VecLike
{
    let split = 2usize.pow(level as u32);
    for n in index.iter_mut() {
        *n = (*n / split as f64).floor();
    }
    index
}

#[test]
fn getting_parent_works() {
    let divides = 4;
    let cells_per_cell = 2usize.pow(divides as u32);
    let testee = ::na::Vec2::new(1., 2.);
    assert_eq!(testee,
               get_parent((testee * cells_per_cell as f64) + ::na::Vec2::new(0., 15.),
                          divides));
}

pub fn is_valid<V>(radius: f64, poisson_type: PoissonType, samples: &[V], sample: V) -> bool
    where V: VecLike
{
    let sqradius = (2. * radius).powi(2);
    samples
        .iter()
        .all(|t| sqdist(t.clone(), sample.clone(), poisson_type) >= sqradius)
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
                *n = self.choices[rem as usize];
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
