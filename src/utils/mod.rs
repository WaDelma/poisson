//! Helper functions that poisson uses.

use {PoissonDisk, PoissonType, VecLike, FloatLike};

use num::{NumCast, Float};

use rand::{Rand, Rng};

use modulo::Mod;

use std::marker::PhantomData;

pub mod math;

#[derive(Clone)]
pub struct Grid<F, V>
    where F: FloatLike,
          V: VecLike<F>
{
    data: Vec<Vec<V>>,
    side: usize,
    cell: F,
    poisson_type: PoissonType,
    _marker: PhantomData<F>,
}

impl<F, V> Grid<F, V>
    where F: FloatLike,
          V: VecLike<F>
{
    pub fn new(radius: F, poisson_type: PoissonType) -> Grid<F, V> {
        let dim = F::f(V::dim(None));
        let cell = (F::f(2) * radius) / dim.sqrt();
        let side = (F::f(1) / cell)
                       .to_usize()
                       .expect("Expected that diving 1 by cell width would be legal.");
        Grid {
            cell: cell,
            side: side,
            data: vec![vec![]; side.pow(dim.to_u32().expect("Dimension should be always be castable to u32."))],
            poisson_type: poisson_type,
            _marker: PhantomData,
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

    pub fn side(&self) -> usize {
        self.side
    }

    pub fn cell(&self) -> F {
        self.cell
    }
}

pub fn encode<F, V>(v: &V, side: usize, poisson_type: PoissonType) -> Option<usize>
    where F: FloatLike,
          V: VecLike<F>
{
    use PoissonType::*;
    let mut index = 0;
    for &n in v.iter() {
        let cur = match poisson_type {
            Perioditic => {
                n.to_isize()
                 .expect("Expected that all scalars of the index vector should be castable to \
                          isize.")
                 .modulo(side as isize) as usize
            }
            Normal => {
                if n < F::f(0) || n >= F::f(side) {
                    return None;
                }
                n.to_usize()
                 .expect("Expected that all scalars of the index vector should be castable to \
                          usize.")
            }
        };
        index = (index + cur) * side;
    }
    Some(index / side)
}

pub fn decode<F, V>(index: usize, side: usize) -> Option<V>
    where F: FloatLike,
          V: VecLike<F>
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
        *n = F::f(last - cur * side);
        last = cur;
    }
    Some(result)
}

#[test]
fn encoding_decoding_works() {
    let n = ::na::Vec2::new(10., 7.);
    assert_eq!(n,
               decode(encode(&n, 15, PoissonType::Normal).unwrap(), 15).unwrap());
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = ::na::Vec2::new(14., 14.);
    assert_eq!(n,
               decode(encode(&n, 15, PoissonType::Normal).unwrap(), 15).unwrap());
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
    assert_eq!(None, decode::<f64, ::na::Vec2<_>>(100, 10));
}

pub fn choose_random_sample<F, V, R>(rng: &mut R, grid: &Grid<F, V>, index: V, level: usize) -> V
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / F::f(side);
    (index + V::rand(rng)) * spacing
}

#[test]
fn random_point_is_between_right_values_top_lvl() {
    use num::Zero;
    use rand::{SeedableRng, XorShiftRng};
    use na::Vec2;
    let mut rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let radius = 0.2;
    let grid = Grid::<f64, ::na::Vec2<_>>::new(radius, PoissonType::Normal);
    for _ in 0..1000 {
        let result = choose_random_sample(&mut rand, &grid, Vec2::<f64>::zero(), 0);
        assert!(result.x >= 0.);
        assert!(result.x < grid.cell);
        assert!(result.y >= 0.);
        assert!(result.y < grid.cell);
    }
}

pub fn sample_to_index<F, V>(value: &V, side: usize) -> V
    where F: FloatLike,
          V: VecLike<F>
{
    let mut cur = value.clone();
    for c in cur.iter_mut() {
        *c = (*c * F::f(side)).floor();
    }
    cur
}

pub fn is_disk_free<F, V>(grid: &Grid<F, V>,
                          poisson: &PoissonDisk<F, V>,
                          index: V,
                          level: usize,
                          sample: V,
                          outside: &[V])
                          -> bool
    where F: FloatLike,
          V: VecLike<F>
{
    let parent = get_parent(index, level);
    let sqradius = (F::f(2) * poisson.radius).powi(2);
    // NOTE: This does unnessary checks for corners, but it doesn't affect much in higher dimensions: 5^d vs 5^d - 2d
    each_combination(&[-2, -1, 0, 1, 2])
        .filter_map(|t| grid.get(parent.clone() + t))
        .flat_map(|t| t)
        .all(|v| sqdist(v.clone(), sample.clone(), poisson.poisson_type) >= sqradius) &&
    is_valid(poisson, outside, sample)
}

pub fn is_valid<F, V>(poisson: &PoissonDisk<F, V>, samples: &[V], sample: V) -> bool
    where F: FloatLike,
          V: VecLike<F>
{
    let sqradius = (F::f(2) * poisson.radius).powi(2);
    samples.iter()
           .all(|t| sqdist(t.clone(), sample.clone(), poisson.poisson_type) >= sqradius)
}


pub fn sqdist<F, V>(v1: V, v2: V, poisson_type: PoissonType) -> F
    where F: FloatLike,
          V: VecLike<F>
{
    use PoissonType::*;
    let diff = v2 - v1;
    match poisson_type {
        Perioditic => {
            each_combination(&[-1, 0, 1])
                .map(|v| (diff.clone() + v).sqnorm())
                .fold(F::max_value(), |a, b| a.min(b))
        }
        Normal => diff.sqnorm(),
    }
}

pub fn get_parent<F, V>(mut index: V, level: usize) -> V
    where F: FloatLike,
          V: VecLike<F>
{
    let split = 2usize.pow(level as u32);
    for n in index.iter_mut() {
        *n = (*n / F::f(split)).floor();
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

pub struct CombiIter<'a, F, FF, V>
    where F: FloatLike,
          FF: NumCast + 'a,
          V: VecLike<F>
{
    cur: usize,
    choices: &'a [FF],
    _marker: PhantomData<(F, V)>,
}

impl<'a, F, FF, V> Iterator for CombiIter<'a, F, FF, V>
    where F: FloatLike,
          FF: NumCast + Clone,
          V: VecLike<F>
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
                let choice = self.choices[rem as usize].clone();
                *n = NumCast::from(choice)
                         .expect("Expected that all choices were castable to float without \
                                  problems.");
            }
            Some(result)
        }
    }
}

/// Iterates through all combinations of vectors with allowed values as scalars.
pub fn each_combination<'a, F, FF, V>(choices: &[FF]) -> CombiIter<F, FF, V>
    where F: FloatLike + 'a,
          FF: NumCast,
          V: VecLike<F>
{
    CombiIter {
        cur: 0,
        choices: choices,
        _marker: PhantomData,
    }
}

/// Trait that allows flat mapping inplace.
pub trait Inplace<T> {
    /// Does flat map inplace without maintaining order of elements.
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
