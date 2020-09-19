//! Helper functions that poisson uses.

use crate::{Builder, Float, Type, Vector};

use num_traits::{Float as NumFloat, NumCast};

use rand::distributions::{Distribution, Standard};
use rand::Rng;

use modulo::Mod;

use std::marker::PhantomData;

pub mod math;

#[derive(Clone)]
pub struct Grid<F, V>
where
    F: Float,
    V: Vector<F>,
{
    data: Vec<Vec<V>>,
    side: usize,
    cell: F,
    poisson_type: Type,
    _marker: PhantomData<F>,
}

impl<F, V> Grid<F, V>
where
    F: Float,
    V: Vector<F>,
{
    pub fn new(radius: F, poisson_type: Type) -> Grid<F, V> {
        let dim = F::cast(V::dimension());
        let cell = (F::cast(2) * radius) / NumFloat::sqrt(dim);
        let side = (F::cast(1) / cell)
            .to_usize()
            .expect("Expected that dividing 1 by cell width would be legal.");
        Grid {
            cell: cell,
            side: side,
            data: vec![
                vec![];
                side.pow(
                    dim.to_u32()
                        .expect("Dimension should be always be castable to u32.")
                )
            ],
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

pub fn encode<F, V>(v: &V, side: usize, poisson_type: Type) -> Option<usize>
where
    F: Float,
    V: Vector<F>,
{
    use crate::Type::*;
    let mut index = 0;
    for n in 0..V::dimension() {
        let n = v[n];
        let cur = match poisson_type {
            Perioditic => n
                .to_isize()
                .expect(
                    "Expected that all scalars of the index vector should be castable to \
                     isize.",
                )
                .modulo(side as isize) as usize,
            Normal => {
                if n < F::cast(0) || n >= F::cast(side) {
                    return None;
                }
                n.to_usize().expect(
                    "Expected that all scalars of the index vector should be castable to \
                     usize.",
                )
            }
        };
        index = (index + cur) * side;
    }
    Some(index / side)
}

pub fn decode<F, V>(index: usize, side: usize) -> Option<V>
where
    F: Float,
    V: Vector<F>,
{
    let dim = V::dimension();
    if index >= side.pow(dim as u32) {
        return None;
    }
    let mut result = V::zero();
    let mut last = index;
    for n in (0..V::dimension()).rev() {
        let cur = last / side;
        result[n] = F::cast(last - cur * side);
        last = cur;
    }
    Some(result)
}

#[test]
fn encoding_decoding_works() {
    let n = nalgebra::Vector2::new(10., 7.);
    assert_eq!(
        n,
        decode::<_, nalgebra::Vector2<_>>(encode(&n, 15, Type::Normal).unwrap(), 15).unwrap(),
    );
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = nalgebra::Vector2::new(14., 14.);
    assert_eq!(
        n,
        decode::<_, nalgebra::Vector2<_>>(encode(&n, 15, Type::Normal).unwrap(), 15).unwrap()
    );
}

#[test]
fn encoding_outside_of_area_fails() {
    let n = nalgebra::Vector2::new(9., 7.);
    assert_eq!(None, encode(&n, 9, Type::Normal));
    let n = nalgebra::Vector2::new(7., 9.);
    assert_eq!(None, encode(&n, 9, Type::Normal));
}

#[test]
fn decoding_outside_of_area_fails() {
    assert_eq!(None, decode::<f64, nalgebra::Vector2<_>>(100, 10));
}

pub fn choose_random_sample<F, V, R>(rng: &mut R, grid: &Grid<F, V>, index: V, level: usize) -> V
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    Standard: Distribution<V>,
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / F::cast(side);
    (index + rng.gen()) * spacing
}

#[test]
fn random_point_is_between_right_values_top_lvl() {
    use num_traits::Zero;
    use rand::{rngs::SmallRng, SeedableRng};
    let mut rand = SmallRng::from_seed([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]);
    let radius = 0.2;
    let grid = Grid::<f64, nalgebra::Vector2<_>>::new(radius, Type::Normal);
    for _ in 0..1000 {
        let result = choose_random_sample(&mut rand, &grid, nalgebra::Vector2::<f64>::zero(), 0);
        assert!(result.x >= 0.);
        assert!(result.x < grid.cell);
        assert!(result.y >= 0.);
        assert!(result.y < grid.cell);
    }
}

pub fn sample_to_index<F, V>(value: &V, side: usize) -> V
where
    F: Float,
    V: Vector<F>,
{
    let mut cur = value.clone();
    for n in 0..V::dimension() {
        cur[n] = NumFloat::floor(cur[n] * F::cast(side));
    }
    cur
}

pub fn index_to_sample<F, V>(value: &V, side: usize) -> V
where
    F: Float,
    V: Vector<F>,
{
    let mut cur = value.clone();
    for n in 0..V::dimension() {
        cur[n] = cur[n] / F::cast(side);
    }
    cur
}

pub fn is_disk_free<F, V>(
    grid: &Grid<F, V>,
    poisson: &Builder<F, V>,
    index: V,
    level: usize,
    sample: V,
    outside: &[V],
) -> bool
where
    F: Float,
    V: Vector<F>,
{
    let parent = get_parent(index, level);
    let sqradius = NumFloat::powi(F::cast(2) * poisson.radius, 2);
    // NOTE: This does unnessary checks for corners, but it doesn't affect much in higher dimensions: 5^d vs 5^d - 2d
    each_combination(&[-2, -1, 0, 1, 2])
        .filter_map(|t| grid.get(parent.clone() + t))
        .flat_map(|t| t)
        .all(|v| sqdist(v.clone(), sample.clone(), poisson.poisson_type) >= sqradius)
        && is_valid(poisson, outside, sample)
}

pub fn is_valid<F, V>(poisson: &Builder<F, V>, samples: &[V], sample: V) -> bool
where
    F: Float,
    V: Vector<F>,
{
    let sqradius = NumFloat::powi(F::cast(2) * poisson.radius, 2);
    samples
        .iter()
        .all(|t| sqdist(t.clone(), sample.clone(), poisson.poisson_type) >= sqradius)
}

pub fn sqdist<F, V>(v1: V, v2: V, poisson_type: Type) -> F
where
    F: Float,
    V: Vector<F>,
{
    use crate::Type::*;
    let diff = v2 - v1;
    match poisson_type {
        Perioditic => each_combination(&[-1, 0, 1])
            .map(|v| (diff.clone() + v).norm_squared())
            .fold(NumFloat::max_value(), |a, b| NumFloat::min(a, b)),
        Normal => diff.norm_squared(),
    }
}

pub fn get_parent<F, V>(mut index: V, level: usize) -> V
where
    F: Float,
    V: Vector<F>,
{
    let split = 2usize.pow(level as u32);
    for n in 0..V::dimension() {
        index[n] = NumFloat::floor(index[n] / F::cast(split));
    }
    index
}

#[test]
fn getting_parent_works() {
    let divides = 4;
    let cells_per_cell = 2usize.pow(divides as u32);
    let testee = nalgebra::Vector2::new(1., 2.);
    assert_eq!(
        testee,
        get_parent(
            (testee * cells_per_cell as f64) + nalgebra::Vector2::new(0., 15.),
            divides
        )
    );
}

pub struct CombiIter<'a, F, FF, V>
where
    F: Float,
    FF: NumCast + 'a,
    V: Vector<F>,
{
    cur: usize,
    choices: &'a [FF],
    _marker: PhantomData<(F, V)>,
}

impl<'a, F, FF, V> Iterator for CombiIter<'a, F, FF, V>
where
    F: Float,
    FF: NumCast + Clone,
    V: Vector<F>,
{
    type Item = V;
    fn next(&mut self) -> Option<Self::Item> {
        let dim = V::dimension();
        let len = self.choices.len();
        if self.cur >= len.pow(dim as u32) {
            None
        } else {
            let mut result = V::zero();
            let mut div = self.cur;
            self.cur += 1;
            for n in 0..V::dimension() {
                let rem = div % len;
                div /= len;
                let choice = self.choices[rem as usize].clone();
                result[n] = NumCast::from(choice).expect(
                    "Expected that all choices were castable to float without \
                     problems.",
                );
            }
            Some(result)
        }
    }
}

/// Iterates through all combinations of vectors with allowed values as scalars.
pub fn each_combination<'a, F, FF, V>(choices: &[FF]) -> CombiIter<F, FF, V>
where
    F: Float + 'a,
    FF: NumCast,
    V: Vector<F>,
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
    where
        I: IntoIterator<Item = T>,
        F: FnMut(T) -> I;
}

impl<T> Inplace<T> for Vec<T> {
    fn flat_map_inplace<F, I>(&mut self, mut f: F)
    where
        I: IntoIterator<Item = T>,
        F: FnMut(T) -> I,
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
