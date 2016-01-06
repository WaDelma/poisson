//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//! * For each point there is disk of certain radius which doesn't intersect
//! with any other disk of other points
//! * Samples fill the space uniformly
//!
extern crate modulo;
use modulo::Mod;

extern crate sphere;
use sphere::sphere_volume;

extern crate rand;
use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

extern crate num;
use num::{Zero, One};

extern crate nalgebra as na;
use na::{Dim, Norm};

#[macro_use]
extern crate lazy_static;

use std::cmp::PartialEq;
use std::ops::{Sub, Mul, Add, Div, IndexMut};
use std::marker::PhantomData;
use std::mem::replace;
use std::f64;

use utils::{each_combination, Inplace};

mod math;
mod utils;

/// Describes what traits the algorithm needs to be able to work.
pub trait VecLike:
    IndexMut<usize, Output = f64> +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<f64, Output = Self> +
    Div<f64, Output = Self> +
    Norm<f64> +
    PartialEq +
    Zero +
    One +
    Dim +
    Copy {}
impl<T> VecLike for T where T:
    IndexMut<usize, Output = f64> +
    Add<Output = T> +
    Sub<Output = T> +
    Mul<f64, Output = T> +
    Div<f64, Output = T> +
    Norm<f64> +
    PartialEq +
    Zero +
    One +
    Dim +
    Copy {}

/// Builds PoissonGen with wanted properties.
///
/// # Examples
///
/// ```¨
/// extern crate poisson;
/// extern crate rand;
/// extern crate nalgebra as na;
/// type Vec2 = na::Vec2<f64>;
/// use poisson::PoissonDisk;
///
/// fn main() {
///     let mut poisson = PoissonDisk::new(rand::weak_rng()).build_radius::<Vec2>(0.1);
///     let vecs = poisson.generate();
///     println!("{:?}", vecs);
/// }
/// ```
pub struct PoissonDisk<R>
    where R: Rng
{
    rand: R,
    periodicity: bool,
}

impl<R> PoissonDisk<R> where R: Rng
{
    /// Creates new Poisson-disk generator builder with random generator specified.
    pub fn new(rand: R) -> Self {
        PoissonDisk {
            rand: rand,
            periodicity: false,
        }
    }

    /// Sets the generator to generate perioditic Poisson-disk distributions.
    pub fn perioditic(mut self) -> Self {
        self.periodicity = true;
        self
    }

    /// Builds the generator with relative radius specified.
    /// Radius should be ]0, 1]
    pub fn build_relative_radius<V>(self, radius: f64) -> PoissonGen<R, V>
        where V: VecLike
    {
        assert!(0. < radius);
        assert!(radius <= 1.);
        PoissonGen::new(self.rand, radius * (2f64.sqrt() / 2.), self.periodicity)
    }

    /// Builds the generator with radius specified.
    /// Radius should be ]0, √2 / 2]
    pub fn build_radius<V>(self, radius: f64) -> PoissonGen<R, V>
        where V: VecLike
    {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        PoissonGen::new(self.rand, radius, self.periodicity)
    }

    /// Builds the generator with radius calculated so that approximately specified number of samples are generated.
    /// Amount of samples should be larger than 0.
    /// Relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    pub fn build_samples<V>(self, samples: u32, relative_radius: f64) -> PoissonGen<R, V>
        where V: VecLike
    {
        assert!(self.periodicity || V::dim(None) < 5);
        assert!(samples > 0);
        assert!(relative_radius >= 0.);
        assert!(relative_radius <= 1.);
        PoissonGen::new(self.rand,
                        math::calc_radius::<V>(samples, relative_radius, self.periodicity),
                        self.periodicity)
    }
}

struct Grid<V>
    where V: VecLike
{
    data: Vec<Option<V>>,
    side: usize,
    cell: f64,
    periodicity: bool,
}

impl<V> Grid<V> where V: VecLike
{
    fn new(radius: f64, periodicity: bool) -> Grid<V> {
        let dim = V::dim(None);
        let cell = (2. * radius) / (dim as f64).sqrt();
        let side = (1. / cell) as usize;
        Grid {
            cell: cell,
            side: side,
            data: vec![None; side.pow(dim as u32)],
            periodicity: periodicity,
        }
    }

    fn get(&self, index: V) -> Option<&Option<V>> {
        encode(&index, self.side, self.periodicity).map(|t| &self.data[t])
    }

    fn get_mut(&mut self, index: V) -> Option<&mut Option<V>> {
        encode(&index, self.side, self.periodicity).map(move |t| &mut self.data[t])
    }

    fn cells(&self) -> usize {
        self.data.len()
    }
}

/// Generates Poisson-disk distribution in [0, 1]² area with O(N) time and space complexity relative to the number of samples generated.
/// Based on Ebeida, Mohamed S., et al. "A Simple Algorithm for Maximal Poisson‐Disk Sampling in High Dimensions." Computer Graphics Forum. Vol. 31. No. 2pt4. Blackwell Publishing Ltd, 2012.
pub struct PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    dim: PhantomData<V>,
    rand: R,
    radius: f64,
    periodicity: bool,
    a: f64,
}

impl<R, V> PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    /// Sets the radius of the generator.
    pub fn set_radius(&mut self, radius: f64) {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        self.radius = radius;
    }

    /// Gets the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Generates a vector with Poisson-disk distribution in area [0, 1]²
    pub fn generate(&mut self) -> Vec<V> {
        self.into_iter().collect()
    }
}

impl<R, V> PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    fn new(rand: R, radius: f64, periodicity: bool) -> PoissonGen<R, V> {
        let dim = V::dim(None);
        PoissonGen {
            dim: PhantomData,
            radius: radius,
            rand: rand,
            periodicity: periodicity,
            a: match dim {
                2 => 0.3,
                3 => 0.3,
                4 => 0.6,
                5 => 10.,
                6 => 700.,
                // TODO: Figure out what are optimal values beyond 6 dimensions
                _ => 700. + 100. * dim as f64,
            },
        }
    }
}

impl<'a, R, V> IntoIterator for &'a mut PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    type IntoIter = PoissonIter<'a, R, V>;
    type Item = V;
    fn into_iter(self) -> PoissonIter<'a, R, V> {
        let dim = V::dim(None);
        let grid = Grid::new(self.radius, self.periodicity);
        let capacity = grid.cells() * dim;
        let mut indices = Vec::with_capacity(capacity);
        let choices = (0..grid.side).map(|i| i as f64).collect::<Vec<_>>();
        indices.extend(each_combination::<V>(&choices));
        let range = Range::new(0, indices.len());
        let throws = (self.a * indices.len() as f64).ceil() as usize;
        PoissonIter {
            poisson: self,
            grid: grid,
            indices: indices,
            level: 0,
            range: range,
            throws: throws,
            success: 0,
        }
    }
}

pub struct PoissonIter<'a, R, V>
    where R: Rng + 'a,
          V: VecLike + 'a
{
    poisson: &'a mut PoissonGen<R, V>,
    grid: Grid<V>,
    indices: Vec<V>,
    level: usize,
    range: Range<usize>,
    throws: usize,
    success: usize,
}

impl<'a, R, V> Iterator for PoissonIter<'a, R, V>
    where R: Rng,
          V: VecLike
{
    type Item = V;
    fn next(&mut self) -> Option<V> {
        while !self.indices.is_empty() && self.level < f64::MANTISSA_DIGITS as usize {
            while self.throws > 0 {
                self.throws -= 1;
                let index = self.range.ind_sample(&mut self.poisson.rand);
                let cur = self.indices[index];
                let parent = get_parent(cur, self.level);
                if self.grid
                       .get(parent)
                       .expect("Indexing base grid by valid parent failed.")
                       .is_some() {
                    self.indices.swap_remove(index);
                    if self.indices.is_empty() {
                        return None;
                    }
                    self.range = Range::new(0, self.indices.len());
                } else {
                    let sample = choose_random_sample(&mut self.poisson.rand,
                                                      &self.grid,
                                                      cur,
                                                      self.level);
                    if is_disk_free(&self.grid, &self.poisson, cur, self.level, sample) {
                        replace(self.grid
                                    .get_mut(parent)
                                    .expect("Indexing base grid by already indexed valid parent \
                                             failed."),
                                Some(sample));
                        self.indices.swap_remove(index);
                        if !self.indices.is_empty() {
                            self.range = Range::new(0, self.indices.len());
                        }
                        self.success += 1;
                        return Some(sample);
                    }
                }
            }
            subdivide(&mut self.indices, &self.grid, &self.poisson, self.level);
            if self.indices.is_empty() {
                return None;
            }
            self.range = Range::new(0, self.indices.len());
            self.throws = (self.poisson.a * self.indices.len() as f64).ceil() as usize;
            self.level += 1;
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dim(None);
        let side = 2usize.pow(self.level as u32);
        let spacing = self.grid.cell / side as f64;
        let cell_volume = spacing.powi(dim as i32);
        let sphere_volume = sphere_volume(2. * self.poisson.radius, dim as u64);
        let mut lower = ((self.indices.len() as f64 * cell_volume) /
                         sphere_volume)
                            .floor() as usize;
        if lower > 0 {
            lower -= 1;
        }
        // Calculating upper bound should work because there is just this many places left in the grid and no more can fit into it.
        let upper = self.grid.cells() - self.success;
        (lower, Some(upper))
    }
}

fn subdivide<R, V>(indices: &mut Vec<V>, grid: &Grid<V>, poisson: &PoissonGen<R, V>, level: usize)
    where R: Rng,
          V: VecLike
{
    let choices = &[0., 1.];
    indices.flat_map_inplace(|i| {
        each_combination::<V>(choices)
            .map(move |n| n + i * 2.)
            .filter(|c| !covered(grid, poisson, *c, level + 1))
    });
}

fn covered<R, V>(grid: &Grid<V>, poisson: &PoissonGen<R, V>, index: V, level: usize) -> bool
    where R: Rng,
          V: VecLike
{
    let parent = get_parent(index, level);
    each_combination(&[-2., -1., 0., 1., 2.])
        .filter_map(|t| grid.get(parent + t))
        .filter_map(|t| *t)
        .any(|v| is_cell_covered(&v, grid, poisson, index, level))
}

fn is_cell_covered<R, V>(v: &V,
                         grid: &Grid<V>,
                         poisson: &PoissonGen<R, V>,
                         index: V,
                         level: usize)
                         -> bool
    where R: Rng,
          V: VecLike
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / side as f64;
    let sqradius = (2. * poisson.radius).powi(2);
    each_combination(&[0., 1.])
        .map(|t| (index + t) * spacing)
        .all(|t| sqdist(t, *v, poisson.periodicity) < sqradius)
}

fn is_disk_free<R, V>(grid: &Grid<V>,
                      poisson: &PoissonGen<R, V>,
                      index: V,
                      level: usize,
                      c: V)
                      -> bool
    where R: Rng,
          V: VecLike
{
    let parent = get_parent(index, level);
    let sqradius = (2. * poisson.radius).powi(2);
    // TODO: This does unnessary checks at corners...
    each_combination(&[-2., -1., 0., 1., 2.])
        .filter_map(|t| grid.get(parent + t))
        .filter_map(|t| *t)
        .all(|v| sqdist(v, c, poisson.periodicity) >= sqradius)
}

fn sqdist<V>(v1: V, v2: V, periodicity: bool) -> f64
    where V: VecLike
{
    let diff = v2 - v1;
    if periodicity {
        each_combination(&[-1., 0., 1.])
            .map(|v| (diff + v).sqnorm())
            .fold(f64::MAX, |a, b| a.min(b))
    } else {
        diff.sqnorm()
    }
}

fn choose_random_sample<V, R>(rand: &mut R, grid: &Grid<V>, index: V, level: usize) -> V
    where V: VecLike,
          R: Rng
{
    let dim = V::dim(None);
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / side as f64;
    let mut result = index * spacing;
    for n in 0..dim {
        let place = f64::rand(rand);
        result[n] += place * spacing;
    }
    result
}

#[test]
fn random_point_is_between_right_values_top_lvl() {
    use rand::{SeedableRng, XorShiftRng};
    let mut rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let radius = 0.2;
    let grid = Grid::<na::Vec2<f64>>::new(radius, false);
    for _ in 0..1000 {
        let result = choose_random_sample(&mut rand, &grid, na::Vec2::<f64>::zero(), 0);
        assert!(result.x >= 0.);
        assert!(result.x < grid.cell);
        assert!(result.y >= 0.);
        assert!(result.y < grid.cell);
    }
}

fn encode<V>(v: &V, side: usize, periodicity: bool) -> Option<usize>
    where V: VecLike
{
    let mut index = 0;
    for n in 0..V::dim(None) {
        let mut cur = v[n] as usize;
        if periodicity {
            cur = (v[n] as isize).modulo(side as isize) as usize;
        } else if v[n] < 0. || v[n] >= side as f64 {
            return None;
        }
        index = (index + cur) * side;
    }
    Some(index / side)
}

#[cfg(test)]
fn decode<V>(index: usize, side: usize) -> Option<V>
    where V: VecLike
{
    let dim = V::dim(None);
    if index >= side.pow(dim as u32) {
        return None;
    }
    let mut result = V::zero();
    let mut last = index;
    for n in (0..dim).rev() {
        let cur = last / side;
        let value = (last - cur * side) as f64;
        result[n] = value;
        last = cur;
    }
    Some(result)
}

#[test]
fn encoding_decoding_works() {
    let n = na::Vec2::new(10., 7.);
    assert_eq!(n, decode(encode(&n, 15, false).unwrap(), 15).unwrap());
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = na::Vec2::new(14., 14.);
    assert_eq!(n, decode(encode(&n, 15, false).unwrap(), 15).unwrap());
}

#[test]
fn encoding_outside_of_area_fails() {
    let n = na::Vec2::new(9., 7.);
    assert_eq!(None, encode(&n, 9, false));
    let n = na::Vec2::new(7., 9.);
    assert_eq!(None, encode(&n, 9, false));
}

#[test]
fn decoding_outside_of_area_fails() {
    assert_eq!(None, decode::<na::Vec2<f64>>(100, 10));
}

fn get_parent<V>(mut index: V, level: usize) -> V
    where V: VecLike
{
    let dim = V::dim(None);
    let split = 2usize.pow(level as u32);
    for n in 0..dim {
        index[n] = (index[n] / split as f64).floor();
    }
    index
}

#[test]
fn getting_parent_works() {
    let divides = 4;
    let cells_per_cell = 2usize.pow(divides as u32);
    let testee = na::Vec2::new(1., 2.);
    assert_eq!(testee,
               get_parent((testee * cells_per_cell as f64) + na::Vec2::new(0., 15.),
                          divides));
}
