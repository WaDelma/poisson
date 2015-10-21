//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//!    * For each point there is disk of certain radius which doesn't intersect with any disk of any other point
//!    * Nodes fill the space uniformly
//!

extern crate image;

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

use std::cmp::PartialEq;
use std::ops::{Sub, Mul, Add, Div, IndexMut};
use std::marker::PhantomData;

use utils::{each_combination, Inplace};

mod math;
#[cfg(test)]
mod test;
mod debug;
mod utils;

/// Generates poisson-disk distribution in [0, 1]² area with O(N log N) time and space complexity relative to the number of samples generated.
/// Based on Gamito, Manuel N., and Steve C. Maddock. "Accurate multidimensional Poisson-disk sampling." ACM Transactions on Graphics (TOG) 29.1 (2009): 8.
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
///     let mut vecs = vec![];
///     poisson.generate(&mut vecs);
///     println!("{:?}", vecs);
/// }
/// ```
pub struct PoissonDisk<R: Rng> {
    rand: R,
    periodicity: bool,
}

impl<R: Rng> PoissonDisk<R> {
    /// Creates new poisson-disk generator builder with random generator specified.
    pub fn new(rand: R) -> Self {
        PoissonDisk{rand: rand, periodicity: false}
    }

    /// Sets the generator to generate perioditic poisson-disk distributions.
    pub fn perioditic(mut self) -> Self {
        self.periodicity = true;
        self
    }

    /// Builds the generator with relative radius specified.
    /// Radius should be ]0, 1]
    pub fn build_relative_radius<V: VecLike>(self, radius: f64) -> PoissonGen<R, V> {
        assert!(0. < radius);
        assert!(radius <= 1.);
        PoissonGen {
            dim: PhantomData,
            radius: radius * (2f64.sqrt() / 2.),
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }

    /// Builds the generator with radius specified.
    /// Radius should be ]0, √2 / 2]
    pub fn build_radius<V: VecLike>(self, radius: f64) -> PoissonGen<R, V> {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        PoissonGen {
            dim: PhantomData,
            radius: radius,
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }

    /// Builds the generator with radius calculated so that approximately specified number of samples are generated.
    /// Amount of samples should be larger than 0.
    /// Relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    pub fn build_samples<V: VecLike>(self, samples: u32, relative_radius: f64) -> PoissonGen<R, V> {
        assert!(self.periodicity || V::dim(None) < 5);
        assert!(samples > 0);
        assert!(relative_radius >= 0.);
        assert!(relative_radius <= 1.);
        PoissonGen {
            dim: PhantomData,
            radius: math::calc_radius::<V>(samples, relative_radius, self.periodicity),
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }
}

pub struct PoissonGen<R: Rng, V: VecLike> {
    dim: PhantomData<V>,
    rand: R,
    radius: f64,
    periodicity: bool,
}

impl<R: Rng, V: VecLike> PoissonGen<R, V> {
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

    /// Populates given vector with poisson-disk distribution [0, 1]²
    /// Resulting samples will be a poisson-disk distribution iff given samples were already valid poisson-disk distribution.
    /// Resulting samples will be a maximal poisson-disk distribution [0, 1]² iff given samples have same radius and are already valid poisson-disk distribution.
    pub fn generate(&mut self, points: &mut Vec<Sample<V>>) {
        let dim = V::dim(None);
        let top_lvl_cell = self.radius / (dim as f64).sqrt();
        let top_lvl_side = (1. / top_lvl_cell) as usize;
        let mut grid = vec![None; top_lvl_side.pow(dim as u32)];
        let capacity = grid.len() * dim;
        let mut indices = Vec::with_capacity(capacity);
        indices.extend((0..grid.len()));
        //let child_amount = 2usize.pow(V::dim(None) as u32);
        let mut level = 0;
        let max_level = 63;
        let a = 0.3;
        while !indices.is_empty() && level < max_level {
            // if level > 6 {
            //     panic!();
            // }
            println!("{}/63, {}/{}, {}/{}", level, indices.len(), (top_lvl_side * 2usize.pow(level as u32)).pow(dim as u32), grid.iter().filter(|n| n.is_some()).count(), grid.len());
            let mut range = Range::new(0, indices.len());
            let throws = (a * indices.len() as f64) as usize;
            for _ in 0..throws {
                let index = range.ind_sample(&mut self.rand);
                let cur = indices[index];
                assert!(get_parent::<V>(cur, level, top_lvl_side).is_some());
                let parent_v = get_parent::<V>(cur, level, top_lvl_side).unwrap();
                let parent = encode(&parent_v, top_lvl_side).unwrap();
                if grid[parent].is_some() {
                    indices.swap_remove(index);
                    if indices.is_empty() {
                        break;
                    }
                    range = Range::new(0, indices.len());
                } else {
                    let c = choose_random_point(&mut self.rand, cur, level, top_lvl_side, top_lvl_cell);
                    if self.is_disk_free(&grid, cur, level, c, top_lvl_side) {
                        grid[parent] = Some(c);
                        indices.swap_remove(index);
                        range = Range::new(0, indices.len());
                    }
                }
            }

            debug::visualise(level, &grid, top_lvl_side, &indices, top_lvl_cell, self.radius);

            let cells_per_cell = 2usize.pow(level as u32);
            let side = cells_per_cell * top_lvl_side;
            level += 1;
            let next_cells_per_cell = 2usize.pow(level as u32);
            let next_side = next_cells_per_cell * top_lvl_side;
            let choices = &[0., 1.];
            indices.flat_map_inplace(|i| {
                assert!(decode::<V>(i, side).is_some());
                let ind = decode::<V>(i, side).unwrap();
                each_combination::<V>(choices)
                    .map(move |n| encode(&(n + ind * 2.), next_side).unwrap())
                    .filter(|&c| !self.covered(&grid, c, level, top_lvl_side, top_lvl_cell))
            });
            // If this assert fails then a is too small or subdivide code is broken
            // assert_eq!(capacity, indices.capacity());
        }
        points.extend(grid.into_iter().filter_map(|v| v).map(|v| Sample::new(v, self.radius)));
    }
}

impl <R: Rng, V: VecLike> PoissonGen<R, V> {

    fn is_disk_free(&self, grid: &Vec<Option<V>>, index: usize, level: usize, c: V, top_lvl_side: usize) -> bool {
        assert!(get_parent::<V>(index, level, top_lvl_side).is_some());
        let parent = get_parent::<V>(index, level, top_lvl_side).unwrap();
        let sqradius = self.radius.powi(2);
        //TODO: Does unnessary checking...
        for t in each_combination(&[-2., -1., 0., 1., 2.]) {
            if let Some(i) = encode(&(parent + t), top_lvl_side) {
                if let Some(s) = grid[i] {
                    if (s - c).sqnorm() < sqradius {
                        return false;
                    }
                }
            }
        }
        // for cc in 0..grid.len() {
        //     if let Some(s) = grid[cc] {
        //         if (s - c).sqnorm() < sqradius {
        //             return false;
        //         }
        //     }
        // }
        true
    }

    fn covered(&self, grid: &Vec<Option<V>>, index: usize, level: usize, top_lvl_side: usize, top_lvl_cell: f64) -> bool {
        assert!(get_parent::<V>(index, level, top_lvl_side).is_some());
        let parent = get_parent::<V>(index, level, top_lvl_side).unwrap();
        for t in each_combination(&[-1., 0., 1.]) {
            if let Some(i) = encode(&(parent + t), top_lvl_side) {
                if let Some(s) = grid[i] {
                    if self.is_cell_covered(&s, index, top_lvl_cell, top_lvl_side, level) {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn is_cell_covered(&self, v: &V, index: usize, top_lvl_cell: f64, top_lvl_side: usize, level: usize) -> bool {
        let side = 2usize.pow(level as u32);
        let spacing = top_lvl_cell / side as f64;
        let sqradius = self.radius.powi(2);
        assert!(decode::<V>(index, side * top_lvl_side).is_some());
        let base = decode::<V>(index, side * top_lvl_side).unwrap();
        for t in each_combination(&[-1., 1.]) {
            let vec = (base + t) * spacing;
            if (vec - *v).sqnorm() >= sqradius {
                return false;
            }
        }
        true
    }

}

fn choose_random_point<V: VecLike, R: Rng>(rand: &mut R, index: usize, level: usize, top_lvl_side: usize, top_lvl_cell: f64) -> V {
    let dim = V::dim(None);
    let side = 2usize.pow(level as u32);
    let spacing = top_lvl_cell / side as f64;
    assert!(decode::<V>(index, side * top_lvl_side).is_some());
    let mut result = decode::<V>(index, side * top_lvl_side).unwrap() * spacing;
    for n in 0..dim {
        let place = f64::rand(rand);
        result[n] += place * spacing;//mul_add
    }
    result
}

#[test]
fn random_point_is_between_right_values_top_lvl() {
    use rand::{SeedableRng, XorShiftRng};
    let mut rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let dim  = 2;
    let radius = 0.2;
    let top_lvl_cell = radius / (dim as f64).sqrt();
    let top_lvl_side = (1. / top_lvl_cell) as usize;
    for _ in 0..1000 {
        let result = choose_random_point::<na::Vec2<f64>, _>(&mut rand, 0, 0, top_lvl_side, top_lvl_cell);
        assert!(result.x >= 0.);
        assert!(result.x < top_lvl_cell);
        assert!(result.y >= 0.);
        assert!(result.y < top_lvl_cell);
    }
}

fn encode<V: VecLike>(v: &V, side: usize) -> Option<usize> {
    let mut index = 0;
    for n in 0..V::dim(None) {
        if v[n] < 0. || v[n] >= side as f64 {
            return None;
        }
        index = (index + v[n] as usize) * side;
    }
    Some(index / side as usize)
}

fn decode<V: VecLike>(index: usize, side: usize) -> Option<V> {
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
    assert_eq!(n, decode(encode(&n, 15).unwrap(), 15).unwrap());
}

#[test]
fn encoding_decoding_at_edge_works() {
    let n = na::Vec2::new(14., 14.);
    assert_eq!(n, decode(encode(&n, 15).unwrap(), 15).unwrap());
}

#[test]
fn encoding_outside_of_area_fails() {
    let n = na::Vec2::new(9., 7.);
    assert_eq!(None, encode(&n, 9));
    let n = na::Vec2::new(7., 9.);
    assert_eq!(None, encode(&n, 9));
}

#[test]
fn decoding_outside_of_area_fails() {
    assert_eq!(None, decode::<na::Vec2<f64>>(100, 10));
}

fn get_parent<V: VecLike>(index: usize, level: usize, top_lvl_side: usize) -> Option<V> {
    let dim = V::dim(None);
    let split = 2usize.pow(level as u32);
    //assert!(decode::<V>(index, split * top_lvl_side).is_some());
    decode::<V>(index, split * top_lvl_side)
        .and_then(|mut r| {
            for n in 0..dim {
                if r[n] >= top_lvl_side as f64 {
                    //println!("gp: {} {}", r[n], n);
                    //return None;
                }
                r[n] = (r[n] / split as f64).floor();
            }
            Some(r)
        })
}

#[test]
fn getting_parent_works() {
    let cells_per_side = 3;
    let divides = 4;
    let cells_per_cell = 2usize.pow(divides as u32);
    let cells_per_side_divided = cells_per_side * cells_per_cell;
    let testee = na::Vec2::new(1., 2.);
    let index = encode(&((testee * cells_per_cell as f64) + na::Vec2::new(0., 15.)), cells_per_side_divided).unwrap();
    assert_eq!(Some(testee), get_parent(index, divides, cells_per_side));
}

#[test]
fn getting_parent_outside_of_area_fails() {
    let cells_per_side = 3;
    let divides = 4;
    let cells_per_cell = 2usize.pow(divides as u32);
    let cells_per_side_divided = cells_per_side * cells_per_cell;
    let testee = na::Vec2::new(1., 3.);
    let mut index = 0;
    for n in 0..2 {
        index = (index + testee[n] as usize) * cells_per_side_divided;
    }
    index  = index / cells_per_side_divided as usize;
    assert_eq!(None::<na::Vec2<f64>>, get_parent(index, divides, cells_per_side));
}

/// Describes position of sample and radius of disk around it.
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Sample<T: VecLike> {
    pub pos: T,
    radius: f64,
}

impl<T: VecLike> Sample<T> {
    pub fn new(pos: T, radius: f64) -> Self {
        Sample{pos: pos, radius: radius}
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }
}
