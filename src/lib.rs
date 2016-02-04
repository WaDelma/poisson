//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//! * For each point there is disk of certain radius which doesn't intersect
//! with any other disk of other points
//! * Samples fill the space uniformly
//!
extern crate modulo;

extern crate sphere;

extern crate rand;
use rand::{Rand, Rng};
extern crate num;
use num::{Zero, One};

extern crate nalgebra as na;
use na::{Dim, Norm, Iterable, IterableMut};

#[macro_use]
extern crate lazy_static;

use std::cmp::PartialEq;
use std::ops::{Sub, Mul, Add, Div};
use std::marker::PhantomData;

pub use algo::*;

mod algo;
mod math;
mod utils;
mod debug;

/// Describes what traits the algorithm needs to be able to work.
pub trait VecLike:
    IterableMut<f64> +
    Iterable<f64> +
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
    IterableMut<f64> +
    Iterable<f64> +
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
#[derive(Default, Clone, Debug)]
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

/// Generates Poisson-disk distribution in [0, 1]² area with O(N) time and space complexity relative to the number of samples generated.
/// Based on Ebeida, Mohamed S., et al. "A Simple Algorithm for Maximal Poisson‐Disk Sampling in High Dimensions." Computer Graphics Forum. Vol. 31. No. 2pt4. Blackwell Publishing Ltd, 2012.
#[derive(Clone, Debug)]
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

    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.radius
    }
    /// Returns the periodicity of the generator.
    pub fn periodicity(&self) -> bool {
        self.periodicity
    }

    /// Generates a vector with Poisson-disk distribution in area [0, 1]²
    pub fn generate(&mut self) -> Vec<V> {
        self.into_iter().collect()
    }

    pub fn iter_mut<'a>(&'a mut self) -> PoissonIterMut<'a, R, V> {
        self.into_iter()
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

impl<'a, R, V> IntoIterator for PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    type IntoIter = PoissonIter<R, V>;
    type Item = V;
    fn into_iter(mut self) -> PoissonIter<R, V> {
        let algo = PoissonAlgo::new(&mut self);
        PoissonIter {
            poisson: self,
            algo: algo,
        }
    }
}

#[derive(Clone)]
pub struct PoissonIter<R, V>
    where R: Rng,
          V: VecLike
{
    poisson: PoissonGen<R, V>,
    algo: PoissonAlgo<V>,
}

impl<R, V> PoissonIter<R, V>
    where R: Rng,
          V: VecLike
{
    pub fn insert(&mut self, value: V) {
        self.algo.insert(value);
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.poisson.radius()
    }
    /// Returns the periodicity of the generator.
    pub fn periodicity(&self) -> bool {
        self.poisson.periodicity()
    }

    pub fn stays_legal(&self, value: V) -> bool {
        self.algo.stays_legal(&self.poisson, value)
    }
}

impl<R, V> Iterator for PoissonIter<R, V>
    where R: Rng,
          V: VecLike
{
    type Item = V;
    fn next(&mut self) -> Option<V> {
        self.algo.next(&mut self.poisson)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.algo.size_hint(&self.poisson)
    }
}

impl<R, V> Extend<V> for PoissonIter<R, V>
    where R: Rng,
          V: VecLike
{
    fn extend<T>(&mut self, iterable: T) where T: IntoIterator<Item=V> {
        for t in iterable {
            self.insert(t);
        }
    }
}

impl<'a, R, V> IntoIterator for &'a mut PoissonGen<R, V>
    where R: Rng,
          V: VecLike
{
    type IntoIter = PoissonIterMut<'a, R, V>;
    type Item = V;
    fn into_iter(self) -> PoissonIterMut<'a, R, V> {
        let algo = PoissonAlgo::new(self);
        PoissonIterMut {
            poisson: self,
            algo: algo,
        }
    }
}

pub struct PoissonIterMut<'a, R, V>
    where R: Rng + 'a,
          V: VecLike + 'a
{
    poisson: &'a mut PoissonGen<R, V>,
    algo: PoissonAlgo<V>,
}

impl<'a, R, V> PoissonIterMut<'a, R, V>
    where R: Rng,
          V: VecLike
{
    pub fn insert(&mut self, value: V) {
        self.algo.insert(value);
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.poisson.radius()
    }
    /// Returns the periodicity of the generator.
    pub fn periodicity(&self) -> bool {
        self.poisson.periodicity()
    }

    pub fn stays_legal(&self, value: V) -> bool {
        self.algo.stays_legal(&self.poisson, value)
    }
}

impl<'a, R, V> Iterator for PoissonIterMut<'a, R, V>
    where R: Rng,
          V: VecLike
{
    type Item = V;
    fn next(&mut self) -> Option<V> {
        self.algo.next(self.poisson)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.algo.size_hint(&self.poisson)
    }
}

impl<'a, R, V> Extend<V> for PoissonIterMut<'a, R, V>
    where R: Rng,
          V: VecLike
{
    fn extend<T>(&mut self, iterable: T) where T: IntoIterator<Item=V> {
        for t in iterable {
            self.insert(t);
        }
    }
}
