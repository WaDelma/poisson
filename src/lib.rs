//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//! * For each point there is disk of certain radius which doesn't intersect
//! with any other disk of other points
//! * Samples fill the space uniformly
//!
//! # Examples
//!
//! ```
//! extern crate poisson;
//! use poisson::PoissonDisk;
//!
//! extern crate rand;
//!
//! extern crate nalgebra as na;
//! type Vec2 = na::Vec2<f64>;
//!
//! fn main() {
//!     let poisson =
//!         PoissonDisk::<Vec2>::with_radius(0.1, PoissonType::Normal)
//!             .build(rand::weak_rng());
//!     let samples = poisson.generate();
//!     println!("{:?}", samples);
//! }
//! ```
extern crate image;
extern crate modulo;

extern crate sphere;

extern crate rand;
use rand::Rng;

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

pub static mut SEED: usize = 0;

/// Describes what traits the algorithms need to be able to work.
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
    Clone {}
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
    Clone {}

/// Enum for determining the type of poisson-disk distribution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoissonType {
    Normal,
    Perioditic,
}

impl Default for PoissonType {
    fn default() -> PoissonType {
        PoissonType::Normal
    }
}


/// Builder for PoissonGen.
#[derive(Default, Clone, Debug, PartialEq)]
pub struct PoissonDisk<V>
    where V: VecLike,
{
    radius: f64,
    poisson_type: PoissonType,
    _marker: PhantomData<V>,
}

impl<V> PoissonDisk<V>
    where V: VecLike,
{
    /// New PoissonDisk with type of distribution and radius specified.
    /// The radius should be ]0, √2 / 2]
    pub fn with_radius(radius: f64, poisson_type: PoissonType) -> Self {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        PoissonDisk {
            radius: radius,
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New PoissonDisk with type of distribution and relative radius specified.
    /// The relative radius should be ]0, 1]
    pub fn with_relative_radius(relative: f64, poisson_type: PoissonType) -> Self {
        assert!(relative >= 0.);
        assert!(relative <= 1.);
        PoissonDisk {
            radius: relative * (2f64.sqrt() / 2.),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New PoissonDisk with type of distribution, approximate amount of samples and relative radius specified.
    /// The amount of samples should be larger than 0.
    /// The relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    pub fn with_samples(samples: u32, relative: f64, poisson_type: PoissonType) -> Self {
        assert!(PoissonType::Perioditic == poisson_type || V::dim(None) < 5);
        assert!(samples > 0);
        assert!(relative >= 0.);
        assert!(relative <= 1.);
        PoissonDisk {
            radius: math::calc_radius::<V>(samples, relative, poisson_type),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> PoissonType {
        self.poisson_type
    }

    /// Builds PoissonGen with random number generator and default algorithm.
    pub fn build<R>(self, rng: R) -> PoissonGen<V, R, algo::EbeidaAlgorithm<V>>
        where R: Rng,
    {
        PoissonGen::new(self, rng)
    }

    /// Builds PoissonGen with random number generator and algorithm specified.
    pub fn build_with_algo<R, A>(self, rng: R) -> PoissonGen<V, R, A>
        where R: Rng,
              A: PoissonAlgorithm<V>,
    {
        PoissonGen::new(self, rng)
    }
}

/// Generates poisson-disk distribution for [0, 1]² area.
#[derive(Clone, Debug)]
pub struct PoissonGen<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    poisson: PoissonDisk<V>,
    rng: R,
    _algo: PhantomData<A>,
}

impl<V, R, A> PoissonGen<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    fn new(poisson: PoissonDisk<V>, rng: R) -> Self {
        PoissonGen {
            rng: rng,
            poisson: poisson,
            _algo: PhantomData,
        }
    }

    /// Sets the radius of the generator.
    pub fn set_radius(&mut self, radius: f64) {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        self.poisson.radius = radius;
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.poisson.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> PoissonType {
        self.poisson.poisson_type
    }

    /// Generates Poisson-disk distribution.
    pub fn generate(self) -> Vec<V> {
        self.into_iter().collect()
    }
}

impl<V, R, A> IntoIterator for PoissonGen<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    type IntoIter = PoissonIter<V, R, A>;
    type Item = V;

    fn into_iter(self) -> Self::IntoIter {
        let algo = A::new(&self.poisson);
        PoissonIter {
            rng: self.rng,
            poisson: self.poisson,
            algo: algo,
        }
    }
}

/// Iterator for generating poisson-disk distribution.
#[derive(Clone)]
pub struct PoissonIter<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    poisson: PoissonDisk<V>,
    rng: R,
    algo: A,
}

impl<V, R, A> Iterator for PoissonIter<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    type Item = V;
    fn next(&mut self) -> Option<Self::Item> {
        self.algo.next(&mut self.poisson, &mut self.rng)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.algo.size_hint(&self.poisson)
    }
}

impl<V, R, A> PoissonIter<V, R, A>
    where V: VecLike,
          R: Rng,
          A: PoissonAlgorithm<V>,
{
    /// Returns the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.poisson.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> PoissonType {
        self.poisson.poisson_type
    }

    /// Restricts the poisson algorithm with arbitary sample.
    pub fn restrict(&mut self, value: V) {
        self.algo.restrict(value);
    }

    /// Checks legality of sample for currrent distribution.
    pub fn stays_legal(&self, value: V) -> bool {
        self.algo.stays_legal(&self.poisson, value)
    }
}
