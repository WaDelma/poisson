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
//! use poisson::{PoissonDisk, PoissonType, EbeidaAlgorithm};
//!
//! extern crate rand;
//!
//! extern crate nalgebra as na;
//! type Vec2 = na::Vec2<f64>;
//!
//! fn main() {
//!     let poisson =
//!         PoissonDisk::<_, Vec2>::with_radius(0.1, PoissonType::Normal)
//!             .build(rand::weak_rng(), EbeidaAlgorithm::new);
//!     let samples = poisson.generate();
//!     println!("{:?}", samples);
//! }
//! ```
extern crate image;
extern crate modulo;

extern crate sphere;

extern crate rand;
use rand::{Rand, Rng};

extern crate num;
use num::{NumCast, Float, Zero, One};

extern crate nalgebra as na;
use na::{BaseFloat, Dim, Norm, Iterable, IterableMut};

#[macro_use]
extern crate lazy_static;

use std::cmp::PartialEq;
use std::ops::{Sub, Mul, Add, Div};
use std::marker::PhantomData;

pub use algo::{Ebeida, Bridson};
use algo::{AlgorithmCreator, Algorithm};

pub mod algo;
pub mod math;
pub mod utils;
mod debug;

pub static mut SEED: usize = 0;

/// Describes what traits floats have.
pub trait FloatLike:
    BaseFloat +
    Float +
    Rand
{
    /// Casts usize to float.
    fn f(n: usize) -> Self {
        NumCast::from(n).expect("Casting usize to float should always succeed.")
    }
}
impl<T> FloatLike for T
    where T: BaseFloat + Float + Rand
{}

/// Describes what traits the algorithms need to be able to work.
pub trait VecLike<F>:
    IterableMut<F> +
    Iterable<F> +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<F, Output = Self> +
    Div<F, Output = Self> +
    Norm<F> +
    PartialEq +
    Zero +
    One +
    Dim +
    Clone
    where F: FloatLike
{}
impl<T, F> VecLike<F> for T
    where F: FloatLike,
          T: IterableMut<F> +
             Iterable<F> +
             Add<Output = T> +
             Sub<Output = T> +
             Mul<F, Output = T> +
             Div<F, Output = T> +
             Norm<F> +
             PartialEq +
             Zero +
             One +
             Dim +
             Clone
{}

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
pub struct PoissonDisk<F, V>
    where F: FloatLike,
          V: VecLike<F>,
{
    radius: F,
    poisson_type: PoissonType,
    _marker: PhantomData<V>,
}

impl<V, F> PoissonDisk<F, V>
    where F: FloatLike,
          V: VecLike<F>,
{
    /// New PoissonDisk with type of distribution and radius specified.
    /// The radius should be ]0, √2 / 2]
    pub fn with_radius(radius: F, poisson_type: PoissonType) -> Self {
        assert!(F::f(0) < radius);
        assert!(radius <= NumCast::from(2f64.sqrt() / 2.).unwrap());
        PoissonDisk {
            radius: radius,
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New PoissonDisk with type of distribution and relative radius specified.
    /// The relative radius should be ]0, 1]
    pub fn with_relative_radius(relative: F, poisson_type: PoissonType) -> Self {
        assert!(relative >= F::f(0));
        assert!(relative <= F::f(1));
        PoissonDisk {
            radius: relative * NumCast::from(2f64.sqrt() / 2.).unwrap(),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New PoissonDisk with type of distribution, approximate amount of samples and relative radius specified.
    /// The amount of samples should be larger than 0.
    /// The relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    pub fn with_samples(samples: usize, relative: F, poisson_type: PoissonType) -> Self {
        PoissonDisk {
            radius: math::calc_radius::<F, V>(samples, relative, poisson_type),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
        self.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> PoissonType {
        self.poisson_type
    }

    /// Builds PoissonGen with random number generator and algorithm specified.
    pub fn build<R, A>(self, rng: R, _algo: A) -> PoissonGen<F, V, R, A>
        where R: Rng,
              A: AlgorithmCreator<F, V>,
    {
        PoissonGen::new(self, rng)
    }
}

/// Generates poisson-disk distribution for [0, 1]² area.
#[derive(Clone, Debug)]
pub struct PoissonGen<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: AlgorithmCreator<F, V>,
{
    poisson: PoissonDisk<F, V>,
    rng: R,
    _algo: PhantomData<A>,
}

impl<F, V, R, A> PoissonGen<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: AlgorithmCreator<F, V>,
{
    fn new(poisson: PoissonDisk<F, V>, rng: R) -> Self {
        PoissonGen {
            rng: rng,
            poisson: poisson,
            _algo: PhantomData,
        }
    }

    /// Sets the radius of the generator.
    pub fn set_radius(&mut self, radius: F) {
        assert!(F::f(0) < radius);
        assert!(radius <= NumCast::from(2f64.sqrt() / 2.).unwrap());
        self.poisson.radius = radius;
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
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

impl<F, V, R, A> IntoIterator for PoissonGen<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: AlgorithmCreator<F, V>,
{
    type IntoIter = PoissonIter<F, V, R, A::Algo>;
    type Item = V;

    fn into_iter(self) -> Self::IntoIter {
        let algo = A::create(&self.poisson);
        PoissonIter {
            rng: self.rng,
            poisson: self.poisson,
            algo: algo,
        }
    }
}

/// Iterator for generating poisson-disk distribution.
#[derive(Clone)]
pub struct PoissonIter<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: Algorithm<F, V>,
{
    poisson: PoissonDisk<F, V>,
    rng: R,
    algo: A,
}

impl<F,V, R, A> Iterator for PoissonIter<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: Algorithm<F, V>,
{
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        self.algo.next(&mut self.poisson, &mut self.rng)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.algo.size_hint(&self.poisson)
    }
}

impl<F, V, R, A> PoissonIter<F, V, R, A>
    where F: FloatLike,
          V: VecLike<F>,
          R: Rng,
          A: Algorithm<F, V>,
{
    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
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
