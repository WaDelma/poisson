//! # Poisson-disk distribution generation
//!
//! Generates distribution of points in [0, 1)<sup>d</sup> where:
//!
//! * For each point there is disk of certain radius which doesn't intersect
//! with any other disk of other points
//! * Samples fill the space uniformly
//!
//! Due it's blue noise properties poisson-disk distribution
//! can be used for object placement in procedural texture/world generation,
//! as source distribution for digital stipling,
//! as distribution for sampling in rendering or for (re)meshing.
//!
//! # Examples
//!
//! Generate non-tiling poisson-disk distribution in [0, 1)<sup>2</sup> with disk radius 0.1
//! using slower but more accurate algorithm.
//!
//! ````rust
//! extern crate poisson;
//! extern crate rand;
//! extern crate nalgebra as na;
//!
//! use rand::SeedableRng;
//! use rand::rngs::SmallRng;
//!
//! use poisson::{Builder, Type, algorithm};
//!
//! fn main() {
//!     let poisson =
//!         Builder::<_, na::Vector2<f64>>::with_radius(0.1, Type::Normal)
//!             .build(SmallRng::from_entropy(), algorithm::Ebeida);
//!     let samples = poisson.generate();
//!     println!("{:?}", samples);
//! }
//! ````
//!
//! Generate tiling poisson-disk distribution in [0, 1)<sup>3</sup> with approximately 100 samples
//! and relative disk radius 0.9 using faster but less accurate algorithm.
//!
//! ````rust
//! # extern crate nalgebra as na;
//! # use poisson::{Builder, Type, algorithm};
//! # use rand::SeedableRng;
//! # use rand::rngs::SmallRng;
//!
//! fn main() {
//!     let poisson =
//!         Builder::<_, na::Vector3<f32>>::with_samples(100, 0.9, Type::Perioditic)
//!             .build(SmallRng::from_entropy(), algorithm::Bridson);
//!     for sample in poisson {
//!         println!("{:?}", sample)
//!     }
//! }
//! ````

use rand::Rng;

use num_traits::Float as NumFloat;
use num_traits::{NumCast, Zero};

use alga::general::RealField;
use alga::linear::{FiniteDimVectorSpace, NormedSpace};

#[macro_use]
extern crate lazy_static;

use std::marker::PhantomData;
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

use crate::algorithm::{Algorithm, Creator};
use crate::utils::math::calc_radius;

pub mod algorithm;
mod utils;

/// Describes what floats are.
pub trait Float: NumFloat + RealField + AddAssign + SubAssign + MulAssign + DivAssign {
    /// Casts usize to float.
    fn cast(n: usize) -> Self {
        NumCast::from(n).expect("Casting usize to float should always succeed.")
    }
}
impl<T> Float for T where T: NumFloat + RealField + AddAssign + SubAssign + MulAssign + DivAssign
{}

/// Describes what vectors are.
pub trait Vector<F>:
    Zero
    + FiniteDimVectorSpace
    + NormedSpace<RealField = F, ComplexField = F>
where
    F: Float,
{
}
impl<T, F> Vector<F> for T
where
    F: Float,
    T: Zero
        + FiniteDimVectorSpace
        + NormedSpace<RealField = F, ComplexField = F>,
{
}

/// Enum for determining the type of poisson-disk distribution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Type {
    /// Acts like there is void all around the space placing no restrictions to sides.
    Normal,
    /// Makes the space to wrap around on edges allowing tiling of the generated poisson-disk distribution.
    Perioditic,
}

impl Default for Type {
    fn default() -> Type {
        Type::Normal
    }
}

/// Builder for the generator.
#[derive(Default, Clone, Debug, PartialEq)]
pub struct Builder<F, V>
where
    F: Float,
    V: Vector<F>,
{
    radius: F,
    poisson_type: Type,
    _marker: PhantomData<V>,
}

impl<V, F> Builder<F, V>
where
    F: Float,
    V: Vector<F>,
{
    /// New Builder with type of distribution and radius specified.
    /// The radius should be ]0, √2 / 2]
    pub fn with_radius(radius: F, poisson_type: Type) -> Self {
        assert!(F::cast(0) < radius);
        assert!(
            radius
                <= NumCast::from(2f64.sqrt() / 2.).expect("Casting constant should always work.")
        );
        Builder {
            radius: radius,
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New Builder with type of distribution and relative radius specified.
    /// The relative radius should be ]0, 1]
    pub fn with_relative_radius(relative: F, poisson_type: Type) -> Self {
        assert!(relative >= F::cast(0));
        assert!(relative <= F::cast(1));
        Builder {
            radius: relative
                * NumCast::from(2f64.sqrt() / 2.).expect("Casting constant should always work."),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// New Builder with type of distribution, approximate amount of samples and relative radius specified.
    /// The amount of samples should be larger than 0.
    /// The relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    /// For perioditic this is supported up to 8 dimensions.
    pub fn with_samples(samples: usize, relative: F, poisson_type: Type) -> Self {
        Builder {
            radius: calc_radius::<F, V>(samples, relative, poisson_type),
            poisson_type: poisson_type,
            _marker: PhantomData,
        }
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
        self.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> Type {
        self.poisson_type
    }

    /// Builds generator with random number generator and algorithm specified.
    pub fn build<R, A>(self, rng: R, _algo: A) -> Generator<F, V, R, A>
    where
        R: Rng,
        A: Creator<F, V>,
    {
        Generator::new(self, rng)
    }
}

/// Generates poisson-disk distribution in [0, 1]<sup>d</sup> area.
#[derive(Clone, Debug)]
pub struct Generator<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    A: Creator<F, V>,
{
    poisson: Builder<F, V>,
    rng: R,
    _algo: PhantomData<A>,
}

impl<F, V, R, A> Generator<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    A: Creator<F, V>,
{
    fn new(poisson: Builder<F, V>, rng: R) -> Self {
        Generator {
            rng: rng,
            poisson: poisson,
            _algo: PhantomData,
        }
    }

    /// Sets the radius of the generator.
    pub fn set_radius(&mut self, radius: F) {
        assert!(F::cast(0) < radius);
        assert!(
            radius
                <= NumCast::from(2f64.sqrt() / 2.).expect("Casting constant should always work.")
        );
        self.poisson.radius = radius;
    }

    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
        self.poisson.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> Type {
        self.poisson.poisson_type
    }
}

impl<F, V, R, A> Generator<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
    R: Rng + Clone,
    A: Creator<F, V>,
{
    /// Generates Poisson-disk distribution.
    pub fn generate(&self) -> Vec<V> {
        self.clone().into_iter().collect()
    }
}

impl<F, V, R, A> IntoIterator for Generator<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    A: Creator<F, V>,
{
    type IntoIter = PoissonIter<F, V, R, A::Algo>;
    type Item = V;

    fn into_iter(self) -> Self::IntoIter {
        PoissonIter {
            rng: self.rng,
            algo: A::create(&self.poisson),
            poisson: self.poisson,
        }
    }
}

/// Iterator for generating poisson-disk distribution.
#[derive(Clone)]
pub struct PoissonIter<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    A: Algorithm<F, V>,
{
    poisson: Builder<F, V>,
    rng: R,
    algo: A,
}

impl<F, V, R, A> Iterator for PoissonIter<F, V, R, A>
where
    F: Float,
    V: Vector<F>,
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
where
    F: Float,
    V: Vector<F>,
    R: Rng,
    A: Algorithm<F, V>,
{
    /// Returns the radius of the generator.
    pub fn radius(&self) -> F {
        self.poisson.radius
    }

    /// Returns the type of the generator.
    pub fn poisson_type(&self) -> Type {
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
