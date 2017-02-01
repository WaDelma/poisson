//! Module that contains traits that describe poisson-disk distribution generating algorithms.

use {PoissonConfiguration, Vector, Float};

use rand::Rng;

use std::fmt::Debug;

pub use self::bridson::Bridson;
pub use self::ebeida::Ebeida;

mod bridson;
mod ebeida;

/// Trait for building algorithms.
pub trait Creator<F, V>: Copy + Debug
    where F: Float,
          V: Vector<F>,
{
    type Algo: Algorithm<F, V>;

    /// Creates new algorithm.
    fn create(&PoissonConfiguration<F, V>) -> Self::Algo;
}

/// Trait that describes what poisson-disk distribution generating algorithm needs.
pub trait Algorithm<F, V>
    where F: Float,
          V: Vector<F>,
{
    /// Advances algorithm based on PoissonConfiguration and Rng.
    fn next<R>(&mut self, &mut PoissonConfiguration<F, V>, &mut R) -> Option<V> where R: Rng;

    /// Return lower and upper bound of samples remaining for algorithm to generate based on PoissonConfiguration.
    fn size_hint(&self, &PoissonConfiguration<F, V>) -> (usize, Option<usize>);

    /// Restricts the algorithm with arbitary sample.
    fn restrict(&mut self, V);

    /// Checks if sample is valid based on PoissonConfiguration.
    fn stays_legal(&self, &PoissonConfiguration<F, V>, V) -> bool;
}
