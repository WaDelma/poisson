//! Module that contains traits that describe poisson-disk distribution generating algorithms.

use {PoissonDisk, VecLike, FloatLike};

use rand::Rng;

use std::fmt::Debug;

pub use self::bridson::Bridson;
pub use self::ebeida::Ebeida;

mod bridson;
mod ebeida;

/// Trait for building algorithms.
pub trait AlgorithmCreator<F, V>: Copy + Debug
    where F: FloatLike,
          V: VecLike<F>,
{
    type Algo: Algorithm<F, V>;

    /// Creates new algorithm.
    fn create(&PoissonDisk<F, V>) -> Self::Algo;
}

/// Trait that describes what poisson-disk distribution generating algorithm needs.
pub trait Algorithm<F, V>
    where F: FloatLike,
          V: VecLike<F>,
{
    /// Advances algorithm based on PoissonDisk and Rng.
    fn next<R>(&mut self, &mut PoissonDisk<F, V>, &mut R) -> Option<V> where R: Rng;

    /// Return lower and upper bound of samples remaining for algorithm to generate based on PoissonDisk.
    fn size_hint(&self, &PoissonDisk<F, V>) -> (usize, Option<usize>);

    /// Restricts the algorithm with arbitary sample.
    fn restrict(&mut self, V);

    /// Checks if sample is valid based on PoissonDisk.
    fn stays_legal(&self, &PoissonDisk<F, V>, V) -> bool;
}
