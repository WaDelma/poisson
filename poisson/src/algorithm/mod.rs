//! Module that contains traits that describe poisson-disk distribution generating algorithms.

use crate::{Builder, Float, Vector};

use rand::Rng;

use std::fmt::Debug;

pub use self::bridson::Bridson;
pub use self::ebeida::Ebeida;

mod bridson;
mod ebeida;

/// Constructs new instance of the algorithm.
pub trait Creator<F, V>: Copy + Debug
where
    F: Float,
    V: Vector<F>,
{
    /// Algorithm instance associated with the trait
    type Algo: Algorithm<F, V>;

    /// Creates new and empty algorithm instance.
    fn create(_: &Builder<F, V>) -> Self::Algo;
}

/// Trait that describes poisson-disk distribution generating algorithm.
pub trait Algorithm<F, V>
where
    F: Float,
    V: Vector<F>,
{
    /// Generates new sample advancing the algorithm.
    fn next<R>(&mut self, _: &mut Builder<F, V>, _: &mut R) -> Option<V>
    where
        R: Rng;

    /// Returns lower and upper bound of the amount of samples remaining for the algorithm to generate.
    fn size_hint(&self, _: &Builder<F, V>) -> (usize, Option<usize>);

    /// Restricts the algorithm with an arbitary sample.
    fn restrict(&mut self, _: V);

    /// Checks if a sample is valid for the poisson-disk distribution generated thus far by the algorithm.
    fn stays_legal(&self, _: &Builder<F, V>, _: V) -> bool;
}
