use ::{PoissonDisk, VecLike};

use rand::Rng;

pub use self::bridson::BridsonAlgorithm;
pub use self::ebeida::EbeidaAlgorithm;

mod bridson;
mod ebeida;

/// Trait that describes what poisson-disk distribution generating algorithm needs.
pub trait PoissonAlgorithm<V>
    where V: VecLike,
{
    /// Creates new algorithm based on PoissonDisk.
    fn new(poisson: &PoissonDisk<V>) -> Self;

    /// Advances algorithm based on PoissonDisk and Rng.
    fn next<R>(&mut self, poisson: &mut PoissonDisk<V>, rng: &mut R) -> Option<V>
        where R: Rng;

    /// Return lower and upper bound of samples remaining for algorithm to generate based on PoissonDisk.
    fn size_hint(&self, poisson: &PoissonDisk<V>) -> (usize, Option<usize>);

    /// Restricts the algorithm with arbitary sample.
    fn restrict(&mut self, sample: V);

    /// Checks if sample is valid based on PoissonDisk.
    fn stays_legal(&self, poisson: &PoissonDisk<V>, sample: V) -> bool;
}
