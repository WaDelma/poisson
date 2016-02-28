use ::{PoissonDisk, VecLike};

use rand::Rng;

pub use self::bridson::BridsonAlgorithm;
pub use self::ebeida::EbeidaAlgorithm;

mod bridson;
mod ebeida;

pub trait PoissonAlgorithm<V>
    where V: VecLike,
{
    fn new(poisson: &PoissonDisk<V>) -> Self;

    fn next<R>(&mut self, poisson: &mut PoissonDisk<V>, rng: &mut R) -> Option<V>
        where R: Rng;

    fn size_hint(&self, poisson: &PoissonDisk<V>) -> (usize, Option<usize>);

    fn insert(&mut self, sample: V);

    fn stays_legal(&self, poisson: &PoissonDisk<V>, sample: V) -> bool;
}
