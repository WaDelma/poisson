use ::{PoissonDisk, VecLike};
use utils::*;
use algo::PoissonAlgorithm;

use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

use sphere::sphere_volume;

static mut COUNTER: usize = 0;

#[derive(Clone)]
pub struct BridsonAlgorithm<V>
    where V: VecLike
{
    grid: Grid<V>,
    active_samples: Vec<V>,
    outside: Vec<V>,
    success: usize,
}

impl<V> PoissonAlgorithm<V> for BridsonAlgorithm<V>
    where V: VecLike,
{
    fn new(poisson: &PoissonDisk<V>) -> Self {
        BridsonAlgorithm {
            grid: Grid::new(poisson.radius, poisson.poisson_type),
            active_samples: vec![],
            outside: vec![],
            success: 0,
        }
    }

    fn next<R>(&mut self, poisson: &mut PoissonDisk<V>, rng: &mut R) -> Option<V>
        where R: Rng
    {
        while !self.active_samples.is_empty() {
            let index = Range::new(0, self.active_samples.len()).ind_sample(rng);
            let cur = self.active_samples[index];
            for _ in 0..30 {
                let sample = cur + random_point_annulus(rng, 2. * poisson.radius, 4. * poisson.radius);
                if sample.iter().all(|&c| 0. <= c && c <= 1.) {
                    let index = sample_to_index(&sample, self.grid.side);
                    if self.insert_if_valid(poisson, index, sample) {
                        return Some(sample);
                    }
                }
            }
            self.active_samples.swap_remove(index);
            if unsafe{::SEED} == 2 {
                // ::debug::visualise_3d(unsafe{COUNTER += 1; COUNTER}, 0, &self.grid, &vec![], &vec![], 0.5 * poisson.radius);
            }
        }
        if self.success == 0 {
            loop {
                let cell = Range::new(0, self.grid.cells()).ind_sample(rng);
                let index = decode(cell, self.grid.side)
                    .expect("Because we are decoding random index within grid this should work.");
                let sample = choose_random_sample(rng, &self.grid, index, 0);
                if self.insert_if_valid(poisson, index, sample) {
                    return Some(sample);
                }
            }
        }
        None
    }

    fn size_hint(&self, poisson: &PoissonDisk<V>) -> (usize, Option<usize>) {
        // Calculating upper bound should work because there is this many places left in the grid and no more can fit into it.
        let upper = self.grid.cells() - self.success;
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and
        // how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dim(None);
        let spacing = self.grid.cell as f64;
        let grid_volume = upper as f64 * spacing.powi(dim as i32);
        let sphere_volume = sphere_volume(2. * poisson.radius, dim as u64);
        let mut lower = (grid_volume / sphere_volume).floor() as usize;
        if lower > 0 {
            lower -= 1;
        }

        (lower, Some(upper))
    }

    fn insert(&mut self, sample: V) {
        self.success += 1;
        let index = sample_to_index(&sample, self.grid.side);
        if let Some(g) = self.grid.get_mut(index) {
            g.push(sample);
        } else {
            self.outside.push(sample);
        }
    }

    fn stays_legal(&self, poisson: &PoissonDisk<V>, sample: V) -> bool {
        let index = sample_to_index(&sample, self.grid.side);
        is_disk_free(&self.grid, poisson.radius, poisson.poisson_type, index, 0, sample) &&
        is_valid(poisson.radius, poisson.poisson_type, &self.outside, sample)
    }
}

impl<V> BridsonAlgorithm<V>
    where V: VecLike
{
    fn insert_if_valid(&mut self, poisson: &mut PoissonDisk<V>, index: V, sample: V) -> bool
    {
        if is_disk_free(&self.grid, poisson.radius, poisson.poisson_type, index, 0, sample) && is_valid(poisson.radius, poisson.poisson_type, &self.outside, sample) {
            self.active_samples.push(sample);
            self.grid.get_mut(index)
                .expect("Because the sample is [0, 1] indexing it should work.")
                .push(sample);
            self.success += 1;
            true
        } else {
            false
        }
    }
}

fn random_point_annulus<V, R>(rand: &mut R, min: f64, max: f64) -> V
    where V: VecLike,
          R: Rng
{
    loop {
        let mut result = V::zero();
        for c in result.iter_mut() {
            *c = ::rand::distributions::normal::StandardNormal::rand(rand).0;
        }
        let result = result.normalize() * f64::rand(rand) * max;
        if result.norm() >= min {
            return result;
        }
    }
}
