use {Builder, Vector, Float};
use algorithm::{Creator, Algorithm};
use utils::*;

use num::NumCast;

use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;
use rand::distributions::normal::StandardNormal;

use sphere::sphere_volume;

/// Generates non-uniform non-maximal Poisson-disk distribution with O(n) time and O(n) space complexity relative to the number of samples generated.
/// Based on Bridson, Robert. "Fast Poisson disk sampling in arbitrary dimensions." SIGGRAPH Sketches. 2007.
#[derive(Debug, Clone, Copy)]
pub struct Bridson;

impl<F, V> Creator<F, V> for Bridson
    where F: Float,
          V: Vector<F>
{
    type Algo = Algo<F, V>;

    fn create(poisson: &Builder<F, V>) -> Self::Algo {
        Algo {
            grid: Grid::new(poisson.radius, poisson.poisson_type),
            active_samples: vec![],
            outside: vec![],
            success: 0,
        }
    }
}

pub struct Algo<F, V>
    where F: Float,
          V: Vector<F>
{
    grid: Grid<F, V>,
    active_samples: Vec<V>,
    outside: Vec<V>,
    success: usize,
}

impl<F, V> Algorithm<F, V> for Algo<F, V>
    where F: Float,
          V: Vector<F>
{
    fn next<R>(&mut self, poisson: &mut Builder<F, V>, rng: &mut R) -> Option<V>
        where R: Rng
    {
        while !self.active_samples.is_empty() {
            let index = Range::new(0, self.active_samples.len()).ind_sample(rng);
            let cur = self.active_samples[index].clone();
            for _ in 0..30 {
                let min = F::cast(2) * poisson.radius;
                let max = F::cast(4) * poisson.radius;
                let sample = cur.clone() + random_point_annulus(rng, min, max);
                if sample.iter().all(|&c| F::cast(0) <= c && c <= F::cast(1)) {
                    let index = sample_to_index(&sample, self.grid.side());
                    if self.insert_if_valid(poisson, index, sample.clone()) {
                        return Some(sample);
                    }
                }
            }
            self.active_samples.swap_remove(index);
        }
        while self.success == 0 {
            let cell = Range::new(0, self.grid.cells()).ind_sample(rng);
            let index: V = decode(cell, self.grid.side())
                               .expect("Because we are decoding random index within grid \
                                        this should work.");
            let sample = choose_random_sample(rng, &self.grid, index.clone(), 0);
            if self.insert_if_valid(poisson, index, sample.clone()) {
                return Some(sample);
            }
        }
        None
    }

    fn size_hint(&self, poisson: &Builder<F, V>) -> (usize, Option<usize>) {
        // Calculating upper bound should work because there is this many places left in the grid and no more can fit into it.
        let upper = if self.grid.cells() > self.success {
            self.grid.cells() - self.success
        } else {
            0
        };
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and
        // how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dim(None);
        let spacing = self.grid.cell();
        let grid_volume = F::cast(upper) * spacing.powi(dim as i32);
        let sphere_volume = sphere_volume(F::cast(2) * poisson.radius, dim as u64);
        let lower: F = grid_volume / sphere_volume;
        let mut lower = lower.floor()
                             .to_usize()
                             .expect("Grids volume divided by spheres volume should be always \
                                      castable to usize.");
        if lower > 0 {
            lower -= 1;
        }
        (lower, Some(upper))
    }

    fn restrict(&mut self, sample: V) {
        self.success += 1;
        let index = sample_to_index(&sample, self.grid.side());
        if let Some(g) = self.grid.get_mut(index) {
            g.push(sample);
        } else {
            self.outside.push(sample);
        }
    }

    fn stays_legal(&self, poisson: &Builder<F, V>, sample: V) -> bool {
        let index = sample_to_index(&sample, self.grid.side());
        is_disk_free(&self.grid, poisson, index, 0, sample.clone(), &self.outside)
    }
}

impl<F, V> Algo<F, V>
    where F: Float,
          V: Vector<F>
{
    fn insert_if_valid(&mut self, poisson: &mut Builder<F, V>, index: V, sample: V) -> bool {
        if is_disk_free(&self.grid,
                        poisson,
                        index.clone(),
                        0,
                        sample.clone(),
                        &self.outside) {
            self.active_samples.push(sample.clone());
            self.grid
                .get_mut(index)
                .expect("Because the sample is [0, 1] indexing it should work.")
                .push(sample);
            self.success += 1;
            true
        } else {
            false
        }
    }
}

fn random_point_annulus<F, V, R>(rand: &mut R, min: F, max: F) -> V
    where F: Float,
          V: Vector<F>,
          R: Rng
{
    loop {
        let mut result = V::zero();
        for c in result.iter_mut() {
            *c = NumCast::from(StandardNormal::rand(rand).0)
                     .expect("The f64 produced by StandardNormal should be always castable to \
                              float.");
        }
        let result = result.normalize() * F::rand(rand) * max;
        if result.norm() >= min {
            return result;
        }
    }
}
