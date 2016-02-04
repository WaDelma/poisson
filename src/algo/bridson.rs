use ::{PoissonGen, VecLike};

use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

use sphere::sphere_volume;

use utils::*;

#[derive(Clone)]
pub struct PoissonAlgo<V>
    where V: VecLike
{
    grid: Grid<V>,
    active_samples: Vec<V>,
    outside: Vec<V>,
    success: usize,
}

impl<V> PoissonAlgo<V>
    where V: VecLike
{
    pub fn new<R>(poisson: &mut PoissonGen<R, V>) -> PoissonAlgo<V>
        where R: Rng
    {
        PoissonAlgo {
            grid: Grid::new(poisson.radius, poisson.periodicity),
            active_samples: vec![],
            outside: vec![],
            success: 0,
        }
    }

    pub fn next<R>(&mut self, poisson: &mut PoissonGen<R, V>) -> Option<V>
        where R: Rng
    {
        while !self.active_samples.is_empty() {
            let index = Range::new(0, self.active_samples.len()).ind_sample(&mut poisson.rand);
            let cur = self.active_samples[index];
            for _ in 0..30 {
                let sample = cur + random_point_annulus::<V, R>(&mut poisson.rand, 2. * poisson.radius, 4. * poisson.radius);
                if sample.iter().all(|&v| 0. <= v && v <= 1.) {
                    let index = sample_to_index(&sample, self.grid.side);
                    if self.insert_if_valid(poisson, index, sample) {
                        return Some(sample);
                    }
                }
            }
            self.active_samples.swap_remove(index);
        }
        if self.success == 0 {
            loop {
                let cell = Range::new(0, self.grid.cells()).ind_sample(&mut poisson.rand);
                let index = decode(cell, self.grid.side)
                    .expect("Because we are decoding random index within grid this should work.");
                let sample = choose_random_sample(&mut poisson.rand, &self.grid, index, 0);
                if self.insert_if_valid(poisson, index, sample) {
                    return Some(sample);
                }
            }
        }
        None
    }

    fn insert_if_valid<R>(&mut self, poisson: &mut PoissonGen<R, V>, i: V, sample: V) -> bool
        where R: Rng
    {
        if is_disk_free::<R, V>(&self.grid, poisson.radius, poisson.periodicity, i, 0, sample) && is_valid::<R, V>(poisson.radius, poisson.periodicity, &self.outside, sample) {
            self.active_samples.push(sample);
            self.grid.get_mut(i)
                .expect("Because the sample is [0, 1] indexing it should work.")
                .push(sample);
            self.success += 1;
            true
        } else {
            false
        }
    }

    pub fn size_hint<R>(&self, poisson: &PoissonGen<R, V>) -> (usize, Option<usize>)
        where R: Rng
    {
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

    pub fn insert(&mut self, sample: V) {
        self.success += 1;
        let index = sample_to_index(&sample, self.grid.side);
        if let Some(g) = self.grid.get_mut(index) {
            g.push(sample);
        } else {
            self.outside.push(sample);
        }
    }

    pub fn stays_legal<R>(&self, poisson: &PoissonGen<R, V>, sample: V) -> bool
        where R: Rng
    {
        let index = sample_to_index(&sample, self.grid.side);
        is_disk_free::<R, V>(&self.grid, poisson.radius, poisson.periodicity, index, 0, sample) &&
        is_valid::<R, V>(poisson.radius, poisson.periodicity, &self.outside, sample)
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
