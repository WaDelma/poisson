use ::{PoissonGen, VecLike};

use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

use sphere::sphere_volume;

use std::f64;

use utils::*;

#[derive(Clone)]
pub struct PoissonAlgo<V>
    where V: VecLike
{
    grid: Grid<V>,
    indices: Vec<V>,
    level: usize,
    range: Range<usize>,
    throws: usize,
    success: usize,
    outside: Vec<V>,
}

impl<V> PoissonAlgo<V>
    where V: VecLike
{
    pub fn new<R>(poisson: &mut PoissonGen<R, V>) -> PoissonAlgo<V>
        where R: Rng
    {
        let dim = V::dim(None);
        let grid = Grid::new(poisson.radius, poisson.periodicity);
        let capacity = grid.cells() * dim;
        let mut indices = Vec::with_capacity(capacity);
        let choices = (0..grid.side).map(|i| i as f64).collect::<Vec<_>>();
        indices.extend(each_combination::<V>(&choices));
        let range = Range::new(0, indices.len());
        let throws = (poisson.a * indices.len() as f64).ceil() as usize;
        PoissonAlgo {
            grid: grid,
            indices: indices,
            level: 0,
            range: range,
            throws: throws,
            success: 0,
            outside: vec![],
        }
    }

    pub fn next<R>(&mut self, poisson: &mut PoissonGen<R, V>) -> Option<V>
        where R: Rng
    {
        while !self.indices.is_empty() && self.level < f64::MANTISSA_DIGITS as usize {
            while self.throws > 0 {
                self.throws -= 1;
                let index = self.range.ind_sample(&mut poisson.rand);
                let cur = self.indices[index];
                let parent = get_parent(cur, self.level);
                if !self.grid
                       .get(parent)
                       .expect("Indexing base grid by valid parent failed.")
                       .is_empty() {
                    self.indices.swap_remove(index);
                    if self.indices.is_empty() {
                        return None;
                    }
                    self.range = Range::new(0, self.indices.len());
                } else {
                    let sample = choose_random_sample(&mut poisson.rand,
                                                      &self.grid,
                                                      cur,
                                                      self.level);
                    if is_disk_free::<R, V>(&self.grid, poisson.radius, poisson.periodicity, cur, self.level, sample) && is_valid::<R, V>(poisson.radius, poisson.periodicity, &self.outside, sample) {
                        self.grid.get_mut(parent)
                            .expect("Indexing base grid by already indexed valid parent failed.")
                            .push(sample);
                        self.indices.swap_remove(index);
                        if !self.indices.is_empty() {
                            self.range = Range::new(0, self.indices.len());
                        }
                        self.success += 1;
                        return Some(sample);
                    }
                }
            }
            subdivide(&mut self.indices, &self.grid, &self.outside, &poisson, self.level);
            if self.indices.is_empty() {
                return None;
            }
            self.range = Range::new(0, self.indices.len());
            self.throws = (poisson.a * self.indices.len() as f64).ceil() as usize;
            self.level += 1;
        }
        None
    }

    pub fn size_hint<R>(&self, poisson: &PoissonGen<R, V>) -> (usize, Option<usize>)
        where R: Rng
    {
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and
        // how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dim(None);
        let side = 2usize.pow(self.level as u32);
        let spacing = self.grid.cell / side as f64;
        let grid_volume = self.indices.len() as f64 * spacing.powi(dim as i32);
        let sphere_volume = sphere_volume(2. * poisson.radius, dim as u64);
        let mut lower = (grid_volume / sphere_volume).floor() as usize;
        if lower > 0 {
            lower -= 1;
        }
        // Calculating upper bound should work because there is this many places left in the grid and no more can fit into it.
        let upper = self.grid.cells() - self.success;
        (lower, Some(upper))
    }

    pub fn insert(&mut self, sample: V) {
        //TODO: Figure out when the value should be returned by iterator if at all.
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

fn subdivide<R, V>(indices: &mut Vec<V>, grid: &Grid<V>, outside: &[V], poisson: &PoissonGen<R, V>, level: usize)
    where R: Rng,
          V: VecLike
{
    let choices = &[0., 1.];
    indices.flat_map_inplace(|i| {
        each_combination::<V>(choices)
            .map(move |n| n + i * 2.)
            .filter(|c| !covered(grid, poisson, outside, *c, level + 1))
    });
}

fn covered<R, V>(grid: &Grid<V>, poisson: &PoissonGen<R, V>, outside: &[V], index: V, level: usize) -> bool
    where R: Rng,
          V: VecLike
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell / side as f64;
    let sqradius = (2. * poisson.radius).powi(2);
    let parent = get_parent(index, level);
    each_combination(&[0., 1.])
        .map(|t| (index + t) * spacing)
        .all(|t| {
            each_combination(&[-2., -1., 0., 1., 2.])
                .filter_map(|t| grid.get(parent + t))
                .flat_map(|t| t)
                .any(|v| sqdist(*v, t, poisson.periodicity) < sqradius) ||
            !is_valid::<R, V>(poisson.radius, poisson.periodicity, &outside, t)
        })

}
