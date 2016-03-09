use ::{PoissonDisk, VecLike, FloatLike};
use utils::*;
use algo::{AlgorithmCreator, Algorithm};

use num::Float;

use rand::Rng;
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

use sphere::sphere_volume;

use std::f64;

/// Generates poisson-disk distribution with O(N) time and space complexity relative to the number of samples generated.
/// Based on Ebeida, Mohamed S., et al. "A Simple Algorithm for Maximal Poisson‐Disk Sampling in High Dimensions." Computer Graphics Forum. Vol. 31. No. 2pt4. Blackwell Publishing Ltd, 2012.
#[derive(Debug, Clone, Copy)]
pub struct Ebeida;

impl<F, V> AlgorithmCreator<F, V> for Ebeida
    where F: FloatLike,
          V: VecLike<F>,
{
    type Algo = EbeidaAlgorithm<F, V>;

    fn create(poisson: &PoissonDisk<F, V>) -> Self::Algo {
        let dim = V::dim(None);
        let grid = Grid::new(poisson.radius, poisson.poisson_type);
        let capacity = grid.cells() * dim;
        let mut indices = Vec::with_capacity(capacity);
        let choices = (0..grid.side()).collect::<Vec<_>>();
        indices.extend(each_combination(&choices));
        let range = Range::new(0, indices.len());
        let a = match dim {
            2 => 0.3,
            3 => 0.3,
            4 => 0.6,
            5 => 10.,
            6 => 700.,
            // TODO: Figure out what are optimal values beyond 6 dimensions
            _ => 700. + 100. * dim as f64,
        };
        let throws = (a * indices.len() as f64).ceil() as usize;
        EbeidaAlgorithm {
            a: a,
            grid: grid,
            indices: indices,
            level: 0,
            range: range,
            throws: throws,
            success: 0,
            outside: vec![],
        }
    }
}

pub struct EbeidaAlgorithm<F, V>
    where F: FloatLike,
          V: VecLike<F>,
{
    grid: Grid<F, V>,
    indices: Vec<V>,
    level: usize,
    range: Range<usize>,
    throws: usize,
    success: usize,
    outside: Vec<V>,
    a: f64,
}

impl<F, V> Algorithm<F, V> for EbeidaAlgorithm<F, V>
    where F: FloatLike,
          V: VecLike<F>,
{
    fn next<R>(&mut self, poisson: &mut PoissonDisk<F, V>, rng: &mut R) -> Option<V>
        where R: Rng
    {
        while !self.indices.is_empty() && self.level < f64::MANTISSA_DIGITS as usize {
            while self.throws > 0 {
                self.throws -= 1;
                let index = self.range.ind_sample(rng);
                let cur = self.indices[index].clone();
                let parent = get_parent(cur.clone(), self.level);
                if !self.grid
                       .get(parent.clone())
                       .expect("Indexing base grid by valid parent failed.")
                       .is_empty() {
                    self.indices.swap_remove(index);
                    if self.indices.is_empty() {
                        return None;
                    }
                    self.range = Range::new(0, self.indices.len());
                } else {
                    let sample = choose_random_sample(rng,
                                                      &self.grid,
                                                      cur.clone(),
                                                      self.level);
                    if is_disk_free(&self.grid, poisson.radius, poisson.poisson_type, cur.clone(), self.level, sample.clone()) && is_valid(poisson.radius, poisson.poisson_type, &self.outside, sample.clone()) {
                        self.grid.get_mut(parent)
                            .expect("Indexing base grid by already indexed valid parent failed.")
                            .push(sample.clone());
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
            self.throws = (self.a * self.indices.len() as f64).ceil() as usize;
            self.level += 1;
        }
        None
    }

    fn size_hint(&self, poisson: &PoissonDisk<F, V>) -> (usize, Option<usize>) {
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and
        // how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dim(None);
        let side = 2usize.pow(self.level as u32);
        let spacing = self.grid.cell() / F::f(side);
        let grid_volume = F::f(self.indices.len()) * spacing.powi(dim as i32);
        let sphere_volume  = sphere_volume(F::f(2) * poisson.radius, dim as u64);
        let lower = grid_volume / sphere_volume;
        let mut lower = lower.floor().to_usize().unwrap();
        if lower > 0 {
            lower -= 1;
        }
        // Calculating upper bound should work because there is this many places left in the grid and no more can fit into it.
        let upper = self.grid.cells() - self.success;
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

    fn stays_legal(&self, poisson: &PoissonDisk<F, V>, sample: V) -> bool {
        let index = sample_to_index(&sample, self.grid.side());
        is_disk_free(&self.grid, poisson.radius, poisson.poisson_type, index, 0, sample.clone()) &&
        is_valid(poisson.radius, poisson.poisson_type, &self.outside, sample)
    }
}

fn subdivide<F, V>(indices: &mut Vec<V>, grid: &Grid<F, V>, outside: &[V], poisson: &PoissonDisk<F, V>, level: usize)
    where F: FloatLike,
          V: VecLike<F>,
{
    let choices = &[0, 1];
    indices.flat_map_inplace(|i| {
        each_combination(choices)
            .map(move |n: V| n + i.clone() * F::f(2))
            .filter(|c| !covered(grid, poisson, outside, c.clone(), level + 1))
    });
}

fn covered<F, V>(grid: &Grid<F, V>, poisson: &PoissonDisk<F, V>, outside: &[V], index: V, level: usize) -> bool
    where F: FloatLike,
          V: VecLike<F>,
{
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell() / F::f(side);
    let sqradius = (F::f(2) * poisson.radius).powi(2);
    let parent = get_parent(index.clone(), level);
    each_combination(&[0, 1])
        .map(|t| (index.clone() + t) * spacing)
        .all(|t| {
            each_combination(&[-2, -1, 0, 1, 2])
                .filter_map(|t| grid.get(parent.clone() + t))
                .flat_map(|t| t)
                .any(|v| sqdist(v.clone(), t.clone(), poisson.poisson_type) < sqradius) ||
            !is_valid(poisson.radius, poisson.poisson_type, &outside, t)
        })

}
