use crate::algorithm::{Algorithm, Creator};
use crate::utils::*;
use crate::{Builder, Float, Vector};

use num_traits::{Float as NumFloat};
use rand::distributions::{Distribution, Standard, Uniform};
use rand::Rng;

use sphere::sphere_volume;

/// Generates uniform maximal poisson-disk distribution with O(n2<sup>d</sup>) time and O(n2<sup>d</sup>) space complexity relative to the number of samples generated and the dimensionality of the sampling volume.
/// Based on Ebeida, Mohamed S., et al. "A Simple Algorithm for Maximal Poisson‐Disk Sampling in High Dimensions." Computer Graphics Forum. Vol. 31. No. 2pt4. Blackwell Publishing Ltd, 2012.
#[derive(Debug, Clone, Copy)]
pub struct Ebeida;

impl<F, V> Creator<F, V> for Ebeida
where
    F: Float,
    V: Vector<F>,
    Standard: Distribution<F>,
    Standard: Distribution<V>,
{
    type Algo = Algo<F, V>;

    fn create(poisson: &Builder<F, V>) -> Self::Algo {
        let dim = V::dimension();
        let grid = Grid::new(poisson.radius, poisson.poisson_type);
        let mut indices = Vec::with_capacity(grid.cells() * dim);
        let choices = (0..grid.side()).collect::<Vec<_>>();
        indices.extend(each_combination(&choices));
        let a = match dim {
            2 => 0.3,
            3 => 0.3,
            4 => 0.6,
            5 => 10.,
            6 => 700.,
            // TODO: Figure out what are optimal values beyond 6 dimensions
            _ => 700. + 100. * dim as f64,
        };
        Algo {
            a: a,
            grid: grid,
            throws: (a * indices.len() as f64).ceil() as usize,
            range: Uniform::new(0, indices.len()),
            indices: indices,
            level: 0,
            success: 0,
            outside: vec![],
            mantissa_digits: {
                let (mantissa, _, _) = <F as NumFloat>::max_value().integer_decode();
                mantissa.count_ones() as usize
            },
        }
    }
}

/// Implementation for the Ebeida algorithm
pub struct Algo<F, V>
where
    F: Float,
    V: Vector<F>,
{
    grid: Grid<F, V>,
    indices: Vec<V>,
    level: usize,
    range: Uniform<usize>,
    throws: usize,
    success: usize,
    outside: Vec<V>,
    mantissa_digits: usize,
    a: f64,
}

impl<F, V> Algorithm<F, V> for Algo<F, V>
where
    F: Float,
    V: Vector<F>,
    Standard: Distribution<F>,
    Standard: Distribution<V>,
{
    fn next<R>(&mut self, poisson: &mut Builder<F, V>, rng: &mut R) -> Option<V>
    where
        R: Rng,
    {
        if self.indices.is_empty() {
            return None;
        }
        while self.level < self.mantissa_digits {
            while self.throws > 0 {
                self.throws -= 1;
                let index = rng.sample(self.range);
                let cur = self.indices[index].clone();
                let parent = get_parent(cur.clone(), self.level);
                if !self
                    .grid
                    .get(parent.clone())
                    .expect("Indexing base grid by valid parent failed.")
                    .is_empty()
                {
                    self.indices.swap_remove(index);
                    if self.indices.is_empty() {
                        return None;
                    }
                    self.range = Uniform::new(0, self.indices.len());
                } else {
                    let sample = choose_random_sample(rng, &self.grid, cur.clone(), self.level);
                    if is_disk_free(
                        &self.grid,
                        poisson,
                        cur.clone(),
                        self.level,
                        sample.clone(),
                        &self.outside,
                    ) {
                        self.grid
                            .get_mut(parent)
                            .expect("Indexing base grid by already indexed valid parent failed.")
                            .push(sample.clone());
                        self.indices.swap_remove(index);
                        if !self.indices.is_empty() {
                            self.range = Uniform::new(0, self.indices.len());
                        }
                        self.success += 1;
                        return Some(sample);
                    }
                }
            }
            self.subdivide(&poisson);
            if self.indices.is_empty() {
                return None;
            }
            self.range = Uniform::new(0, self.indices.len());
            self.throws = (self.a * self.indices.len() as f64).ceil() as usize;
            self.level += 1;
        }
        let index = rng.sample(self.range);
        let cur = self.indices.swap_remove(index);
        let side = 2usize.pow(self.level as u32);
        let sample = index_to_sample(&cur, side);
        if is_disk_free(
            &self.grid,
            poisson,
            cur.clone(),
            self.level,
            sample.clone(),
            &self.outside,
        ) {
            Some(sample)
        } else {
            None
        }
    }

    fn size_hint(&self, poisson: &Builder<F, V>) -> (usize, Option<usize>) {
        // Calculating lower bound should work because we calculate how much volume is left to be filled at worst case and
        // how much sphere can fill it at best case and just figure out how many fills are still needed.
        let dim = V::dimension();
        let side = 2usize.pow(self.level as u32);
        let spacing = self.grid.cell() / F::cast(side);
        let grid_volume = F::cast(self.indices.len()) * NumFloat::powi(spacing, dim as i32);
        let sphere_volume = sphere_volume(F::cast(2) * poisson.radius, dim as u64);
        let lower = grid_volume / sphere_volume;
        let mut lower = NumFloat::floor(lower).to_usize().expect(
            "Grids volume divided by spheres volume should be always \
             castable to usize.",
        );
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

    fn stays_legal(&self, poisson: &Builder<F, V>, sample: V) -> bool {
        let index = sample_to_index(&sample, self.grid.side());
        is_disk_free(&self.grid, poisson, index, 0, sample.clone(), &self.outside)
    }
}

impl<F, V> Algo<F, V>
where
    F: Float,
    V: Vector<F>,
{
    fn subdivide(&mut self, poisson: &Builder<F, V>) {
        let choices = &[0, 1];
        let (grid, outside, level) = (&self.grid, &self.outside, self.level);
        self.indices.flat_map_inplace(|i| {
            each_combination(choices)
                .map(move |n: V| n + i.clone() * F::cast(2))
                .filter(|c| !covered(grid, poisson, outside, c.clone(), level + 1))
        });
    }
}

fn covered<F, V>(
    grid: &Grid<F, V>,
    poisson: &Builder<F, V>,
    outside: &[V],
    index: V,
    level: usize,
) -> bool
where
    F: Float,
    V: Vector<F>,
{
    // TODO: This does 4^d checking of points even though it could be done 3^d
    let side = 2usize.pow(level as u32);
    let spacing = grid.cell() / F::cast(side);
    let sqradius = NumFloat::powi(F::cast(2) * poisson.radius, 2);
    let parent = get_parent(index.clone(), level);
    each_combination(&[0, 1])
        .map(|t| (index.clone() + t) * spacing)
        .all(|t| {
            each_combination(&[-2, -1, 0, 1, 2])
                .filter_map(|t| grid.get(parent.clone() + t))
                .flat_map(|t| t)
                .any(|v| sqdist(v.clone(), t.clone(), poisson.poisson_type) < sqradius)
                || !is_valid(poisson, &outside, t)
        })
}
