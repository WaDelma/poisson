//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//!    * For each point there is disk of certain radius which doesn't intersect with any disk of any other point
//!    * Nodes fill the space uniformly
//!

extern crate rand;
use rand::{Rand, Rng};
use rand::distributions::range::Range;
use rand::distributions::IndependentSample;

extern crate num;
use num::{Zero, One};

extern crate nalgebra as na;
use na::{Dim, Norm};

#[macro_use]
extern crate lazy_static;

/// Describes what traits the algorithm needs to be able to work.
pub trait VecLike:
    IndexMut<usize, Output = f64> +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<f64, Output = Self> +
    Div<f64, Output = Self> +
    Norm<f64> +
    PartialEq +
    Zero +
    One +
    Dim +
    Copy {}
impl<T> VecLike for T where T:
    IndexMut<usize, Output = f64> +
    Add<Output = T> +
    Sub<Output = T> +
    Mul<f64, Output = T> +
    Div<f64, Output = T> +
    Norm<f64> +
    PartialEq +
    Zero +
    One +
    Dim +
    Copy {}

use std::cmp::PartialEq;
use std::mem::drop;
use std::ops::{Sub, Mul, Add, Div, IndexMut};
use std::marker::PhantomData;

use tree::Node;

use math::Hypercube;
use math::Intersection::*;

mod tree;
mod math;
#[cfg(test)]
mod test;

// This means that last subdivide gives cubes which side is the smallest possible subnormal positive double.
const MAX_DEPTH: u32 = 1074;

/// Generates poisson-disk distribution in [0, 1]² area with O(N log N) time and space complexity relative to the number of samples generated.
/// Based on Gamito, Manuel N., and Steve C. Maddock. "Accurate multidimensional Poisson-disk sampling." ACM Transactions on Graphics (TOG) 29.1 (2009): 8.
///
/// # Examples
///
/// ```¨
/// extern crate poisson;
/// extern crate rand;
/// extern crate nalgebra as na;
/// type Vec2 = na::Vec2<f64>;
/// use poisson::PoissonDisk;
///
/// fn main() {
///     let mut poisson = PoissonDisk::new(rand::weak_rng()).build_radius::<Vec2>(0.1);
///     let mut vecs = vec![];
///     poisson.generate(&mut vecs);
///     println!("{:?}", vecs);
/// }
/// ```
pub struct PoissonDisk<R: Rng> {
    rand: R,
    periodicity: bool,
}

impl<R: Rng> PoissonDisk<R> {
    /// Creates new poisson-disk generator builder with random generator specified.
    pub fn new(rand: R) -> Self {
        PoissonDisk{rand: rand, periodicity: false}
    }

    /// Sets the generator to generate perioditic poisson-disk distributions.
    pub fn perioditic(mut self) -> Self {
        self.periodicity = true;
        self
    }

    /// Builds the generator with relative radius specified.
    /// Radius should be ]0, 1]
    pub fn build_relative_radius<V: VecLike>(self, radius: f64) -> PoissonGen<R, V> {
        assert!(0. < radius);
        assert!(radius <= 1.);
        PoissonGen {
            dim: PhantomData,
            radius: radius * (2f64.sqrt() / 2.),
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }

    /// Builds the generator with radius specified.
    /// Radius should be ]0, √2 / 2]
    pub fn build_radius<V: VecLike>(self, radius: f64) -> PoissonGen<R, V> {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        PoissonGen {
            dim: PhantomData,
            radius: radius,
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }

    /// Builds the generator with radius calculated so that approximately specified number of samples are generated.
    /// Amount of samples should be larger than 0.
    /// Relative radius should be [0, 1].
    /// For non-perioditic this is supported only for 2, 3 and 4 dimensional generation.
    pub fn build_samples<V: VecLike>(self, samples: u32, relative_radius: f64) -> PoissonGen<R, V> {
        assert!(self.periodicity || V::dim(None) < 5);
        assert!(samples > 0);
        assert!(relative_radius >= 0.);
        assert!(relative_radius <= 1.);
        PoissonGen {
            dim: PhantomData,
            radius: math::calc_radius::<V>(samples, relative_radius, self.periodicity),
            rand: self.rand,
            periodicity: self.periodicity,
        }
    }
}

pub struct PoissonGen<R: Rng, V: VecLike> {
    dim: PhantomData<V>,
    rand: R,
    radius: f64,
    periodicity: bool,
}

impl<R: Rng, V: VecLike> PoissonGen<R, V> {
    /// Populates given vector with poisson-disk distribution [0, 1]²
    /// Resulting samples will be a poisson-disk distribution iff given samples were already valid poisson-disk distribution.
    /// Resulting samples will be a maximal poisson-disk distribution [0, 1]² iff given samples have same radius and are already valid poisson-disk distribution.
    pub fn generate(&mut self, points: &mut Vec<Sample<V>>) {
        let dim = V::dim(None);
        let cell_width = self.radius / (dim as f64).sqrt();
        let cells_for_dim = (1. / cell_width) as usize;
        let mut grid = vec![None; cells_for_dim.pow(dim as u32)];
        let capacity = grid.len() * dim;
        let mut indices = Vec::with_capacity(capacity);
        indices.extend((0..grid.len()));
        let mut level = 0;
        let max_level = 63;
        let a = 0.3;
        while !indices.is_empty() && level < max_level {
            let mut range = Range::new(0, indices.len());
            for _ in 0..((a * grid.len() as f64) as usize) {
                let index = range.ind_sample(&mut self.rand);
                if let Some(_) = grid[Self::get_parent(index, level)] {
                    indices.swap_remove(index);
                    range = Range::new(0, indices.len());
                } else {
                    let c = self.choose_random_point(indices[index], level, cell_width);
                    if Self::is_disk_free(&grid, index, level, c) {
                        grid[Self::get_parent(index, level)] = Some(c);
                        indices.swap_remove(index);
                    }
                }
            }
            let mut n = 0;
            let mut len = indices.len();
            let mut added = 0;
            while n < len {
                let index = indices[n];
                let mut first = true;
                for child in Self::childs(index, level) {
                    if !self.covered(&grid, child, level + 1, cells_for_dim, cell_width) {
                        if first {
                            // If inserting first child we can just replace the parent
                            first = false;
                            indices[n] = child;
                        } else {
                            // Otherwise we just push the child to the end and keep tally how many times we have done this
                            indices.push(child);
                            added += 1;
                        }
                    }
                }
                if first {
                    // If all children were covered then we need to kill the parent
                    indices.swap_remove(n);
                    if added > 0 {
                        // If we have added children already to the end, we need to skip the child we replaced the parent with
                        added -= 1;
                        n += 1;
                    } else {
                        // Otherwise we just deal the replaced parent in the next iteration
                        len -= 1;
                    }
                } else {
                    // If even one child was added, we have already replaced the parent and we can move forwards
                    n += 1;
                }
            }
            level += 1;
            // If this assert fails then a is too small or subdivide code is broken
            assert_eq!(capacity, indices.capacity());
        }
	for sample in grid.iter().map(|v| Sample::new(v.unwrap(), self.radius)) {
	    points.push(sample);
	}
        /*let tree = Node::new(Hypercube::new(V::zero(), V::one()));
        for p in points.iter() {
            self.update_with_periodicity(&tree, None, *p);
        }
        while !tree.0.borrow().is_empty() {
            let (volume, sample) = self.generate_sample(&tree, None, 0);
            tree.reduce(None, volume);
            if let Some(s) = sample {
                self.update_with_periodicity(&tree, None, s);
                points.push(s);
            }
        }*/
    }

    /// Sets the radius of the generator.
    pub fn set_radius(&mut self, radius: f64) {
        assert!(0. < radius);
        assert!(radius <= (2f64.sqrt() / 2.));
        self.radius = radius;
    }

    /// Gets the radius of the generator.
    pub fn radius(&self) -> f64 {
        self.radius
    }
}

impl <R: Rng, V: VecLike> PoissonGen<R, V> {

    fn choose_random_point(&mut self, index: usize, level: usize, width: f64) -> V {
        let dim = V::dim(None);
        let cells_for_dim = 2f64.powf(-(level as f64)) as usize;
        let side = cells_for_dim as f64 * width;
        let mut t = Self::decode(index, cells_for_dim);
        for n in 0..dim {
            let place = f64::rand(&mut self.rand);
            t[n] = t[n] + place * side;//mul_add
        }
        t
    }
    
    fn is_disk_free(grid: &Vec<Option<V>>, index: usize, level: usize, c: V) -> bool {
        let to_be_checked = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1), (-1, -1), (-1, 1), (2, 0), (0, 2), (-2, 0), (0, -2), (2, 1), (2, -1), (-2, -1), (-2, 1), (1, 2), (-1, 2), (-1, -2), (1, -2)];
        for (a, b) in to_be_checked {
            //TODO
        }
        false
    }

    fn childs(index: usize, level: usize) -> Vec<usize> {
        let child_amount = 2usize.pow(V::dim(None) as u32);
        let mut childs = Vec::with_capacity(child_amount);
        for i in 0..child_amount {
	    childs.push(index * child_amount + i);
        }
        childs
    }

    fn covered(&self, grid: &Vec<Option<V>>, index: usize, level: usize, cells: usize, width: f64) -> bool {
        let parent = Self::get_parent(index, level);
        Self::for_each_combination(&[-1., 0., 1.], |t| {
            let index = parent + Self::encode(&t, cells);
            if let Some(s) = grid[index] {
                 if self.do_it(&s, index, width, level) {
                     return Some(true);
                 }
            }
        None
        }).unwrap_or(false)
    }

    fn do_it(&self, v: &V, index: usize, width: f64, level: usize) -> bool {
        let cells_for_dim = 2f64.powf(-(level as f64)) as usize;
        let side = cells_for_dim as f64 * width;
        let sqradius = self.radius.powi(2);
        let mut base = Self::decode(index, cells_for_dim);
        Self::for_each_combination(&[-1., 1.], |t| {
             let vec = (base + t) * side;
             if vec.sqnorm() > sqradius {
                 Some(false)
             } else {
                 None
             }
        }).unwrap_or(true)
    }

    fn for_each_combination<T, F: Fn(V) -> Option<T>>(choices: &[f64], func: F) -> Option<T> {
        let dim = V::dim(None);
        for n in 0..choices.len().pow(dim as u32) {
            let mut t = V::zero();
            let mut div = n;
            for n in 0..dim {
                let rem = div % choices.len();
                div /= choices.len();
                t[n] = choices[rem as usize];
            }
            if let s @ Some(_) = (func)(t) {
                return s;
            }
        }
        None
    }
    
    fn get_parent(index: usize, level: usize) -> usize {
	let split = 2usize.pow(V::dim(None) as u32);
	index / ((level + 1) * split)
    }

    fn encode(v: &V, level: usize) -> usize {
        let mut index = 0;
        for n in 0..V::dim(None) {
            index = (index + v[n] as usize) * level;
        }
        index / level
    }

    fn decode(index: usize, width: usize) -> V {
        let mut result = V::zero();
        let mut last = index;
        for n in (0..V::dim(None)).rev() {
            let cur = last / width;
            result[n] = (last - cur * width) as f64;
            last = cur;
        }
        result
    }
/*
    fn generate_sample(&mut self, node: &Node<V>, parent: Option<&Node<V>>, depth: u32) -> (f64, Option<Sample<V>>) {
        let borrow = node.0.borrow();
        if borrow.is_leaf() {
            let sample = Sample::new(borrow.cube.random_point_inside(&mut self.rand), self.radius);
            if borrow.is_sample_valid(sample) {
                (0., Some(sample))
            }else if depth < MAX_DEPTH {//borrow.cube.edge() > std::f64::MIN_POSITIVE {
                drop(borrow);
                (self.subdivide(node), None)
            } else {
                (borrow.volume, None)
            }
        } else {
            drop(borrow);
            let child = self.choose_random_child(&node);
            let result = self.generate_sample(&child, Some(node), depth + 1);
            child.reduce(parent, result.0);
            result
        }
    }

    fn choose_random_child(&mut self, node: &Node<V>) -> Node<V> {
        let borrow = node.0.borrow();
        let random_limit = f64::rand(&mut self.rand) * borrow.volume;
        let mut volume_counter = 0.;
        for child in &borrow.childs {
            volume_counter += child.0.borrow().volume;
            if volume_counter > random_limit {
                return child.clone();
            }
        }
        unreachable!("Either volumes of child nodes combined doesn't equal volume of the node or random doesn't generate number [0, 1[. This should never happen.");
    }

    fn subdivide(&self, node: &Node<V>) -> f64 {
        let mut delta = 0.;
        let mut childs = node.create_childs();
        let mut borrow = node.0.borrow_mut();
        childs.retain(|child| {
            let mut child_borrow = child.0.borrow_mut();
            for sample in &borrow.samples {
                match math::test_intersection(child_borrow.cube, sample.pos, sample.radius + self.radius) {
                    Over => child_borrow.samples.push(*sample),
                    In => {
                        delta += child_borrow.volume;
                        return false;
                    }
                    Out => {}
                }
            }
            true
        });
        borrow.childs = childs;
        borrow.samples.clear();
        delta
    }

    fn update_with_periodicity(&self, node: &Node<V>, parent: Option<&Node<V>>, sample: Sample<V>) {
        if self.periodicity {
            let choices = [-1., 0., 1.];
            let dim = V::dim(None);
            for n in 0..choices.len().pow(dim as u32) {
                let mut t = V::zero();
                let mut div = n;
                for i in 0..dim {
                    let rem = div % choices.len();
                    div /= choices.len();
                    t[i] = choices[rem as usize];
                }
                node.reduce(parent, self.update(node, Sample::new(sample.pos + t, sample.radius)));
            }
        } else {
            node.reduce(parent, self.update(node, sample));
        }
    }

    fn update(&self, node: &Node<V>, sample: Sample<V>) -> f64 {
        // TODO: Determine does this filtering help performance or not
        // if ((T::one() / 2.0) - sample.pos).sqnorm() > (0.5f64.sqrt() + sample.radius).powi(2) {
        //     return 0.0;
        // }
        let mut borrow = node.0.borrow_mut();
        match math::test_intersection(borrow.cube, sample.pos, sample.radius + self.radius) {
            Out => 0.,
            In => borrow.volume,
            Over => {
                let mut result = 0.;
                if borrow.is_leaf() {
                    borrow.samples.push(sample);
                } else {
                    borrow.childs.retain(|child| {
                        let delta = self.update(child, sample);
                        result += delta;
                        let mut child_borrow = child.0.borrow_mut();
                        if delta < child_borrow.volume {
                            child_borrow.volume -= delta;
                            true
                        } else {
                            false
                        }
                    });
                }
                result
            }
        }
    }*/
}

/// Describes position of sample and radius of disk around it.
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Sample<T: VecLike> {
    pub pos: T,
    radius: f64,
}

impl<T: VecLike> Sample<T> {
    pub fn new(pos: T, radius: f64) -> Self {
        Sample{pos: pos, radius: radius}
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }
}
