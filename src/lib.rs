//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//!    * For each point there is disk of certain radius which doesn't intersect with any disk of any other point
//!    * Nodes fill the space uniformly
//!

extern crate rand;
use rand::{Rand, Rng};

extern crate num;
use num::{Zero, One};

extern crate nalgebra as na;
use na::{Dim, Norm};
use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;
pub trait VecLike<T>:
    Index<usize, Output = f64> +
    IndexMut<usize, Output = f64> +
    Dim +
    Add<Output = T> +
    Sub<Output = T> +
    Norm<f64> +
    PartialEq +
    Div<f64, Output = T> +
    Zero +
    One +
    Copy +
    Clone +
    Debug {}
impl<T: Index<usize, Output = f64> + IndexMut<usize, Output = f64> + Dim + Add<Output = T> + Sub<Output = T> + Norm<f64> + PartialEq + Div<f64, Output = T> + Zero + One + Copy + Clone + Debug> VecLike<T> for T {}

use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::PartialEq;
use std::mem::drop;
use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;
use std::ops::{Sub, Add, Div, Index, IndexMut};

use math::{Intersection, Hypercube};

mod math;
#[cfg(test)]
mod test;

// This means that last subdivide gives cubes which side is the smallest possible subnormal positive double.
const MAX_DEPTH: u32 = 1074;

/// Calculates approximately needed radius from amount of points wanted and relative radius specified.
/// Points should be larger than 0.
/// Relative radius should be ]0, 1].
///
/// # Example
///
/// ```
/// extern crate rand;
/// extern crate poisson;
///
/// let radius = poisson::calc_radius(100, 0.8);
/// let mut poisson = poisson::PoissonDisk::new(rand::weak_rng(), radius);
/// let mut vecs = vec![];
/// poisson.create(&mut vecs);
/// println!("{:?}", vecs);
/// ```
pub fn calc_radius(points: u32, relative_radius: f64) -> f64 {
    assert!(0 < points);
    assert!(0f64 <= relative_radius);
    assert!(relative_radius <= 1f64);
    return relative_radius / (2f64 * 3f64.sqrt() * points as f64).sqrt();
}

/// Generates poisson-disk distribution in [0, 1]² area with O(N log N) time and space complexity relative to the number of samples generated.
/// Based on Gamito, Manuel N., and Steve C. Maddock. "Accurate multidimensional Poisson-disk sampling." ACM Transactions on Graphics (TOG) 29.1 (2009): 8.
///
/// # Examples
///
/// ```
/// extern crate rand;
/// extern crate poisson;
///
/// let mut poisson = poisson::PoissonDisk::new(rand::weak_rng(), 0.1);
/// let mut vecs = vec![];
/// poisson.create(&mut vecs);
/// println!("{:?}", vecs);
/// ```
pub struct PoissonDisk<R> where R: Rng {
    rand: R,
    radius: f64,
    periodicity: bool,
}

impl <R> PoissonDisk<R> where R: Rng {
    /// Creates new PoissonDisk with random generator and radius specified.
    /// Radius should be ]0, √2 / 2]
    pub fn new(rand: R, radius: f64) -> PoissonDisk<R> {
        assert!(0f64 < radius);
        assert!(radius <= (2f64.sqrt() / 2f64));
        PoissonDisk{rand: rand, radius: radius, periodicity: false}
    }


    /// Creates new peridiotic PoissonDisk with random generator, radius specified.
    /// Radius should be ]0, 0.5[
    ///
    /// # Examples
    ///
    /// ```
    /// extern crate rand;
    /// extern crate poisson;
    ///
    /// let mut poisson = poisson::PoissonDisk::perioditic(rand::weak_rng(), 0.1);
    /// let mut vecs = vec![];
    /// poisson.create(&mut vecs);
    /// println!("{:?}", vecs);
    /// ```
    pub fn perioditic(rand: R, radius: f64) -> PoissonDisk<R> {
        assert!(0f64 < radius);
        assert!(radius < 0.5);
        PoissonDisk{rand: rand, radius: radius, periodicity: true}
    }

    /// Populates given vector with poisson-disk distribution [0, 1]²
    /// Resulting samples will be a poisson-disk distribution iff given samples were already valid poisson-disk distribution.
    /// Resulting samples will be a maximal poisson-disk distribution [0, 1]² iff given samples have same radius and are already valid poisson-disk distribution.
    pub fn create<T: VecLike<T>>(&mut self, points: &mut Vec<Sample<T>>) {
        let tree = Node::new(None, Hypercube::new(T::zero(), T::one()));
        for p in points.iter() {
            self.update_with_periodicity(&tree, *p);
        }
        while !tree.0.borrow().is_empty() {
            let (volume, sample) = self.generate(&tree, 0);
            tree.reduce(volume);
            if let Some(s) = sample {
                self.update_with_periodicity(&tree, s);
                points.push(s);
            }
        }
    }

    fn generate<T: VecLike<T>>(&mut self, node: &Node<T>, depth: u32) -> (f64, Option<Sample<T>>) {
        let borrow = node.0.borrow();
        if borrow.is_leaf() {
            let sample = Sample::new(borrow.cube.random_point_inside(&mut self.rand), self.radius);
            if borrow.is_sample_valid(sample) {
                (0f64, Some(sample))
            }else if depth < MAX_DEPTH {
                drop(borrow);
                (self.subdivide(node), None)
            } else {
                (borrow.volume, None)
            }
        } else {
            drop(borrow);
            let child = self.choose_random_child(&node);
            let result = self.generate(&child, depth + 1);
            child.reduce(result.0);
            result
        }
    }

    fn choose_random_child<T: VecLike<T>>(&mut self, node: &Node<T>) -> Node<T> {
        let mut borrow = node.0.borrow_mut();
        let random_limit = f64::rand(&mut self.rand) * borrow.volume;
        let mut volume_counter = 0f64;
        for child in &mut borrow.childs {
            volume_counter += child.0.borrow().volume;
            if volume_counter > random_limit {
                return child.clone();
            }
        }
        unreachable!("Either volumes of child nodes combined doesn't equal volume of the node or random doesn't generate number [0, 1[. This should never happen.");
    }

    fn subdivide<T: VecLike<T>>(&self, node: &Node<T>) -> f64 {
		let mut delta = 0f64;
        let mut childs = node.create_childs();
        let mut borrow = node.0.borrow_mut();
        childs.retain(|child| {
            let mut child_borrow = child.0.borrow_mut();
        	for sample in &borrow.samples {
                match math::test_intersection(child_borrow.cube, sample.pos, sample.radius + self.radius) {
                    Intersection::Over => {
                        child_borrow.samples.push(*sample);
                    }
                    Intersection::In => {
    			        delta += child_borrow.volume;
    				    return false;
                    }
                    Intersection::Out => {}
                }
		    }
            true
        });
        borrow.childs = childs;
		borrow.samples.clear();
		return delta;
	}

    fn update_with_periodicity<T: VecLike<T>>(&self, node: &Node<T>, sample: Sample<T>) {
        if self.periodicity {
            let dim = T::dim(None);
            // println!("uwp:");
            for n in 0..3u64.pow(dim as u32) {
                let mut t = T::zero();
                for i in 0..dim {
                    // println!("{}:{}:{}", n, i, ((n / (1 + i as u64 * 2)) % 3));
                    t[i] = match (n / (1 + i as u64 * 2)) % 3 {
                        0 => -1.0,
                        1 => 0.0,
                        2 => 1.0,
                        j @ _ => unreachable!("This shouldn't be possible, but Rust cannot figure it out. Failed with: {}", j)
                    }
                }
                // println!("{:?}", t);
                node.reduce(self.update(node, sample + t));
            }
            // for x in &[-1.0, 0.0, 1.0] {
            //     for y in &[-1.0, 0.0, 1.0] {
            //         node.reduce(self.update(node, sample + Vec2::new(*x, *y)));
            //     }
            // }
        } else {
            node.reduce(self.update(node, sample));
        }
    }

    fn update<T: VecLike<T>>(&self, node: &Node<T>, sample: Sample<T>) -> f64 {
        let mut borrow = node.0.borrow_mut();
        match math::test_intersection(borrow.cube, sample.pos, sample.radius + self.radius) {
            Intersection::Out => {
                0.0
            }
            Intersection::In => {
                borrow.volume
            }
            Intersection::Over => {
                if borrow.is_leaf() {
                    borrow.samples.push(sample);
                    return 0.0;
                }
                let mut result = 0.0;
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
                result
            }
        }
    }
}

/// Describes position of sample and radius of disk around it.
#[derive(Clone, Copy)]
pub struct Sample<T: VecLike<T>> {
    pub pos: T,
    radius: f64,
}

impl<T: VecLike<T>> Sample<T> {
    pub fn new(pos: T, radius: f64) -> Self {
        Sample{pos: pos, radius: radius}
    }

    pub fn get_radius(&self) -> f64 {
        self.radius
    }
}

impl<T: VecLike<T>> Add<T> for Sample<T> {
    type Output = Sample<T>;
    fn add(self, other: T) ->  Self::Output {
        Sample{pos: self.pos + other, .. self}
    }
}

impl<T: VecLike<T>> Debug for Sample<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{:}:{:?}", self.radius, self.pos)
    }
}

#[derive(Clone)]
struct Node<T: VecLike<T>>(Rc<RefCell<InnerNode<T>>>);

struct InnerNode<T: VecLike<T>> {
    childs: Vec<Node<T>>,
    parent: Option<Node<T>>,
    samples: Vec<Sample<T>>,
    cube: Hypercube<T>,
    volume: f64,
}

impl<T: VecLike<T>> Debug for Node<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{:?}", self.0.borrow().childs)
    }
}

impl<T: VecLike<T>> Node<T> {
    fn new(parent: Option<Node<T>>, cube: Hypercube<T>) -> Node<T> {
        Node(Rc::new(RefCell::new(InnerNode{
                childs: vec![],
                samples: vec![],
                volume: cube.volume(),
                cube: cube,
                parent: parent
            })))
    }

    fn reduce(&self, amount: f64) {
        let mut borrow = self.0.borrow_mut();
        if amount < borrow.volume {
            borrow.volume -= amount;
        } else {
            let parent = borrow.parent.clone();
            borrow.childs.clear();
            borrow.samples.clear();
            borrow.parent = None;
            borrow.volume = 0f64;
            drop(borrow);
            if let Some(p) = parent {
                p.0.borrow_mut().childs.retain(|a| a.0.borrow().parent.is_some());
            }
        }
    }

    fn create_childs(&self) -> Vec<Node<T>> {
        let borrow = self.0.borrow();
        let mut result = vec![];
        let dim = T::dim(None);
        for n in 0..2u64.pow(dim as u32) {
            result.push(Node::new(Some(self.clone()), borrow.child_cube(n)));
        }
        result
        // vec![
        // Node::new(Some(self.clone()), borrow.child_cube(0, 0)),
        // Node::new(Some(self.clone()), borrow.child_cube(0, 1)),
        // Node::new(Some(self.clone()), borrow.child_cube(1, 0)),
        // Node::new(Some(self.clone()), borrow.child_cube(1, 1))]
    }
}

impl<T: VecLike<T>> PartialEq for InnerNode<T> {
    fn eq(&self, other: &Self) -> bool {
        self.cube == other.cube
    }
}

impl<T: VecLike<T>> InnerNode<T> {
    fn is_empty(&self) -> bool {
        self.childs.is_empty() && self.samples.is_empty() && self.volume == 0.0
    }

    fn is_leaf(&self) -> bool {
        self.childs.is_empty()
    }

    fn is_sample_valid(&self, sample: Sample<T>) -> bool {
        self.samples.iter().all(|s| (s.pos - sample.pos).sqnorm() >= (s.radius + sample.radius).powi(2))
    }

    fn child_cube(&self, n: u64) -> Hypercube<T> {
        Hypercube::new(self.calc_min(n), self.calc_max(n))
    }

    fn calc_min(&self, bits: u64) -> T {
        let center = self.cube.center();
        let mut t = T::zero();
        let dim = T::dim(None);
        for n in 0..dim {
            let bit = (bits >> n) & 1;
            t[n] = if bit == 0 {
                self.cube.min[n]
            } else {
                center[n]
            };
        }
        t
        // Vec2::new(
        //     if x == 0 {
        //         self.cube.min[0]
        //     }else{
        //         center[0]
        //     },
        //     if y == 0 {
        //         self.cube.min[1]
        //     }else{
        //         center[1]
        //     })
    }

    fn calc_max(&self, bits: u64) -> T {
        let center = self.cube.center();
        let mut t = T::zero();
        let dim = T::dim(None);
        for n in 0..dim {
            let bit = (bits >> n) & 1;
            t[n] = if bit == 0 {
                center[n]
            } else {
                self.cube.max[n]
            };
        }
        t
        // Vec2::new(
        //     if x == 0 {
        //         center[0]
        //     }else{
        //         self.cube.max[0]
        //     },
        //     if y == 0 {
        //         center[1]
        //     }else{
        //         self.cube.max[1]
        //     })
    }
}
