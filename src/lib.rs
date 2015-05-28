//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//!    * For each point there is disk of certain radius which doesn't intersect with any disk of any other point
//!    * Nodes fill the space uniformly
//!

extern crate rand;
use rand::Rng;

extern crate nalgebra as na;
use na::Norm;
use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::PartialEq;
use std::mem::drop;
use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;

use math::{Intersection, Rect};

mod math;
#[cfg(test)]
mod test;
#[cfg(test)]
#[cfg(feature = "visualise")]
mod visualise;

// This means that last subdivide gives rects which side is the smallest possible subnormal positive double.
const MAX_DEPTH: u32 = 1074;

/// Calculates approximately needed radius from amount of points wanted and relative radius specified.
/// Relative radius should be [0, 1].
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
    /// Radius should be ]0, √2 / 2]
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
        assert!(radius <= (2f64.sqrt() / 2f64));
        PoissonDisk{rand: rand, radius: radius, periodicity: true}
    }

    /// Populates given vector with poisson-disk distribution [0, 1]²
    pub fn create(&mut self, points: &mut Vec<Vec2>) {
        let tree = Node::new(None, Rect::new(Vec2::new(0f64, 0f64), Vec2::new(1f64, 1f64)));
        for p in points.iter() {
            self.update_with_periodicity(&tree, *p);
        }
        while !tree.0.borrow().is_empty() {
            let (area, sample) = self.generate(&tree, 0);
            tree.reduce(area);
            if let Some(s) = sample {
                self.update_with_periodicity(&tree, s);
                points.push(s);
            }
        }
    }

    fn generate(&mut self, node: &Node, depth: u32) -> (f64, Option<Vec2>) {
        let borrow = node.0.borrow();
        if borrow.is_leaf() {
            let sample = borrow.rect.random_point_inside(&mut self.rand);
            if self.is_sample_valid(node, sample) {
                (0f64, Some(sample))
            }else if depth < MAX_DEPTH {
                drop(borrow);
                (self.subdivide(node), None)
            } else {
                (borrow.area, None)
            }
        } else {
            drop(borrow);
            let child = self.choose_random_child(&node);
            let result = self.generate(&child, depth + 1);
            child.reduce(result.0);
            result
        }
    }

    fn choose_random_child(&mut self, node: &Node) -> Node {
        let mut borrow = node.0.borrow_mut();
        let random_limit = self.rand.gen::<f64>() * borrow.area;
        let mut area_counter = 0f64;
        for child in &mut borrow.childs {
            area_counter += child.0.borrow().area;
            if area_counter > random_limit {
                return child.clone();
            }
        }
        unreachable!("Either areas of child nodes combined doesn't equal area of the node or random doesn't generate number [0, 1[. This should never happen.");
    }

    fn is_sample_valid(&self, node: &Node, sample: Vec2) -> bool {
        let diameter = 2f64 * self.radius;
        let d2 = diameter * diameter;
        node.0.borrow().samples.iter().all(|s| (*s - sample).sqnorm() >= d2)
    }

    fn subdivide(&self, node: &Node) -> f64 {
		let mut delta = 0f64;
        let mut childs = node.create_childs();
        let mut borrow = node.0.borrow_mut();
        childs.retain(|child| {
            let mut child_borrow = child.0.borrow_mut();
        	for sample in &borrow.samples {
                match math::test_intersection(child_borrow.rect, *sample, 2f64 * self.radius) {
                    Intersection::Over => {
                        child_borrow.samples.push(*sample);
                    }
                    Intersection::In => {
    			        delta += child_borrow.area;
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

    fn update_with_periodicity(&self, node: &Node, sample: Vec2) {
        if self.periodicity {
            for x in &[-1, 0, 1] {
                for y in &[-1, 0, 1] {
                    node.reduce(self.update(node, sample + Vec2::new(*x as f64, *y as f64)));
                }
            }
        } else {
            node.reduce(self.update(node, sample));
        }
    }

	fn update(&self, node: &Node, sample: Vec2) -> f64 {
        let mut borrow = node.0.borrow_mut();
        match math::test_intersection(borrow.rect, sample, 2f64 * self.radius) {
            Intersection::Out => {
                0f64
            }
            Intersection::In => {
                // prune node and all its descendants
                // node.discard(); TODO
                borrow.area
            }
            Intersection::Over => {
                if borrow.is_leaf() {
        			borrow.samples.push(sample);
        			return 0f64;
        		}
        		let mut result = 0f64;
                borrow.childs.retain(|child| {
                	let delta = self.update(child, sample);
                    result += delta;
                    let mut child_borrow = child.0.borrow_mut();
        			if delta < child_borrow.area {
                        child_borrow.area -= delta;
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

#[derive(Clone)]
struct Node(Rc<RefCell<InnerNode>>);

struct InnerNode {
    childs: Vec<Node>,
    parent: Option<Node>,
    samples: Vec<Vec2>,
    rect: Rect,
    area: f64,
}

impl Debug for Node {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{:?}", self.0.borrow().childs)
    }
}

impl Node {
    fn new(parent: Option<Node>, rect: Rect) -> Node {
        Node(Rc::new(RefCell::new(InnerNode{
                childs: vec![],
                samples: vec![],
                area: rect.area(),
                rect: rect,
                parent: parent
            })))
    }

    fn reduce(&self, amount: f64) {
        let mut borrow = self.0.borrow_mut();
        if amount < borrow.area {
            borrow.area -= amount;
        } else {
            let parent = borrow.parent.clone();
            borrow.childs.clear();
            borrow.samples.clear();
            borrow.parent = None;
            borrow.area = 0f64;
            drop(borrow);
            if let Some(p) = parent {
                p.0.borrow_mut().childs.retain(|a| a.0.borrow().parent.is_some());
            }
        }
    }

    fn create_childs(&self) -> Vec<Node> {
        let borrow = self.0.borrow();
        vec![
        Node::new(Some(self.clone()), borrow.child_rect(0, 0)),
        Node::new(Some(self.clone()), borrow.child_rect(0, 1)),
        Node::new(Some(self.clone()), borrow.child_rect(1, 0)),
        Node::new(Some(self.clone()), borrow.child_rect(1, 1))]
    }
}

impl PartialEq for InnerNode {
    fn eq(&self, other: &Self) -> bool {
        self.rect == other.rect
    }
}

impl InnerNode {
    fn is_empty(&self) -> bool {
        self.childs.is_empty() && self.samples.is_empty() && self.area == 0f64
    }

    fn is_leaf(&self) -> bool {
        self.childs.is_empty()
    }

    fn child_rect(&self, x: u32, y: u32) -> Rect {
        Rect::new(self.calc_min(x, y), self.calc_max(x, y))
    }

    fn calc_min(&self, x: u32, y: u32) -> Vec2 {
        Vec2::new(
            if x == 0 {
                self.rect.min.x
            }else{
                self.rect.center_x()
            },
            if y == 0 {
                self.rect.min.y
            }else{
                self.rect.center_y()
            })
    }

    fn calc_max(&self, x: u32, y: u32) -> Vec2 {
        Vec2::new(
            if x == 0 {
                self.rect.center_x()
            }else{
                self.rect.max.x
            },
            if y == 0 {
                self.rect.center_y()
            }else{
                self.rect.max.y
            })
    }
}
