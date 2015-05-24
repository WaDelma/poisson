//! # Poisson-disk distribution generation
//!
//! Generates distribution of points where:
//!
//!    * For each point there is disk of certain radius which doesn't intercect with any disk of any other point
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

use math::{Intercection, Rect};

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
}

impl <R> PoissonDisk<R> where R: Rng {
    /// Creates new PoissonDisk with random generator and radius specified.
    /// Radius should be ]0, √2 / 2]
    pub fn new(rand: R, radius: f64) -> PoissonDisk<R> {
        assert!(0f64 < radius);
        assert!(radius <= (2f64.sqrt() / 2f64));
        PoissonDisk{rand: rand, radius: radius}
    }

    /// Populates given vector with poisson-disk distribution [0, 1]²
    pub fn create(&mut self, points: &mut Vec<Vec2>) {
        let tree = Node::new(None, Rect::new(Vec2::new(0f64, 0f64), Vec2::new(1f64, 1f64)));
        //TODO: Fill the tree with points
        while !tree.0.borrow().is_empty() {
            let (area, sample) = self.generate(&tree, 0);
            tree.reduce(area);
            if let Some(s) = sample {
                tree.reduce(self.update(&tree, s));
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
        let n = node.0.borrow();
        let diameter = 2f64 * self.radius;
        for s in &n.samples {
            let diff = *s - sample;
            if diff.sqnorm() < diameter * diameter {
                return false;
            }
        }
        true
	}

    fn subdivide(&self, node: &Node) -> f64 {
		let mut delta = 0f64;
        let mut borrow = node.0.borrow_mut();
    	let mut childs = vec![
        Node::new(Some(node.clone()), Rect::new(borrow.calc_min(0, 0), borrow.calc_max(0, 0))),
        Node::new(Some(node.clone()), Rect::new(borrow.calc_min(0, 1), borrow.calc_max(0, 1))),
        Node::new(Some(node.clone()), Rect::new(borrow.calc_min(1, 0), borrow.calc_max(1, 0))),
        Node::new(Some(node.clone()), Rect::new(borrow.calc_min(1, 1), borrow.calc_max(1, 1)))];
        childs.retain(|child| {
            let mut child_borrow = child.0.borrow_mut();
        	for sample in &borrow.samples {
                match math::test_intercection(child_borrow.rect, *sample, 2f64 * self.radius) {
                    Intercection::Over => {
                        child_borrow.samples.push(*sample);
                    }
                    Intercection::In => {
    			        delta += child_borrow.area;
    				    return false;
                    }
                    Intercection::Out => {}
                }
		    }
            true
        });
        borrow.childs = childs;
		borrow.samples.clear();
		return delta;
	}

	fn update(&self, node: &Node, sample: Vec2) -> f64 {
        let mut borrow = node.0.borrow_mut();
        match math::test_intercection(borrow.rect, sample, 2f64 * self.radius) {
            Intercection::Out => {
                0f64
            }
            Intercection::In => {
                // prune node and all its descendants
                // node.discard(); TODO
                borrow.area
            }
            Intercection::Over => {
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

    fn calc_min(&self, quadrant_x: u32, quadrant_y: u32) -> Vec2 {
        Vec2::new(
            if quadrant_x == 0 {
                self.rect.min.x
            }else{
                self.rect.center_x()
            },
            if quadrant_y == 0 {
                self.rect.min.y
            }else{
                self.rect.center_y()
            })
    }

    fn calc_max(&self, quadrant_x: u32, quadrant_y: u32) -> Vec2 {
        Vec2::new(
            if quadrant_x == 0 {
                self.rect.center_x()
            }else{
                self.rect.max.x
            },
            if quadrant_y == 0 {
                self.rect.center_y()
            }else{
                self.rect.max.y
            })
    }
}
