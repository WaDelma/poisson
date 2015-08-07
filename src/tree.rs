use std::cell::RefCell;
use std::fmt::{Debug, Formatter};
use std::fmt::Result as FmtResult;
use std::rc::Rc;
use math::Hypercube;
use {Sample, VecLike};

#[derive(Clone)]
pub struct Node<T: VecLike>(pub Rc<RefCell<InnerNode<T>>>);

struct InnerNode<T: VecLike> {
    pub childs: Vec<Node<T>>,
    pub samples: Vec<Sample<T>>,
    pub cube: Hypercube<T>,
    pub volume: f64,
}

impl<T: Debug + VecLike> Debug for Node<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{:?}", self.0.borrow().childs)
    }
}

impl<T: VecLike> Node<T> {
    pub fn new(cube: Hypercube<T>) -> Node<T> {
        Node(Rc::new(RefCell::new(InnerNode{
                childs: Vec::with_capacity(0),
                samples: Vec::with_capacity(0),
                volume: cube.volume(),
                cube: cube,
            })))
    }

    #[allow(unused_variables)]
    pub fn with_samples(cube: Hypercube<T>, samples: usize) -> Node<T> {
        Node(Rc::new(RefCell::new(InnerNode{
                childs: Vec::with_capacity(0),
                samples: Vec::with_capacity(samples),
                volume: cube.volume(),
                cube: cube,
            })))
    }

    //TODO: If there is only one child left for us we could compress by removing child and taking it's children? Max depth should then be checked via side and not by depth.
    pub fn reduce(&self, parent: Option<&Node<T>>, amount: f64) {
        let mut borrow = self.0.borrow_mut();
        if amount < borrow.volume {
            borrow.volume -= amount;
            // if borrow.childs.len() == 1 {
            //     let child = borrow.childs.remove(0);
            //     {
            //         let mut child_borrow = child.0.borrow_mut();
            //         if child_borrow.childs.len() > 0 {
            //             child_borrow.parent = None;
            //             for c in &mut child_borrow.childs {
            //                 let mut c_borrow = c.0.borrow_mut();
            //                 c_borrow.parent = Some(self.clone());
            //                 borrow.childs.push(c.clone());
            //             }
            //             child_borrow.childs.clear();
            //         }
            //     }
            //     if borrow.childs.is_empty() {
            //         borrow.childs.push(child);
            //     }
            // }
        } else {
            borrow.childs.clear();
            borrow.samples.clear();
            borrow.volume = 0.;
            drop(borrow);
            if let Some(p) = parent {
                p.0.borrow_mut().childs.retain(|a| a.0.borrow().volume > 0.);
            }
        }
    }

    pub fn create_childs(&self) -> Vec<Node<T>> {
        let borrow = self.0.borrow();
        let dim = T::dim(None);
        let childs = 2usize.pow(dim as u32);
        let mut result = Vec::with_capacity(childs);
        for n in 0..childs {
            result.push(Node::with_samples(borrow.child_cube(n), borrow.samples.len()));
        }
        result
    }
}

impl<T: VecLike> PartialEq for InnerNode<T> {
    fn eq(&self, other: &Self) -> bool {
        self.cube == other.cube
    }
}

impl<T: VecLike> InnerNode<T> {
    pub fn is_empty(&self) -> bool {
        self.is_leaf() && self.samples.is_empty() && self.volume == 0.
    }

    pub fn is_leaf(&self) -> bool {
        self.childs.is_empty()
    }

    pub fn is_sample_valid(&self, sample: Sample<T>) -> bool {
        self.samples.iter().all(|s| (s.pos - sample.pos).sqnorm() >= (s.radius + sample.radius).powi(2))
    }

    pub fn child_cube(&self, n: usize) -> Hypercube<T> {
        let center = self.cube.center();
        let min = Self::calc(n, self.cube.min, center);
        let max = Self::calc(n, center, self.cube.max);
        Hypercube::new(min, max)
    }

    fn calc(bits: usize, a: T, b: T) -> T {
        let mut t = T::zero();
        let dim = T::dim(None);
        for n in 0..dim {
            let bit = (bits >> n) & 1;
            t[n] = if bit == 0 { a } else { b }[n];
        }
        t
    }
}
