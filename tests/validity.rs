extern crate poisson;

extern crate rand;
use rand::{Rand, Rng, SeedableRng, XorShiftRng};

extern crate sphere;

extern crate nalgebra as na;
use na::{Norm, IterableMut};
pub type Vect = na::Vec2<f64>;

extern crate num;
use num::Zero;

use std::cell::RefCell;
use std::rc::Rc;
use std::mem::replace;

use helper::When::*;

mod helper;

#[test]
fn multiple_too_close_invalid() {
    let samples = 100;
    let relative_radius = 0.8;
    let prefiller = |radius| {
        let last = Rc::new(RefCell::new(None::<Vect>));
        let rand = Rc::new(RefCell::new(XorShiftRng::from_seed([9, 8, 7, 6])));
        move |v| {
            let mut borrow = last.borrow_mut();
            if let Some(_) = v {
                if *borrow == v {
                    None
                } else {
                    *borrow = v;
                    let rng = &mut *rand.borrow_mut();
                    v.map(|v| v + sphere_uniform_point(rng) * f64::rand(rng) * radius)
                }
            } else {
                None
            }
        }
    };
    helper::test_with_samples_prefilled(samples, relative_radius, 1000, false, prefiller, Never);
}

pub fn sphere_uniform_point<R: Rng>(rng: &mut R) -> Vect {
    let mut result = Vect::zero();
    for c in result.iter_mut() {
        replace(c, rand::distributions::normal::StandardNormal::rand(rng).0);
    }
    result.normalize()
}
