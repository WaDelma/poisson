extern crate poisson;

extern crate rand;
use rand::{Rand, Rng, SeedableRng, XorShiftRng};

extern crate sphere;

extern crate nalgebra as na;
use na::{Norm, IterableMut};
pub type Vect = na::Vec2<f64>;

extern crate num;
use num::Zero;

use std::mem::replace;

use helper::When::*;

mod helper;

#[test]
fn multiple_too_close_invalid() {
    let samples = 100;
    let relative_radius = 0.8;
    let prefiller = |radius| {
        let mut last = None::<Vect>;
        let mut rand = XorShiftRng::from_seed([9, 8, 7, 6]);
        move |v| {
            if let Some(_) = v {
                if last == v {
                    None
                } else {
                    last = v;
                    let vec = sphere_uniform_point(&mut rand);
                    v.map(|v| v + vec * f64::rand(&mut rand) * radius)
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
