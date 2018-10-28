extern crate poisson;
use poisson::Type;

extern crate rand;
use rand::{Rng, SeedableRng, XorShiftRng};
use rand::distributions::normal::StandardNormal;

extern crate sphere;

extern crate nalgebra as na;
pub type Vect = na::Vector2<f64>;

extern crate alga;
use self::alga::linear::FiniteDimVectorSpace;

extern crate num_traits;
use num_traits::Zero;

use helper::When::*;

mod helper;

#[test]
fn multiple_too_close_invalid() {
    let samples = 100;
    let relative_radius = 0.8;
    let prefiller = |radius| {
        let mut last = None::<Vect>;
        let mut rand = XorShiftRng::from_seed([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]);
        move |v| {
            if let Some(_) = v {
                if last == v {
                    None
                } else {
                    last = v;
                    let vec = sphere_uniform_point(&mut rand);
                    v.map(|v| v + vec * rand.gen::<f64>() * radius)
                }
            } else {
                None
            }
        }
    };
    // TODO: At 10 the test suddenly takes forever and takes all of the memory resulting into getting killed by oom killer
    helper::test_with_samples_prefilled(samples, relative_radius, 5, Type::Normal, prefiller, Never);
}

pub fn sphere_uniform_point<R: Rng>(rng: &mut R) -> Vect {
    let mut result = Vect::zero();
    for c in 0..Vect::dimension() {
        result[c] = rng.sample(StandardNormal);
    }
    result.normalize()
}
