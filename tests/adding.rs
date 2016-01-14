extern crate poisson;
use poisson::{PoissonDisk};

extern crate rand;
use rand::{SeedableRng, XorShiftRng};

extern crate sphere;

extern crate nalgebra as na;
pub type Vect = na::Vec2<f64>;

mod helper;

#[test]
fn adding_valid_at_start_works() {
    let samples = 100;
    let relative_radius = 0.8;
    let rand = XorShiftRng::from_seed([1, 1, 2, 3]);
    let prefiller = PoissonDisk::new(rand).build_samples::<Vect>(samples, relative_radius);
    helper::test_with_samples_prefilled(samples, relative_radius, 100, false, prefiller.into_iter().take(10));
}
