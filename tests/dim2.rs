extern crate poisson;
use poisson::{Builder, algorithm};
use poisson::Type::*;

extern crate rand;
use rand::{XorShiftRng, SeedableRng};

extern crate sphere;

extern crate nalgebra as na;
pub type Vect = na::Vec2<f64>;

mod helper;
use helper::test_with_samples;

#[test]
fn test_one_sample_works() {
    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let builder = Builder::<_, Vect>::with_samples(1, 0.8, Normal);
    let builder = builder.build(rand, algorithm::Ebeida);
    builder.into_iter().collect::<Vec<Vect>>();

    let rand = XorShiftRng::from_seed([1, 2, 3, 4]);
    let builder = Builder::<_, Vect>::with_samples(1, 0.8, Normal);
    let builder = builder.build(rand, algorithm::Bridson);
    builder.into_iter().collect::<Vec<Vect>>();
}

#[test]
fn test_2d_1_80_normal() {
    test_with_samples::<Vect>(1, 0.8, 1600, Normal);
}

#[test]
fn test_2d_1_80_perioditic() {
    test_with_samples::<Vect>(1, 0.8, 800, Perioditic);
}

#[test]
fn test_2d_10_80_normal() {
    test_with_samples::<Vect>(10, 0.8, 800, Normal);
}

#[test]
fn test_2d_10_80_perioditic() {
    test_with_samples::<Vect>(10, 0.8, 400, Perioditic);
}

#[test]
fn test_2d_100_80_normal() {
    test_with_samples::<Vect>(100, 0.8, 400, Normal);
}

#[test]
fn test_2d_100_80_perioditic() {
    test_with_samples::<Vect>(100, 0.8, 200, Perioditic);
}
