use poisson::Type::*;
use poisson::{algorithm, Builder};

use rand::{rngs::SmallRng, SeedableRng};

extern crate nalgebra as na;
pub type Vect = na::Vector2<f64>;

mod helper;
use crate::helper::test_with_samples;

#[test]
fn test_one_sample_works() {
    let rand = SmallRng::from_seed([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]);
    let builder = Builder::<_, Vect>::with_samples(1, 0.8, Normal);
    let builder = builder.build(rand, algorithm::Ebeida);
    builder.into_iter().for_each(drop);

    let rand = SmallRng::from_seed([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]);
    let builder = Builder::<_, Vect>::with_samples(1, 0.8, Normal);
    let builder = builder.build(rand, algorithm::Bridson);
    builder.into_iter().for_each(drop);
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
