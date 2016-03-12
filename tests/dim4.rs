extern crate poisson;
use poisson::Type::*;

extern crate rand;

extern crate sphere;

extern crate nalgebra as na;
pub type Vect = na::Vec4<f64>;

mod helper;
use helper::test_with_samples;

#[test]
fn test_4d_1_80_normal() {
    test_with_samples::<Vect>(1, 0.8, 200, Normal);
}

#[test]
fn test_4d_1_80_perioditic() {
    test_with_samples::<Vect>(1, 0.8, 100, Perioditic);
}

#[test]
fn test_4d_10_80_normal() {
    test_with_samples::<Vect>(10, 0.8, 50, Normal);
}

#[test]
fn test_4d_10_80_perioditic() {
    test_with_samples::<Vect>(10, 0.8, 15, Perioditic);
}

#[test]
fn test_4d_100_80_normal() {
    test_with_samples::<Vect>(100, 0.8, 10, Normal);
}

#[test]
fn test_4d_100_80_perioditic() {
    test_with_samples::<Vect>(100, 0.8, 1, Perioditic);
}
