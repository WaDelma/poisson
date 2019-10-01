use poisson::Type::*;

extern crate nalgebra as na;
pub type Vect = na::Vector3<f64>;

mod helper;
use crate::helper::test_with_samples;

#[test]
fn test_3d_1_80_normal() {
    test_with_samples::<Vect>(1, 0.8, 1600, Normal);
}

#[test]
fn test_3d_1_80_perioditic() {
    test_with_samples::<Vect>(1, 0.8, 800, Perioditic);
}

#[test]
fn test_3d_10_80_normal() {
    test_with_samples::<Vect>(10, 0.8, 400, Normal);
}

#[test]
fn test_3d_10_80_perioditic() {
    test_with_samples::<Vect>(10, 0.8, 200, Perioditic);
}

#[test]
fn test_3d_100_80_normal() {
    test_with_samples::<Vect>(100, 0.8, 100, Normal);
}

#[test]
fn test_3d_100_80_perioditic() {
    test_with_samples::<Vect>(100, 0.8, 50, Perioditic);
}
