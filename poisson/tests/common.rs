use poisson::{Builder, Type};

extern crate nalgebra as na;
use crate::na::Vector2 as naVec2;
pub type Vec2 = naVec2<f64>;

#[test]
#[should_panic]
fn test_normal_too_small_radius() {
    let _ = Builder::<_, Vec2>::with_radius(0.0, Type::Normal);
}

#[test]
#[should_panic]
fn test_normal_too_large_radius() {
    let _ = Builder::<_, Vec2>::with_radius(2f64.sqrt() / 2.0 + 0.0001, Type::Normal);
}

// #[test]
// #[should_panic]
// fn test_calc_radius_too_small_alpha() {
//     calc_radius(1, -0.000001);
// }
//
// #[test]
// #[should_panic]
// fn test_calc_radius_too_small_points() {
//     calc_radius(0, 0.5);
// }
//
// #[test]
// fn test_calc_radius_smallest_alpha_calculates_valid() {
//     let r = calc_radius(1, 0.000001);
//     Builder::new(rand::weak_rng(), r);
// }
//
//
// #[test]
// fn test_calc_radius_largest_alpha_calculates_valid() {
//     let r = calc_radius(1, 1.0);
//     Builder::new(rand::weak_rng(), r);
// }
//
//
// #[test]
// #[should_panic]
// fn test_calc_radius_too_large_alpha() {
//     calc_radius(1, 1.000001);
// }
