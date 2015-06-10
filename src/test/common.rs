use ::PoissonDisk;

use rand;

use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

#[test]
#[should_panic]
fn test_normal_too_small_radius() {
    let _ = PoissonDisk::<_, Vec2>::with_radius(rand::weak_rng(), 0.0, false);
}

#[test]
#[should_panic]
fn test_normal_too_large_radius() {
    let _ = PoissonDisk::<_, Vec2>::with_radius(rand::weak_rng(), 2f64.sqrt() / 2.0 + 0.0001, false);
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
//     PoissonDisk::new(rand::weak_rng(), r);
// }
//
//
// #[test]
// fn test_calc_radius_largest_alpha_calculates_valid() {
//     let r = calc_radius(1, 1.0);
//     PoissonDisk::new(rand::weak_rng(), r);
// }
//
//
// #[test]
// #[should_panic]
// fn test_calc_radius_too_large_alpha() {
//     calc_radius(1, 1.000001);
// }
