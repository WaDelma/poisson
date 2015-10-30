extern crate poisson;
use poisson::PoissonDisk;

extern crate rand;
use rand::{SeedableRng, XorShiftRng};

extern crate nalgebra as na;
pub type Vect = na::Vec3<f64>;

mod helper;
use helper::{test_with_samples, test_with_seeds_prefill};

#[test]
fn test_3d_1_80_normal() {
    test_with_samples::<Vect>(1, 0.8, 1600, false);
}

#[test]
fn test_3d_1_80_perioditic() {
    test_with_samples::<Vect>(1, 0.8, 800, true);
}

#[test]
fn test_3d_10_80_normal() {
    test_with_samples::<Vect>(10, 0.8, 800, false);
}

#[test]
fn test_3d_10_80_perioditic() {
    test_with_samples::<Vect>(10, 0.8, 400, true);
}

#[test]
fn test_3d_100_80_normal() {
    test_with_samples::<Vect>(100, 0.8, 400, false);
}

#[test]
fn test_3d_100_80_perioditic() {
    test_with_samples::<Vect>(100, 0.8, 200, true);
}

#[test]
#[ignore]
fn test_3d_2th_prefilled_1th_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vect, _>(radius / 2f64,
                                       800,
                                       false,
                                       &mut |ref mut v, i| {
                                           let rand = XorShiftRng::from_seed([i * 2 + 1,
                                                                              i * 1 + 1, i + 1, 2]);
                                           let mut poisson = PoissonDisk::new(rand)
                                                                 .build_radius(radius);
                                           poisson.generate(v);
                                       });
}

#[test]
#[ignore]
fn test_3d_8th_prefilled_4th_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vect, _>(radius / 8f64,
                                       100,
                                       false,
                                       &mut |ref mut v, i| {
                                           let rand = XorShiftRng::from_seed([i * 2 + 1,
                                                                              i * 1 + 1, i + 1, 2]);
                                           let mut poisson = PoissonDisk::new(rand)
                                                                 .build_radius(radius / 4f64);
                                           poisson.generate(v);
                                       });
}

#[test]
#[ignore]
fn test_3d_2th_prefilled_1th_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vect, _>(radius / 2f64,
                                       200,
                                       true,
                                       &mut |ref mut v, i| {
                                           let rand = XorShiftRng::from_seed([i * 2 + 1,
                                                                              i * 1 + 1, i + 1, 2]);
                                           let mut poisson = PoissonDisk::new(rand)
                                                                 .perioditic()
                                                                 .build_radius(0.499999999);
                                           poisson.generate(v);
                                       });
}

#[test]
#[ignore]
fn test_3d_8th_prefilled_4th_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vect, _>(radius / 8f64,
                                       25,
                                       true,
                                       &mut |ref mut v, i| {
                                           let rand = XorShiftRng::from_seed([i * 2 + 1,
                                                                              i * 1 + 1, i + 1, 2]);
                                           let mut poisson = PoissonDisk::new(rand)
                                                                 .perioditic()
                                                                 .build_radius(radius / 4f64);
                                           poisson.generate(v);
                                       });
}
