use ::PoissonDisk;

use test::{test_with_seeds, test_with_seeds_prefill};

use rand::{SeedableRng, XorShiftRng};

use na::Vec4 as naVec4;
pub type Vec4 = naVec4<f64>;

#[test]
fn test_4d_max_radius_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds::<Vec4>(radius, 16, false);
}

#[test]
fn test_4d_2nd_max_radius_normal() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    test_with_seeds::<Vec4>(radius, 8, false);
}

#[test]
fn test_4d_4th_of_max_radius_normal() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    test_with_seeds::<Vec4>(radius, 4, false);
}

#[test]
fn test_4d_8th_of_max_radius_normal() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    test_with_seeds::<Vec4>(radius, 2, false);
}

#[test]
fn test_4d_16th_of_max_radius_normal() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    test_with_seeds::<Vec4>(radius, 1, false);
}

#[test]
fn test_4d_max_radius_perioditic() {
    let radius = 0.499999999;
    test_with_seeds::<Vec4>(radius, 4, true);
}

#[test]
fn test_4d_2nd_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    test_with_seeds::<Vec4>(radius, 2, true);
}

#[test]
fn test_4d_4th_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    test_with_seeds::<Vec4>(radius, 1, true);
}

#[test]
fn test_4d_8th_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    test_with_seeds::<Vec4>(radius, 1, true);
}

#[test]
fn test_4d_16th_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    test_with_seeds::<Vec4>(radius, 1, true);
}

#[test]
fn test_4d_2th_of_max_radius_prefilled_with_max_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 2f64, 8, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius);
            poisson.create(v);
        });
}

#[test]
fn test_4d_4th_of_max_radius_prefilled_with_2rd_of_max_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 4f64, 4, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 2f64);
            poisson.create(v);
        });
}

#[test]
fn test_4d_8th_of_max_radius_prefilled_with_4th_of_max_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 8f64, 2, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 4f64);
            poisson.create(v);
        });
}

#[test]
fn test_4d_16th_of_max_radius_prefilled_with_8th_of_max_normal() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 16f64, 1, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 8f64);
            poisson.create(v);
        });
}

#[test]
fn test_4d_2th_of_max_radius_prefilled_with_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 2f64, 2, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, 0.499999999);
            poisson.create(v);
        });
}

#[test]
fn test_4d_4th_of_max_radius_prefilled_with_2rd_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 4f64, 1, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 2f64);
            poisson.create(v);
        });
}

#[test]
fn test_4d_8th_of_max_radius_prefilled_with_4th_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 8f64, 1, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 4f64);
            poisson.create(v);
        });
}

#[test]
fn test_4d_16th_of_max_radius_prefilled_with_8th_of_max_radius_perioditic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill::<Vec4, _>(radius / 16f64, 1, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 8f64);
            poisson.create(v);
        });
}
