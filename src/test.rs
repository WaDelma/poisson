use ::{calc_radius, PoissonDisk, Sample};

use rand::{self, SeedableRng, XorShiftRng};

use na::Norm;
use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

#[test]
#[should_panic]
fn test_too_small_radius() {
    let _ = PoissonDisk::new(rand::weak_rng(), 0.0);
}

#[test]
#[should_panic]
fn test_too_large_radius() {
    let _ = PoissonDisk::new(rand::weak_rng(), 2f64.sqrt() / 2.0 + 0.0001);
}

#[test]
#[should_panic]
fn test_perioditic_too_small_radius() {
    let _ = PoissonDisk::perioditic(rand::weak_rng(), 0.0);
}

#[test]
#[should_panic]
fn test_perioditic_too_large_radius() {
    let _ = PoissonDisk::perioditic(rand::weak_rng(), 0.5);
}

#[test]
#[should_panic]
fn test_calc_radius_too_small_alpha() {
    calc_radius(1, -0.000001);
}

#[test]
#[should_panic]
fn test_calc_radius_too_small_points() {
    calc_radius(0, 0.5);
}

#[test]
fn test_calc_radius_smallest_alpha_calculates_valid() {
    let r = calc_radius(1, 0.000001);
    PoissonDisk::new(rand::weak_rng(), r);
}


#[test]
fn test_calc_radius_largest_alpha_calculates_valid() {
    let r = calc_radius(1, 1.0);
    PoissonDisk::new(rand::weak_rng(), r);
}


#[test]
#[should_panic]
fn test_calc_radius_too_large_alpha() {
    calc_radius(1, 1.000001);
}

#[test]
fn test_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds(radius, 1600, false);
}

#[test]
fn test_2nd_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    test_with_seeds(radius, 800, false);
}

#[test]
fn test_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    test_with_seeds(radius, 400, false);
}

#[test]
fn test_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    test_with_seeds(radius, 200, false);
}

#[test]
fn test_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    test_with_seeds(radius, 100, false);
}

#[test]
fn test_max_radius_periodic() {
    let radius = 0.499999999;
    test_with_seeds(radius, 1600, true);
}

#[test]
fn test_2nd_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    test_with_seeds(radius, 800, true);
}

#[test]
fn test_4th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    test_with_seeds(radius, 400, true);
}

#[test]
fn test_8th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    test_with_seeds(radius, 200, true);
}

#[test]
fn test_16th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    test_with_seeds(radius, 100, true);
}

#[test]
fn test_2th_of_max_radius_prefilled_with_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 2f64, 800, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius);
            poisson.create(v);
        });
}

#[test]
fn test_4th_of_max_radius_prefilled_with_2rd_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 4f64, 400, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 2f64);
            poisson.create(v);
        });
}

#[test]
fn test_8th_of_max_radius_prefilled_with_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 8f64, 200, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 4f64);
            poisson.create(v);
        });
}

#[test]
fn test_16th_of_max_radius_prefilled_with_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 16f64, 100, false, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::new(rand, radius / 8f64);
            poisson.create(v);
        });
}

#[test]
fn test_2th_of_max_radius_prefilled_with_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 2f64, 800, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, 0.499999999);
            poisson.create(v);
        });
}

#[test]
fn test_4th_of_max_radius_prefilled_with_2rd_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 4f64, 400, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 2f64);
            poisson.create(v);
        });
}

#[test]
fn test_8th_of_max_radius_prefilled_with_4th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 8f64, 200, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 4f64);
            poisson.create(v);
        });
}

#[test]
fn test_16th_of_max_radius_prefilled_with_8th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds_prefill(radius / 16f64, 100, true, &mut |ref mut v, i| {
            let rand = XorShiftRng::from_seed([i * 2 + 1, i * 1 + 1, i + 1, 2]);
            let mut poisson = PoissonDisk::perioditic(rand, radius / 8f64);
            poisson.create(v);
        });
}


fn test_with_seeds(radius: f64, seeds: u32, periodicity: bool) {
    test_with_seeds_prefill(radius, seeds, periodicity, &mut |_, _|{});
}

fn test_with_seeds_prefill<F>(radius: f64, seeds: u32, periodicity: bool, filler: &mut F) where F: FnMut(&mut Vec<Sample>, u32) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = if periodicity {
            PoissonDisk::perioditic(rand, radius)
        } else {
            PoissonDisk::new(rand, radius)
        };
        let mut vecs = vec![];
        filler(&mut vecs, i);
        poisson.create(&mut vecs);
        let vecs = if periodicity {
            let mut vecs2 = vec![];
            for x in &[-1.0, 0.0, 1.0] {
                for y in &[-1.0, 0.0, 1.0] {
                    for v in &vecs {
                        vecs2.push(*v + Vec2::new(*x, *y));
                    }
                }
            }
            vecs2
        } else {
            vecs
        };
        assert_legal_poisson(&vecs);
    }
}

fn assert_legal_poisson(vecs: &Vec<Sample>) {
    for &v1 in vecs {
        for &v2 in vecs {
            if v1.pos == v2.pos {
                continue;
            }
            let diff = v1.pos - v2.pos;
            let dist = diff.norm();
            let allowed_dist = v1.get_radius() + v2.get_radius();
            assert!(dist >= allowed_dist, "Poisson-disk distribution requirement not met: There exists 2 vectors with distance to each other of {} which is smaller than smallest allowed one {}. The samples: [{:?}, {:?}]", dist, allowed_dist, v1, v2);
        }
    }
}
