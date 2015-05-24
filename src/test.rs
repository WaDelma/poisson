use ::PoissonDisk;

use rand::{SeedableRng, XorShiftRng};

use na::Norm;
use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

#[test]
fn test_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    test_with_seeds(radius, 1600);
}

#[test]
fn test_2nd_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    test_with_seeds(radius, 800);
}

#[test]
fn test_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    test_with_seeds(radius, 400);
}

#[test]
fn test_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    test_with_seeds(radius, 200);
}

#[test]
fn test_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    test_with_seeds(radius, 100);
}

fn test_with_seeds(radius: f64, seeds: u32) {
    for i in 0..seeds {
        let rand = XorShiftRng::from_seed([i + 1, seeds - i + 1, (i + 1) * (i + 1), 1]);
        let mut poisson = PoissonDisk::new(rand, radius);
        let mut vecs = vec![];
        poisson.create(&mut vecs);
        assert_legal_poisson(&vecs, radius);
    }
}

fn assert_legal_poisson(vecs: &Vec<Vec2>, radius: f64) {
    for &v1 in vecs {
        for &v2 in vecs {
            if v1 == v2 {
                continue;
            }
            let diff = v1 - v2;
            let dist = diff.norm();
            assert!(dist >= 2f64 * radius, "Poisson-disk distribution requirement not met: There exists 2 vectors with distance to each other of {} which is smaller than smallest allowed one {}. The vectors: [{:?}, {:?}]", dist, 2f64 * radius, v1, v2);
        }
    }
}
