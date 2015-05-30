use ::{calc_radius, PoissonDisk};

use rand;

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
