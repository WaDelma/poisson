extern crate poisson;
use poisson::{PoissonDisk};

extern crate rand;
use rand::{SeedableRng, XorShiftRng};

extern crate sphere;

extern crate nalgebra as na;
use na::Iterable;
pub type Vect = na::Vec2<f64>;

extern crate num;

use std::iter::repeat;

use helper::When::*;

mod helper;

#[test]
fn adding_valid_start_works() {
    let samples = 100;
    let relative_radius = 0.8;
    let rand = XorShiftRng::from_seed([0, 1, 1, 2]);
    let prefiller = |_| {
        let mut pre = PoissonDisk::new(rand.clone())
            .build_samples::<Vect>(samples, relative_radius)
            .into_iter()
            .take(25)
            .map(Some);
        move |_| pre.next().and_then(|s| s)
    };
    helper::test_with_samples_prefilled(samples, relative_radius, 100, false, prefiller, Always);
}

#[test]
fn adding_valid_middle_works() {
    let samples = 100;
    let relative_radius = 0.8;
    let rand = XorShiftRng::from_seed([1, 1, 2, 3]);
    let prefiller = |_| {
        let prefiller = PoissonDisk::new(rand.clone())
            .build_samples::<Vect>(samples, relative_radius);
        let mut pre = repeat(None)
            .take(25)
            .chain(prefiller
                .into_iter()
                .take(25)
                .map(Some));
        move |_| pre.next().and_then(|s| s)
    };

    helper::test_with_samples_prefilled(samples, relative_radius, 100, false, prefiller, Sometimes);
}

#[test]
fn adding_to_edges_start_works() {
    let samples = 100;
    let relative_radius = 0.8;
    let prefiller = [
        Vect::new(0.0, 0.0), Vect::new(0.0, 0.5),
        Vect::new(0.0, 1.0), Vect::new(0.5, 0.0),
        Vect::new(1.0, 0.0), Vect::new(0.5, 1.0),
        Vect::new(1.0, 0.5), Vect::new(1.0, 1.0),
        ];
    let prefiller = |_| {
        let mut pre = prefiller.iter().cloned().map(Some as fn(_) -> _);
        move |_| pre.next().and_then(|s| s)
    };
    helper::test_with_samples_prefilled(samples, relative_radius, 100, false, prefiller, Always);
}

#[test]
fn adding_to_outside_of_edges_start_works() {
    let samples = 100;
    let relative_radius = 0.8;
    let prefiller = [
        Vect::new(-0.1, -0.1), Vect::new(-0.1, 0.5),
        Vect::new(-0.1, 1.1), Vect::new(0.5, -0.1),
        Vect::new(1.1, -0.1), Vect::new(0.5, 1.1),
        Vect::new(1.1, 0.5), Vect::new(1.1, 1.1),
        ];
    let prefiller = |_| {
        let mut pre = prefiller.iter().cloned().map(Some as fn(_) -> _);
        move |_| pre.next().and_then(|s| s)
    };
    helper::test_with_samples_prefilled(samples, relative_radius, 100, false, prefiller, Always);
}
