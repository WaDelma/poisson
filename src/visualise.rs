use ::PoissonDisk;
use ::math;

extern crate image;

extern crate rand;

use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

use std::fs::File;
use std::path::Path;
use std::collections::HashMap;

#[test]
fn visualise_64th_radius() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_64th.png"));
}

#[test]
fn visualise_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_32th.png"));
}

#[test]
fn visualise_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_16th.png"));
}

#[test]
fn visualise_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_8th.png"));
}

#[test]
fn visualise_4th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_4th.png"));
}

#[test]
fn visualise_2th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson_2th.png"));
}

#[test]
fn visualise_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, false, &Path::new("poisson.png"));
}

#[test]
fn visualise_64th_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_64th.png"));
}

#[test]
fn visualise_32th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_32th.png"));
}

#[test]
fn visualise_16th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_16th.png"));
}

#[test]
fn visualise_8th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_8th.png"));
}

#[test]
fn visualise_4th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 4f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_4th.png"));
}

#[test]
fn visualise_2th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson_2th.png"));
}

#[test]
fn visualise_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("periodic_poisson.png"));
}

#[test]
fn visualise_64th_of_max_radius_prefilled_with_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 32f64);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    let vecs2 = vecs.clone();
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius / 64f64);
    poisson.create(&mut vecs);
    let mut dists = HashMap::new();
    dists.insert((radius / 32f64 * size as f64) as i32, vecs2);
    dists.insert((radius / 64f64 * size as f64) as i32, vecs);
    visualise(&dists, size, true, &Path::new("prefilled_poisson_64th_32th.png"));
}

fn visualise(distributions: &HashMap<i32, Vec<Vec2>>, size: u32, periodicity: bool, path: &Path) {
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    for (radius, vecs) in distributions.iter() {
        for vec in vecs {
            let xx = (vec.x * size as f64) as i32;
            let yy = (vec.y * size as f64) as i32;
            for x in -*radius..*radius {
                for y in -*radius..*radius {
                    if x * x + y * y < radius * radius {
                        let color = if x == 0 || y == 0 {
                            middle
                        } else {
                            color
                        };
                        let cur_x = if periodicity {
                            math::modulo(xx + x, size as i32)
                        } else {
                            let cur_x = xx + x;
                            if cur_x < 0 || cur_x >= size as i32 {
                                continue;
                            }
                            cur_x
                        };
                        let cur_y = if periodicity {
                            math::modulo(yy + y, size as i32)
                        } else {
                            let cur_y = yy + y;
                            if cur_y < 0 || cur_y >= size as i32 {
                                continue;
                            }
                            cur_y
                        };
                        imgbuf.put_pixel(cur_x as u32, cur_y as u32, color);
                    }
                }
            }
        }
    }
    let ref mut fout = File::create(path).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}
