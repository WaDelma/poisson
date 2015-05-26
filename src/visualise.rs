use ::PoissonDisk;

extern crate image;

extern crate rand;

use na::Vec2 as naVec2;
pub type Vec2 = naVec2<f64>;

use std::fs::File;
use std::path::Path;

#[test]
fn visualise_64th_radius() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_64th.png"));
}

#[test]
fn visualise_32th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_32th.png"));
}

#[test]
fn visualise_16th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_16th.png"));
}

#[test]
fn visualise_8th_of_max_radius() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::new(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_8th.png"));
}

#[test]
fn visualise_64th_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 64f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_64th_periodic.png"));
}

#[test]
fn visualise_32th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 32f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_32th_periodic.png"));
}

#[test]
fn visualise_16th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 16f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_16th_periodic.png"));
}

#[test]
fn visualise_8th_of_max_radius_periodic() {
    let radius = 2f64.sqrt() / 2f64 / 8f64;
    let size = 2 << 9;
    let mut poisson = PoissonDisk::perioditic(rand::weak_rng(), radius);
    let mut vecs = vec![];
    poisson.create(&mut vecs);
    visualise(&vecs, radius, size, &Path::new("poisson_8th_periodic.png"));
}


fn visualise(vecs: &Vec<Vec2>, radius: f64, size: u32, path: &Path) {
    let r = (radius * size as f64) as i32;
    let mut imgbuf = image::ImageBuffer::new(size, size);
    for vec in vecs {
        let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
        let xx = (vec.x * size as f64) as i32;
        let yy = (vec.y * size as f64) as i32;
        for x in -r..r {
            for y in -r..r {
                if x * x + y * y < r * r {
                    let cur_x = xx + x;
                    if cur_x < 0 || cur_x >= size as i32 {
                        continue;
                    }
                    let cur_y = yy + y;
                    if cur_y < 0 || cur_y >= size as i32 {
                        continue;
                    }
                    imgbuf.put_pixel(cur_x as u32, cur_y as u32, color);
                }
            }
        }
    }
    let ref mut fout = File::create(path).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}
