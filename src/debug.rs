#![allow(unused)]
use image;

use VecLike;
use utils::{encode, Grid};

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::ops::{Deref, DerefMut};

use modulo::Mod;

pub fn print_v<V: VecLike>(v: V) -> String {
    let mut result = "(".to_owned();
    for i in v.iter() {
        result.push_str(&format!("{}, ", i));
    }
    if V::dim(None) != 0 {
        result.pop();
    }
    result.push(')');
    result
}

pub fn visualise_3d<V: VecLike>(n: usize, level: usize, g: &Grid<V>, outside: &Vec<V>, indices: &Vec<V>, radius: f64) {
    let radius = 2. * radius;
    let periodicity = false;
    let samples: Vec<V> = g.data.iter().flat_map(|v| v).cloned().collect();
    let size = 2u32.pow(8);
    for z in 0..size {
        let mut imgbuf = image::ImageBuffer::new(size, size);
        let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
        let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
        for sample in &samples {
            let (xxx, yyy, zzz) = extract_3d(sample);
            let xx = (xxx * size as f64) as i32;
            let yy = (yyy * size as f64) as i32;
            let range = if periodicity {
                -1..(1 + 1)
            } else {
                0..1
            };
            for z_add in range {
                let zz = ((zzz + z_add as f64) * size as f64) as i32;
                let z2 = (zz - z as i32).pow(2);
                let radius = (radius * size as f64) as i32;
                for x in -radius..(radius + 1) {
                    for y in -radius..(radius + 1) {
                        if x * x + y * y + z2 < radius * radius {
                            let color = if x == 0 || y == 0 {
                                middle
                            } else {
                                color
                            };
                            let cur_x = if periodicity {
                                (xx + x).modulo(size as i32)
                            } else {
                                let cur_x = xx + x;
                                if cur_x < 0 || cur_x >= size as i32 {
                                    continue;
                                }
                                cur_x
                            };
                            let cur_y = if periodicity {
                                (yy + y).modulo(size as i32)
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
        let mut p = PathBuf::new();
        p.push("visualise");
        p.push("3d");
        let _ = fs::create_dir_all(p.clone());
        p.push(&*format!("{}_{}", n, z));
        let ref mut fout = File::create(p.with_extension("png")).unwrap();
        let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
    }
}

// pub fn visualise_3d<V: VecLike>(n: usize, level: usize, g: &Grid<V>, outside: &Vec<V>, indices: &Vec<V>, r: f64) {
//     let r = 2. * r;
//     let size = 2u32.pow(7);
//     let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
//     let samples: Vec<V> = g.data.iter().flat_map(|v| v).cloned().collect();
//     for zz in 0..size {
//         let mut imgbuf = image::ImageBuffer::new(size, size);
//         for sample in &samples {
//             let (x, y, z) = extract_3d(sample);
//             let xx = (x * size as f64) as i32;
//             let yy = (y * size as f64) as i32;
//             let dist = (z - (zz as f64 / size as f64)).abs();
//             if dist > r {
//                 continue;
//             }
//             let r = r * (1. - (dist / r).powi(2));
//             let radius = (r * size as f64) as i32;
//             draw_filled_circle(&mut imgbuf, size, xx, yy, radius, color);
//         }
//         let mut p = PathBuf::new();
//         p.push("visualise");
//         p.push("3d");
//         let _ = fs::create_dir(p.clone());
//         p.push(format!("{}_{}", n, zz));
//         let ref mut fout = File::create(p.with_extension("png")).unwrap();
//         let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
//     }
// }

pub fn visualise<V: VecLike>(n: usize, level: usize, g: &Grid<V>, outside: &Vec<V>, indices: &Vec<V>, r: f64) {
    let r = 2. * r;
    let size = 2u32.pow(8);//512;
    let samples: Vec<V> = g.data.iter().flat_map(|v| v).cloned().collect();
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let half = image::Rgb([127 as u8, 127 as u8, 127 as u8]);
    let purple = image::Rgb([127 as u8, 0 as u8, 255 as u8]);
    let yellow = image::Rgb([255 as u8, 255 as u8, 0 as u8]);
    let green = image::Rgb([0 as u8, 255 as u8, 0 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    let cells_per_cell = 2usize.pow(level as u32);

    let grid = (g.cell / cells_per_cell as f64) * size as f64;
    for i in indices {
        let (x, y) = extract(i);
        let x_start = (x * grid) as i32;
        let y_start = (y * grid) as i32;
        let x_end = ((x + 1.) * grid) as i32;
        let y_end = ((y + 1.) * grid) as i32;
        for x in x_start..x_end {
            for y in y_start..y_end {
                draw_pixel(&mut imgbuf, size, x, y, half);
            }
        }
    }

    let grid = (g.cell * size as f64) as u32;
    for x in 0..size {
        for y in 0..size {
            if x % grid == 0 || y % grid == 0 {
                let xx = x / grid;
                let yy = y / grid;
                let i = create(xx as f64, yy as f64);
                if let Some(s) = g.get(i) {
                    if !s.is_empty() {
                        let h = grid as u32 / 2;
                        imgbuf.put_pixel(xx * grid + h, yy * grid + h, yellow);
                    }
                }
                imgbuf.put_pixel(x as u32, y as u32, green);
            }
        }
    }

    for sample in samples {
        let (x, y) = extract(&sample);
        let xx = (x * size as f64) as i32;
        let yy = (y * size as f64) as i32;
        let radius = (r * size as f64) as i32;
        draw_pixel(&mut imgbuf, size, xx, yy, middle);
        draw_circle(&mut imgbuf, size, xx, yy, radius, color);
    }

    for sample in outside {
        let (x, y) = extract(sample);
        let xx = (x * size as f64) as i32;
        let yy = (y * size as f64) as i32;
        let radius = (r * size as f64) as i32;
        draw_pixel(&mut imgbuf, size, xx, yy, middle);
        draw_circle(&mut imgbuf, size, xx, yy, radius, purple);
    }
    let mut p = PathBuf::new();
    p.push("visualise");
    let _ = fs::create_dir(p.clone());
    p.push(format!("{}", n));
    let ref mut fout = File::create(p.with_extension("png")).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}

fn create<V: VecLike>(x: f64, y: f64) -> V {
    let mut i = V::zero();
    {
        let mut iter = i.iter_mut();
        *iter.next().unwrap() = x;
        *iter.next().unwrap() = y;
    }
    i
}

fn extract<V: VecLike>(v: &V) -> (f64, f64) {
    let mut iter = v.iter();
    (*iter.next().unwrap(), *iter.next().unwrap())
}

fn extract_3d<V: VecLike>(v: &V) -> (f64, f64, f64) {
    let mut iter = v.iter();
    (*iter.next().unwrap(), *iter.next().unwrap(), *iter.next().unwrap())
}

fn draw_filled_circle<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>, size: u32, x0: i32, y0: i32, radius: i32, color: image::Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
    for x in -radius..radius {
        for y in -radius..radius {
            if x.pow(2) + y.pow(2) < radius.pow(2) {
                draw_pixel(imgbuf, size, x + x0,  y + y0, color);
            }
        }
    }
}

fn draw_circle<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>, size: u32, x0: i32, y0: i32, radius: i32, color: image::Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
    let mut x = radius;
    let mut y = 0;
    let mut decision_over_2 = 1 - x;

    while y <= x {
        draw_pixel(imgbuf, size, x + x0,  y + y0, color);
        draw_pixel(imgbuf, size, y + x0,  x + y0, color);
        draw_pixel(imgbuf, size,-x + x0,  y + y0, color);
        draw_pixel(imgbuf, size,-y + x0,  x + y0, color);
        draw_pixel(imgbuf, size,-x + x0, -y + y0, color);
        draw_pixel(imgbuf, size,-y + x0, -x + y0, color);
        draw_pixel(imgbuf, size, x + x0, -y + y0, color);
        draw_pixel(imgbuf, size, y + x0, -x + y0, color);
        y += 1;
        if decision_over_2 <= 0 {
            decision_over_2 += 2 * y + 1;
        } else {
            x -= 1;
            decision_over_2 += 2 * (y - x) + 1;
        }
    }
}

fn draw_pixel<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>, size: u32, x: i32, y: i32, color: image::Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
    if x < 0 || x >= size as i32 {
        return;
    }
    if y < 0 || y >= size as i32 {
        return;
    }
    imgbuf.put_pixel(x as u32, y as u32, color);
}
