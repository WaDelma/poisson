use image;
use VecLike;
use decode;

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::ops::{Deref, DerefMut};

pub fn print_v<V: VecLike>(v: V) -> String {
    let mut result = "(".to_owned();
    for i in 0..V::dim(None) {
        result.push_str(&format!("{}, ", v[i]));
    }
    if V::dim(None) != 0 {
        result.pop();
    }
    result.push(')');
    result
}

pub fn visualise<V: VecLike>(level: usize, grid: &Vec<Option<V>>, side: usize, indices: &Vec<usize>, top_lvl_cell: f64, r: f64) {
    let size = 512;
    let samples: Vec<V> = grid.iter().filter_map(|v| *v).collect();
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let half = image::Rgb([127 as u8, 127 as u8, 127 as u8]);
    let green = image::Rgb([0 as u8, 255 as u8, 0 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    let cells_per_cell = 2usize.pow(level as u32);
    let s = cells_per_cell * side;
    for i in indices {
        let sample = decode::<V>(*i, s).unwrap();
        let x_start = (sample[0] / s as f64 * size as f64) as i32;
        let y_start = (sample[1] / s as f64 * size as f64) as i32;
        let x_end = ((sample[0] + 1.) / s as f64 * size as f64) as i32;
        let y_end = ((sample[1] + 1.) / s as f64 * size as f64) as i32;
        for x in x_start..x_end {
            for y in y_start..y_end {
                imgbuf.put_pixel(x as u32, y as u32, half);
            }
        }
    }
    let grid = (top_lvl_cell * size as f64) as u32;
    for x in 0..size {
        for y in 0..size {
            if x % grid == 0 || y % grid == 0 {
                imgbuf.put_pixel(x as u32, y as u32, green);
            }
        }
    }
    for sample in samples {
        let xx = (sample[0] * size as f64) as i32;
        let yy = (sample[1] * size as f64) as i32;
        let radius = (r * size as f64) as i32;
        draw_pixel(&mut imgbuf, size, xx, yy, middle);
        draw_circle(&mut imgbuf, size, xx, yy, radius, color);
        /*for x in -radius..(radius + 1) {
            for y in -radius..(radius + 1) {
                if x * x + y * y < radius * radius {
                    let color = if x == 0 || y == 0 {
                        middle
                    } else {
                        color
                    };
                    draw_pixel(&mut imgbuf, size, xx + x, yy + y, color);
                }
            }
        }*/
    }
    let mut p = PathBuf::new();
    p.push("visualise");
    let _ = fs::create_dir(p.clone());
    p.push(format!("{}", level));
    let ref mut fout = File::create(p.with_extension("png")).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}

fn draw_circle<C>(imgbuf: &mut image::ImageBuffer<image::Rgb<u8>, C>, size: u32, x0: i32, y0: i32, radius: i32, color: image::Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
    let mut x = radius;
    let mut y = 0;
    let mut decisionOver2 = 1 - x;

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
        if decisionOver2 <= 0 {
            decisionOver2 += 2 * y + 1;
        } else {
            x -= 1;
            decisionOver2 += 2 * (y - x) + 1;
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
