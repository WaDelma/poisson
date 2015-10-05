use image;
use VecLike;
use decode;

use std::fs::{self, File};
use std::path::{Path, PathBuf};

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
    let size = 256;
    let samples: Vec<V> = grid.iter().filter_map(|v| *v).collect();
    let mut imgbuf = image::ImageBuffer::new(size, size);

    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let half = image::Rgb([127 as u8, 127 as u8, 127 as u8]);
    let green = image::Rgb([0 as u8, 255 as u8, 0 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    let cells_per_cell = 2usize.pow(level as u32);
    let s = cells_per_cell * side;
    let spacing = top_lvl_cell / s as f64;
    for i in indices {
        let sample = decode::<V>(*i, s).unwrap();
        let xx = (sample[0] / s as f64 * size as f64) as i32;
        let yy = (sample[1] / s as f64 * size as f64) as i32;
        for x in 0..(spacing * side as f64 * size as f64) as i32 {
            for y in 0..(spacing * side as f64 * size as f64) as i32 {
                let cur_x = xx + x;
                let cur_y = yy + y;
                imgbuf.put_pixel(cur_x as u32, cur_y as u32, half);
            }
        }
    }
    for sample in samples {
        let xx = (sample[0] * size as f64) as i32;
        let yy = (sample[1] * size as f64) as i32;
        let radius = (r * 0.5 * size as f64) as i32;
        for x in -radius..(radius + 1) {
            for y in -radius..(radius + 1) {
                if x * x + y * y < radius * radius {
                    let color = if x == 0 || y == 0 {
                        middle
                    } else {
                        color
                    };
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
    let mut p = PathBuf::new();
    p.push("visualise");
    let _ = fs::create_dir(p.clone());
    p.push(format!("{}", level));
    let ref mut fout = File::create(p.with_extension("png")).unwrap();
    let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
}
