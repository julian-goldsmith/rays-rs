extern crate cgmath;
extern crate png;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix2, Vector2, Vector3};
use png::HasParameters;

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Vector3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn point_at(&self, dist: f32) -> Vector3<f32> {
        self.origin + self.direction * dist
    }
}

// P(t) * P(t) = radius ^ 2 -- in local coordinates
// We want to find t
// dot(P - center, P - center) = radius ^ 2
#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Vector3<f32>,
    pub radius: f32,
}

impl Sphere {
    pub fn hit(&self, r: &Ray) -> Option<f32> {
        let r_origin_local = r.origin - self.center;
        let a = r.direction.dot(r.direction);
        let b = 2.0 * r_origin_local.dot(r.direction);
        let c = r_origin_local.dot(r_origin_local) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;

        if discriminant > 0.0 {
            Some((-b - discriminant.sqrt()) / (2.0 * a))
        } else {
            None
        }
    }
}

fn color(r: &Ray) -> Vector3<f32> {
    let unit_direction = r.direction.normalize();
    let t = 0.5 * (unit_direction.y + 1.0);
    (1.0 - t) * Vector3::new(1.0, 1.0, 1.0) + t * Vector3::new(0.5, 0.7, 1.0)
}

fn render(width: usize, height: usize) -> Vec<u8> {
    let mut pixels = Vec::new();
    let aspect = (height as f32) / (width as f32);
    let origin = Vector3::zero();
    let screen_space = Vector2::new(4.0, 4.0 * aspect);
    let lower_left_corner = (screen_space / -2.0).extend(-1.0);     // FIXME: Go from top left.

    let s = Sphere {
        center: Vector3::new(0.0, 0.0, -1.0),
        radius: 0.5,
    };

    let screen_space_transform = Matrix2::new(screen_space.x / (width as f32), 0.0, 0.0, screen_space.y / (height as f32));

    for y in 0..height {
        for x in 0..width {
            let uv = screen_space_transform * Vector2::new(x as f32, y as f32);

            let r = Ray {
                origin,
                direction: lower_left_corner + uv.extend(0.0),
            };

            let color = 255.0 * match s.hit(&r) {
                Some(dist) => {
                    let normal = (r.point_at(dist) - s.center).normalize();
                    0.5 * (normal + Vector3::new(1.0, 1.0, 1.0))
                },
                None => color(&r),
            };

            pixels.push(color.x as u8);
            pixels.push(color.y as u8);
            pixels.push(color.z as u8);
        };
    };

    pixels
}

fn write_png(path: &Path, data: &[u8], width: usize, height: usize) {
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width as u32, height as u32);
    encoder.set(png::ColorType::RGB).set(png::BitDepth::Eight);

    let mut writer = encoder.write_header().unwrap();
    writer.write_image_data(data).unwrap();
}

fn main() {
    let width = 900;
    let height = 700;
    let path = Path::new("/var/www/rays.png");

    let data = render(width, height);

    write_png(&path, &data, width, height);
}
