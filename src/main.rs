extern crate cgmath;
extern crate png;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Vector2, Vector3};
use png::HasParameters;

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone)]
pub struct Ray {
    pub origin: Vector3<f32>,
    pub direction: Vector3<f32>,
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
    pub fn hit(&self, r: &Ray) -> bool {
        // t^2(D * D) + t(2E * D) + (E * E - 1) = 0
        //
        //
        // https://nelari.us/post/raytracer_with_rust_and_zig/
        // https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
        // http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf

        let r_origin_local = r.origin - self.center;
        let a = r.direction.dot(r.direction);
        let b = 2.0 * r_origin_local.dot(r.direction);
        let c = r_origin_local.dot(r_origin_local) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;

        discriminant > 0.0
    }
}

pub struct World {
    spheres: Vec<Sphere>,
}

fn color(r: &Ray) -> Vector3<f32> {
    let unit_direction = r.direction.normalize();
    let t = 0.5 * (unit_direction.y + 1.0);
    ((1.0 - t) * Vector3::new(1.0, 1.0, 1.0)) + (t * Vector3::new(0.5, 0.7, 1.0))
}

fn render(width: usize, height: usize) -> Vec<u8> {
    let mut pixels = Vec::new();
    let aspect = (width as f32) / (height as f32);
    let horizontal = Vector3::new(4.0, 0.0, 0.0);
    let vertical = Vector3::new(0.0, 4.0 * aspect, 0.0);
    let lower_left_corner = Vector3::new(horizontal.x / 2.0 - horizontal.x, -vertical.y / 2.0, -1.0);   // FIXME: Go from top left.
    let origin = Vector3::new(0.0, 0.0, 0.0);

    let s = Sphere {
        center: Vector3::new(0.0, 0.0, -1.0),
        radius: 0.5,
    };

    for y in 0..height {
        for x in 0..width {
            let u = (x as f32) / (width as f32);
            let v = (y as f32) / (height as f32);
            let r = Ray {
                origin,
                direction: lower_left_corner + (horizontal * u) + (vertical * v),
            };

            let color = if s.hit(&r) {
                Vector3::new(1.0, 0.0, 0.0)
            } else {
                color(&r)
            };

            pixels.push((color.x * 255.0) as u8);
            pixels.push((color.y * 255.0) as u8);
            pixels.push((color.z * 255.0) as u8);
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
    let data = render(900, 900);
    write_png(&Path::new("/var/www/rays.png"), &data, 900, 900);
}
