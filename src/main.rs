extern crate cgmath;
extern crate png;
extern crate rand;
extern crate rand_distr;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix3, Matrix4, Point2, Point3, Vector2, Vector3, Vector4};
use png::HasParameters;
use rand_distr::{Distribution, UnitSphere};

// https://nelari.us/post/raytracer_with_rust_and_zig/
// https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
// http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
// https://www.tutorialspoint.com/computer_graphics/3d_transformation.htm

// P(t) = E + tD, t >= 0
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Point3<f32>,
    pub direction: Vector3<f32>,
}

impl Ray {
    pub fn point_at_distance(&self, dist: f32) -> Point3<f32> {
        self.origin + self.direction * dist
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Intersection {
    pub pos: Point3<f32>,
    pub dist: f32,
    pub normal: Vector3<f32>,
    pub color: Vector3<f32>,
}

pub trait Intersect {
    fn intersect(&self, r: &Ray) -> Option<Intersection>;
}

#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    pub center: Point3<f32>,
    pub radius: f32,
    pub color: Vector3<f32>,
}

impl Intersect for Sphere {
    fn intersect(&self, r: &Ray) -> Option<Intersection> {
        let oc = r.origin - self.center;
        let b = -oc.dot(r.direction);
        let discriminant = b * b + self.radius * self.radius - oc.magnitude2();

        if discriminant > 0.0 {
            let dist = if b < discriminant.sqrt() {
                b + discriminant.sqrt()
            } else {
                b - discriminant.sqrt()
            };

            if dist > 0.0001 {
                let pos = r.point_at_distance(dist);

                Some(Intersection {
                    pos,
                    dist,
                    normal: (pos - self.center).normalize(),
                    color: self.color,
                })
            } else {
                None
            }
        } else {
            None
        }
    }
}

pub struct World {
    pub origin: Point3<f32>,
    pub look_at: Point3<f32>,

    pub spheres: Vec<Sphere>,
}

impl World {
    fn sample(&self, r: Ray, depth: usize) -> Vector3<f32> {
        let mut rng = rand::thread_rng();

        if depth > 4 {
            return Vector3::new(1.0, 0.0, 0.0);
        };

        let curr_ixn = self.spheres.iter().
            filter_map(|s| s.intersect(&r)).
            min_by(|a, b| a.dist.partial_cmp(&b.dist).unwrap());

        match curr_ixn {
            Some(ixn) if ixn.dist > 0.0001 => {
                let num_bounces = 4;
                let mut color = Vector3::zero();

                for _ in 0..num_bounces {
                    let jitter = Vector3::<f32>::from(UnitSphere.sample(&mut rng)) * 0.5;
                    let tr = Ray {
                        origin: ixn.pos,
                        direction: (ixn.normal + jitter - Vector3::new(0.5, 0.5, 0.5)).normalize(),
                    };

                    let bounce_sample = self.sample(tr, depth + 1);
                    color += ixn.color.mul_element_wise(bounce_sample);
                };

                color / num_bounces as f32
            },
            _ => {
                Vector3::new(0.6, 0.6, 0.6)
            },
        }
    }
}

pub struct Renderer {
    pub size: Vector2<usize>,
    pub aspect: f32,
    pub num_samples: usize,

    pub world: World,

    uvt: Matrix3<f32>,
    persp: Matrix4<f32>,
    view: Matrix4<f32>,
    pv: Matrix4<f32>,
}

impl Renderer {
    fn new(width: usize, height: usize, num_samples: usize, world: World) -> Renderer {
        let aspect = (width as f32) / (height as f32);
        let uvt = Matrix3::new(
            1.0 / width as f32, 0.0, 0.0,
            0.0, 1.0 / height as f32, 0.0,
            -0.5, -0.5, 1.0);

        let persp = cgmath::frustum(-0.5, 0.5, aspect * -0.5, aspect * 0.5, 0.1, 100.0);
        let view = Matrix4::look_at(world.origin, world.look_at, Vector3::new(0.0, 1.0, 0.0));

        Renderer {
            size: Vector2::new(width, height),
            aspect,
            num_samples,
            world,
            uvt,
            persp,
            view,
            pv: persp * view,
        }
    }

    fn render(&self, path: &Path) {
        let color_factor = 255.999 / (self.num_samples as f32).sqrt();
        let world = &self.world;
        let mut pixels = Vec::new();

        for y in (0..self.size.y).rev() {
            for x in 0..self.size.x {
                let uv = Point2::new(x as f32, y as f32);

                let mut color = Vector3::zero();

                for _ in 0..self.num_samples {
                    let uvj = Vector3::new(
                            uv.x + rand::random::<f32>() - 0.5,
                            uv.y + rand::random::<f32>() - 0.5, 1.0);
                    let sample_uv = self.uvt * uvj;
                    let direction = (self.pv * Vector4::new(sample_uv.x, sample_uv.y, 0.0, 1.0)).truncate();

                    let r = Ray {
                        origin: world.origin,
                        direction: direction.normalize(),
                    };

                    color += world.sample(r, 0);
                };

                let color = Vector3::new(
                    color.x.sqrt() * color_factor,
                    color.y.sqrt() * color_factor,
                    color.z.sqrt() * color_factor);
                pixels.push(color.x as u8);
                pixels.push(color.y as u8);
                pixels.push(color.z as u8);
            };
        };

        self.write_png(&path, &pixels);
    }

    fn write_png(&self, path: &Path, data: &[u8]) {
        let file = File::create(path).unwrap();
        let ref mut w = BufWriter::new(file);

        let mut encoder = png::Encoder::new(w, self.size.x as u32, self.size.y as u32);
        encoder.set(png::ColorType::RGB).set(png::BitDepth::Eight);

        let mut writer = encoder.write_header().unwrap();
        writer.write_image_data(data).unwrap();
    }
}

fn main() {
    let width = 1920;
    let height = 1080;
    let path = Path::new("rays.png");

    let world = World {
        origin: Point3::new(0.0, 0.0, 0.0),                             // FIXME: Origin doesn't work properly.
        look_at: Point3::new(0.0, 0.0, -10.0),

        spheres: vec![
            Sphere {
                center: Point3::new(0.0, 0.0, -5.0),
                radius: 0.5,
                color: Vector3::new(0.7, 0.3, 0.3),
            },
            Sphere {
                center: Point3::new(0.0, -100.5, -1.0),
                radius: 100.0,
                color: Vector3::new(0.1, 0.1, 0.9),
            },
        ],
    };
    let renderer = Renderer::new(width, height, 4, world);
    renderer.render(&path);
}
