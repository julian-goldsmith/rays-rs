extern crate cgmath;
extern crate png;
extern crate rand;
extern crate rand_distr;

use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use cgmath::prelude::*;
use cgmath::{Matrix3, Matrix4, Point2, Point3, Vector2, Vector3};
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
    pub fn new(origin: Point3<f32>, direction: Vector3<f32>) -> Ray {
        Ray {
            origin,
            direction: direction.normalize(),
        }
    }

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
                    let jitter = (Vector3::from(UnitSphere.sample(&mut rng)) - Vector3::new(0.5, 0.5, 0.5)) * 0.5;
                    let tr = Ray::new(ixn.pos, ixn.normal + jitter);

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
    #[allow(dead_code)] persp: Matrix4<f32>,
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
        let color_factor = 255.999;
        let mut pixels = vec![0; 3 * self.size.x * self.size.y];

        for block_y in 0..(self.size.y / 16) {
            for block_x in 0..(self.size.x / 16) {
                let mut buf = [[Vector3::zero(); 16]; 16];
                self.render_block(block_x, block_y, &mut buf);

                for y in 0..16 {
                    for x in 0..16 {
                        let color = buf[y][x] * color_factor;

                        let actual_x = block_x * 16 + x;
                        let actual_y = block_y * 16 + y;
                        let pixel_base = 3 * (actual_y * self.size.x + actual_x);

                        pixels[pixel_base + 0] = color.x as u8;
                        pixels[pixel_base + 1] = color.y as u8;
                        pixels[pixel_base + 2] = color.z as u8;
                    };
                };
            };
        };

        self.write_png(&path, &pixels);
    }

    fn render_block(&self, block_x: usize, block_y: usize, buf: &mut [[Vector3<f32>; 16]; 16]) {
        let world = &self.world;
        let origin = self.view.transform_point(Point3::origin());

        for y in 0..16 {
            for x in 0..16 {
                let uv = Point2::new((block_x * 16 + x) as f32, (block_y * 16 + y) as f32);

                let mut color = Vector3::zero();

                for _ in 0..self.num_samples {
                    let jitter = Vector2::new(rand::random::<f32>() - 0.5, rand::random::<f32>() - 0.5);
                    let sample_uv = self.uvt.transform_point(uv + jitter) - Point2::origin();
                    let direction = self.pv.transform_vector(sample_uv.extend(0.1));                            // TODO: Don't hardcode near frustum distance.
                    let r = Ray::new(origin, direction);

                    color += world.sample(r, 0) / (self.num_samples as f32);
                };

                buf[y][x] = Vector3::new(color.x.sqrt(), color.y.sqrt(), color.z.sqrt());                       // Apply gamma.
            };
        };
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
        origin: Point3::new(0.0, 0.0, 3.0),                             // FIXME: Origin doesn't work properly (currently, it's flipped in Z).
        look_at: Point3::new(0.0, 0.0, -5.0),

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
