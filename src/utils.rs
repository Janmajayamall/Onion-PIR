use std::mem::take;

use rand_distr::{Distribution, Normal, Uniform};

pub fn sample_gaussian_vec(size: usize, std_dev: f64) -> Vec<f64> {
    let normal = Normal::new(0 as f64, std_dev).unwrap();
    let mut rng = rand::thread_rng();
    normal.sample_iter(&mut rng).take(size).collect()
}

pub fn sample_uniform_vec(size: usize, max: u64) -> Vec<u64> {
    let uniform = Uniform::new(0, max);
    let mut rng = rand::thread_rng();
    uniform.sample_iter(&mut rng).take(size).collect()
}
