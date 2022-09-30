use std::mem::size_of;

use rand_distr::{num_traits::PrimInt, Distribution, Normal, Uniform};

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

/// Returns number of bits `b` such that `2^b <= value`
///
/// @ref: https://github.com/tlepoint/fhe.rs/blob/27508002ea516d9ba41d0aa756bec7347f8404b2/crates/fhe-util/src/lib.rs#L193
pub fn ilog2<T: PrimInt>(value: T) -> usize {
    size_of::<T>() * 8 - 1 - value.leading_zeros() as usize
}
