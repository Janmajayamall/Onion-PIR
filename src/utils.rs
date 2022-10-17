use num_bigint_dig::{BigUint, ModInverse};
use rand::{CryptoRng, RngCore};
use rand_distr::{
    num_traits::{FromPrimitive, PrimInt, ToPrimitive},
    Distribution, Normal, Uniform,
};
use std::mem::size_of;
pub fn sample_vec_cbd(size: usize, variance: isize) -> Vec<i64> {
    todo!()
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

// Returns multiplicative inverse of `v mod q`
pub fn mod_inverse(v: u64, q: u64) -> Option<u64> {
    let v = BigUint::from_u64(v).unwrap();
    let q = BigUint::from_u64(q).unwrap();
    v.mod_inverse(q)?.to_u64()
}

pub fn num_of_windows(q: u64, base: u64) -> u64 {
    ((q as f64).log2() / (base as f64).log2()).ceil() as u64
}
