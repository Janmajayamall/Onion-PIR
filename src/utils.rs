use num_bigint::BigUint;
use num_bigint_dig::{prime::probably_prime, BigUint as BigUintDig, ModInverse};
use num_traits::{FromPrimitive, One, PrimInt, ToPrimitive, Zero};
use rand::{CryptoRng, RngCore};
use rand_distr::{Distribution, Uniform};
use std::{
    mem::size_of,
    ops::{Not, Shr, ShrAssign},
};

pub fn bit_vector(bit_size: usize, m: u64) -> Vec<bool> {
    let mut bit_vec = vec![false; bit_size];

    let mut value = m.clone();
    let mask = 1u64;
    let mut index = 0;
    while !value.is_zero() {
        let bit = &value & &mask;

        if bit.is_one() {
            bit_vec[index] = true;
        } else {
            bit_vec[index] = false;
        }

        value >>= 1;
        index += 1;
    }

    bit_vec
}

pub fn decompose_bits(m: u64, parent_bits: usize, decomp_bits: usize) -> Vec<u64> {
    let bit_vector = bit_vector(parent_bits, m);

    let decomp_m = bit_vector
        .chunks(decomp_bits)
        .into_iter()
        .map(|chunk| {
            chunk
                .iter()
                .rev()
                .fold(0u64, |acc, b| acc * 2 + (*b as u64))
        })
        .collect();

    decomp_m
}

/// Sample a vector of independent centered binomial distributions of a given
/// variance. Returns an error if the variance is strictly larger than 16.
///
/// Ref:https://github.com/tlepoint/fhe.rs/blob/a0287ba3842fcf19b45fd380c56ba7b5e52a387b/crates/fhe-util/src/lib.rs#L38
pub fn sample_vec_cbd<R: RngCore + CryptoRng>(
    vector_size: usize,
    variance: usize,
    rng: &mut R,
) -> Vec<i64> {
    assert!(variance >= 1 && variance <= 16);

    let mut out = Vec::with_capacity(vector_size);

    let number_bits = 4 * variance;
    let mask_add = ((u64::MAX >> (64 - number_bits)) >> (2 * variance)) as u128;
    let mask_sub = (mask_add << (2 * variance)) as u128;

    let mut current_pool = 0u128;
    let mut current_pool_nbits = 0;

    for _ in 0..vector_size {
        if current_pool_nbits < number_bits {
            current_pool |= (rng.next_u64() as u128) << current_pool_nbits;
            current_pool_nbits += 64;
        }
        debug_assert!(current_pool_nbits >= number_bits);
        out.push(
            ((current_pool & mask_add).count_ones() as i64)
                - ((current_pool & mask_sub).count_ones() as i64),
        );
        current_pool >>= number_bits;
        current_pool_nbits -= number_bits;
    }

    out
}

/// Returns number of bits `b` such that `2^b <= value`
///
/// @ref: https://github.com/tlepoint/fhe.rs/blob/27508002ea516d9ba41d0aa756bec7347f8404b2/crates/fhe-util/src/lib.rs#L193
pub fn ilog2<T: PrimInt>(value: T) -> usize {
    size_of::<T>() * 8 - 1 - value.leading_zeros() as usize
}

// Returns multiplicative inverse of `v mod q`
pub fn mod_inverse(v: u64, q: u64) -> Option<u64> {
    let v = BigUintDig::from_u64(v).unwrap();
    let q = BigUintDig::from_u64(q).unwrap();
    v.mod_inverse(q)?.to_u64()
}

pub fn div_ceil<T: PrimInt>(a: T, b: T) -> T {
    assert!(b != T::zero());
    (a + b - T::one()) / b
}

// PRIME Functions
// Ref: https://github.com/tlepoint/fhe.rs/blob/a0287ba3842fcf19b45fd380c56ba7b5e52a387b/crates/fhe-util/src/lib.rs#L0-L1

/// Returns whether the modulus p is prime; this function is 100% accurate.
pub fn is_prime(p: u64) -> bool {
    probably_prime(&BigUintDig::from(p), 0)
}

/// Generate a `num_bits`-bit prime, congruent to 1 mod `modulo`, strictly
/// smaller than `upper_bound`. Note that `num_bits` must belong to (10..=62),
/// and upper_bound must be <= 1 << num_bits.
pub fn generate_prime(num_bits: usize, modulo: u64, upper_bound: u64) -> Option<u64> {
    if !(10..=62).contains(&num_bits) {
        None
    } else {
        debug_assert!(
            (1u64 << num_bits) >= upper_bound,
            "upper_bound larger than number of bits"
        );

        let leading_zeros = (64 - num_bits) as u32;

        let mut tentative_prime = upper_bound - 1;
        while tentative_prime % modulo != 1 && tentative_prime.leading_zeros() == leading_zeros {
            tentative_prime -= 1
        }

        while tentative_prime.leading_zeros() == leading_zeros
            && !is_prime(tentative_prime)
            && tentative_prime >= modulo
        {
            tentative_prime -= modulo
        }

        if tentative_prime.leading_zeros() == leading_zeros && is_prime(tentative_prime) {
            Some(tentative_prime)
        } else {
            None
        }
    }
}

/// Unsigned 256-bit integer represented as four u64.
///
/// Ref - https://github.com/tlepoint/fhe.rs/blob/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d/crates/fhe-util/src/u256.rs#L0-L1
#[repr(C)]
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct U256(u64, u64, u64, u64);

impl U256 {
    /// Returns the additive identity element, 0.
    pub const fn zero() -> Self {
        Self(0, 0, 0, 0)
    }

    /// Add an U256 to self, wrapping modulo 2^256.
    pub fn wrapping_add_assign(&mut self, other: Self) {
        let (a, c1) = self.0.overflowing_add(other.0);
        let (b, c2) = self.1.overflowing_add(other.1);
        let (c, c3) = b.overflowing_add(c1 as u64);
        let (d, c4) = self.2.overflowing_add(other.2);
        let (e, c5) = d.overflowing_add((c2 | c3) as u64);
        let f = self.3.wrapping_add(other.3);
        let g = f.wrapping_add((c4 | c5) as u64);
        self.0 = a;
        self.1 = c;
        self.2 = e;
        self.3 = g;
    }

    /// Subtract an U256 to self, wrapping modulo 2^256.
    pub fn wrapping_sub_assign(&mut self, other: Self) {
        let (a, b1) = self.0.overflowing_sub(other.0);
        let (b, b2) = self.1.overflowing_sub(other.1);
        let (c, b3) = b.overflowing_sub(b1 as u64);
        let (d, b4) = self.2.overflowing_sub(other.2);
        let (e, b5) = d.overflowing_sub((b2 | b3) as u64);
        let f = self.3.wrapping_sub(other.3);
        let g = f.wrapping_sub((b4 | b5) as u64);
        self.0 = a;
        self.1 = c;
        self.2 = e;
        self.3 = g;
    }

    /// Returns the most significant bit of the unsigned integer.
    pub const fn msb(self) -> u64 {
        self.3 >> 63
    }
}

impl From<[u64; 4]> for U256 {
    fn from(a: [u64; 4]) -> Self {
        Self(a[0], a[1], a[2], a[3])
    }
}

impl From<[u128; 2]> for U256 {
    fn from(a: [u128; 2]) -> Self {
        Self(
            a[0] as u64,
            (a[0] >> 64) as u64,
            a[1] as u64,
            (a[1] >> 64) as u64,
        )
    }
}

impl From<U256> for [u64; 4] {
    fn from(a: U256) -> [u64; 4] {
        [a.0, a.1, a.2, a.3]
    }
}

impl From<U256> for [u128; 2] {
    fn from(a: U256) -> [u128; 2] {
        [
            (a.0 as u128) + ((a.1 as u128) << 64),
            (a.2 as u128) + ((a.3 as u128) << 64),
        ]
    }
}

impl From<&U256> for u128 {
    fn from(v: &U256) -> Self {
        debug_assert!(v.2 == 0 && v.3 == 0);
        (v.0 as u128) + ((v.1 as u128) << 64)
    }
}

impl Not for U256 {
    type Output = Self;

    fn not(self) -> Self {
        Self(!self.0, !self.1, !self.2, !self.3)
    }
}

impl Shr<usize> for U256 {
    type Output = Self;

    fn shr(self, rhs: usize) -> Self {
        let mut r = self;
        r >>= rhs;
        r
    }
}

impl ShrAssign<usize> for U256 {
    fn shr_assign(&mut self, rhs: usize) {
        debug_assert!(rhs < 256);

        if rhs >= 192 {
            self.0 = self.3 >> (rhs - 192);
            self.1 = 0;
            self.2 = 0;
            self.3 = 0;
        } else if rhs > 128 {
            self.0 = (self.2 >> (rhs - 128)) | (self.3 << (192 - rhs));
            self.1 = self.3 >> (rhs - 128);
            self.2 = 0;
            self.3 = 0;
        } else if rhs == 128 {
            self.0 = self.2;
            self.1 = self.3;
            self.2 = 0;
            self.3 = 0;
        } else if rhs > 64 {
            self.0 = (self.1 >> (rhs - 64)) | (self.2 << (128 - rhs));
            self.1 = (self.2 >> (rhs - 64)) | (self.3 << (128 - rhs));
            self.2 = self.3 >> (rhs - 64);
            self.3 = 0;
        } else if rhs == 64 {
            self.0 = self.1;
            self.1 = self.2;
            self.2 = self.3;
            self.3 = 0;
        } else if rhs > 0 {
            self.0 = (self.0 >> rhs) | (self.1 << (64 - rhs));
            self.1 = (self.1 >> rhs) | (self.2 << (64 - rhs));
            self.2 = (self.2 >> rhs) | (self.3 << (64 - rhs));
            self.3 >>= rhs;
        }
    }
}

mod tests {
    use super::*;

    #[test]
    fn test_decompose_bits() {
        let value = u64::MAX - 1;
        let decomposed = decompose_bits(value, 64, 4);
        let base = 1 << 4;
        let recovered = decomposed
            .iter()
            .rev()
            .fold(0u64, |acc, b| (acc * &base) + *b);

        assert_eq!(value, recovered);
    }

    #[test]
    fn test_bit_vector() {
        let m = u64::MAX;
        let bits = bit_vector(64, m);

        let value = bits
            .iter()
            .rev()
            .fold(0u64, |acc, bit| (acc * 2u64) + (*bit as u64));

        assert_eq!(value, m);
    }
}
