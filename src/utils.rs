use num_bigint_dig::{prime::probably_prime, BigUint, ModInverse};
use num_traits::PrimInt;
use rand_distr::{
    num_traits::{FromPrimitive, ToPrimitive},
    Distribution, Uniform,
};
use std::{
    mem::size_of,
    ops::{Not, Shr, ShrAssign},
};

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

pub fn div_ceil<T: PrimInt>(a: T, b: T) -> T {
    assert!(b != T::zero());
    (a + b - T::one()) / b
}

// PRIME Functions
// Ref: https://github.com/tlepoint/fhe.rs/blob/a0287ba3842fcf19b45fd380c56ba7b5e52a387b/crates/fhe-util/src/lib.rs#L0-L1

/// Returns whether the modulus p is prime; this function is 100% accurate.
pub fn is_prime(p: u64) -> bool {
    probably_prime(&BigUint::from(p), 0)
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
