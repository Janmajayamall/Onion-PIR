use itertools::izip;
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use rand::{distributions::Uniform, CryptoRng, Rng, RngCore};
use std::{
    ops::{Add, AddAssign, Deref, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    sync::Arc,
};
#[derive(Clone, PartialEq, Debug)]
pub struct Modulus {
    pub p: u64,
    barret_hi: u64,
    barret_lo: u64,
}

/// Implements Modulus
///
/// TODO:
/// (1) Implement optimized lazy reduction (https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L654)
impl Modulus {
    pub fn new(p: u64) -> Self {
        assert!(p < 2 || (p >> 62 == 0));

        // `r` in barret reduction
        // 2**128 / p
        let r = ((BigUint::from(1u128) << 128usize) / p).to_u128().unwrap();
        Self {
            p,
            barret_hi: (r >> 64) as u64,
            barret_lo: r as u64,
        }
    }

    /// reduces `a` in [0, p ^ 2)
    /// to [0, 2p) in constant time
    ///
    /// x - ((r * x) / 2 ^ 2k) * p
    /// k = 64
    ///
    pub fn lazy_reduce_ct(&self, a: u64) -> u64 {
        let low = (a as u128 * (self.barret_lo as u128)) >> 64;
        let high = a as u128 * (self.barret_hi as u128);
        let num = (high + low) >> 64;

        let val = (a as u128) - (num * (self.p as u128));

        //TODO: Add debug asserts to check whether reduction works

        val as u64
    }

    pub fn lazy_reduce_u128(&self, a: u128) -> u64 {
        let alo = a as u64;
        let ahi = (a >> 64) as u64;

        let alo_lo = ((alo as u128) * (self.barret_lo as u128)) >> 64;
        let alo_hi = (alo as u128) * (self.barret_hi as u128);
        let ahi_lo = (ahi as u128) * (self.barret_lo as u128);
        let ahi_hi = (ahi as u128) * (self.barret_hi as u128);

        let num = ((alo_hi + ahi_lo + alo_lo) >> 64) + ahi_hi;
        let val = (a as u128) - (num * (self.p as u128));

        val as u64
    }

    /// `a` must be in the range [0, 2p)
    ///
    /// runs `a < p ? a : a - p`
    /// in constant time
    ///
    /// Ref:
    /// (1) https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L582
    /// (2) https://github.com/hacspec/rust-secret-integers/blob/master/src/lib.rs#L351-L366.
    pub const fn reduce_ct(x: u64, p: u64) -> u64 {
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        let (y, _) = x.overflowing_sub(p);
        let xp = x ^ p;
        let yp = y ^ p;
        let xy = xp ^ yp;
        let xxy = x ^ xy;
        let xxy = xxy >> 63;
        let (c, _) = xxy.overflowing_sub(1);
        let r = (c & y) | ((!c) & x);

        debug_assert!(r == x % p);
        r
    }

    pub fn reduce(&self, a: u64) -> u64 {
        Self::reduce_ct(self.lazy_reduce_ct(a), self.p)
    }

    pub fn reduce_u128(&self, a: u128) -> u64 {
        Self::reduce_ct(self.lazy_reduce_u128(a), self.p)
    }

    /// **Warning: this isn't constant time**
    pub fn reduce_biguint(&self, a: &BigUint) -> u64 {
        (a % self.p).to_u64().unwrap()
    }

    /// Ref - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L426
    fn reduce_i64(&self, a: i64) -> u64 {
        self.reduce_u128((((self.p as i128) << 64) + (a as i128)) as u128)
    }

    pub fn reduce_vec_u64(&self, a: &[u64]) -> Vec<u64> {
        a.iter().map(|ar| self.reduce(*ar)).collect_vec()
    }

    pub fn reduce_vec_i64(&self, a: &[i64]) -> Vec<u64> {
        a.iter().map(|ar| self.reduce_i64(*ar)).collect()
    }

    pub fn reduce_vec_biguint(&self, a: &[BigUint]) -> Vec<u64> {
        a.iter().map(|ar| self.reduce_biguint(ar)).collect()
    }

    /// Modulus exponentiation
    ///
    /// (a ** r) mod p
    pub fn pow(&self, a: u64, r: u64) -> u64 {
        debug_assert!(a < self.p && r < self.p);

        if r == 0 {
            return 1;
        } else if r == 1 {
            return a;
        }

        let mut bits = (62 - r.leading_zeros()) as isize;

        let mut val = a;
        while bits >= 0 {
            val = self.mul(val, val);

            if (r >> bits) & 1 == 1 {
                val = self.mul(val, a);
            }

            bits -= 1;
        }

        val
    }

    /// Modulus Inverse
    ///
    /// Remember p is prime. Therefore,
    /// a^(m-2) = a^(-1) mod m
    ///
    /// https://cp-algorithms.com/algebra/module-inverse.html#finding-the-modular-inverse-using-binary-exponentiation
    pub fn inv(&self, a: u64) -> u64 {
        debug_assert!(a < self.p && a != 0);
        self.pow(a, self.p - 2)
    }

    pub fn scalar_mul_vec(&self, a: &[u64], b: u64) -> Vec<u64> {
        let b_shoup = self.shoup(b);
        a.iter()
            .map(|ai| self.lazy_mul_shoup(*ai, b, b_shoup))
            .collect()
    }

    fn add(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p);
        debug_assert!(b < self.p);
        Self::reduce_ct(a + b, self.p)
    }

    fn sub(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p);
        debug_assert!(b < self.p);
        Self::reduce_ct(a + (self.p - b), self.p)
    }

    pub fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p);
        debug_assert!(b < self.p);
        self.reduce_u128((a as u128) * (b as u128))
    }

    pub fn neg(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        Self::reduce_ct(self.p - a, self.p)
    }

    pub fn add_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a, b).for_each(|(ab, b)| *ab = self.add(*ab, *b))
    }

    pub fn sub_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a, b).for_each(|(ab, b)| *ab = self.sub(*ab, *b))
    }

    pub fn mul_vec(&self, a: &mut [u64], b: &[u64]) {
        izip!(a, b).for_each(|(ab, b)| *ab = self.mul(*ab, *b))
    }

    /// Shoup representation of value
    ///
    /// (a * 2^64)/ p
    ///
    /// TODO: Understand math behind shoup repr
    pub fn shoup(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        (((a as u128) << 64) / (self.p as u128)) as u64
    }

    pub fn shoup_vec(&self, vals: &[u64]) -> Vec<u64> {
        vals.into_iter().map(|v| self.shoup(*v)).collect()
    }

    /// Lazy shoup multiplication of a, b in ct
    ///
    /// returns product in range [0, 2p)
    pub fn lazy_mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        debug_assert!(b < self.p);
        debug_assert!(b_shoup == self.shoup(b));

        // b_shoup = (b * 2^64) / p
        // q = (a * b_shoup) / 2^64 = (a * b * 2^64) / (p * 2^64) = (a * b) / p
        let q = ((a as u128) * (b_shoup as u128)) >> 64;
        let r = (((a as u128) * (b as u128)) - (q * (self.p as u128))) as u64;

        debug_assert!(r < self.p * 2);

        r
    }

    pub fn random_vec<R: RngCore + CryptoRng>(&self, size: usize, rng: &mut R) -> Vec<u64> {
        let uniform_dist = Uniform::from(0..self.p);
        rng.sample_iter(uniform_dist).take(size).collect()
    }

    pub fn modulus(&self) -> u64 {
        self.p
    }
}

#[derive(Clone, Debug)]
pub struct Context {
    pub moduli: Modulus,
    pub degree: usize,
}

impl Context {
    pub fn new(moduli: Modulus, degree: usize) -> Self {
        Self { moduli, degree }
    }
}

impl PartialEq for Context {
    fn eq(&self, other: &Self) -> bool {
        if self.moduli != other.moduli || self.degree != other.degree {
            return false;
        }
        true
    }
}

// #[derive(Clone, Debug)]
// /// Polynomial Ring: Rq = Zq[x] / (x^n + 1)
// pub struct Poly {
//     pub coeffs: Vec<u64>,
//     pub ctx: Arc<Context>,
// }

// // define operations
// impl Poly {
//     pub fn zero(ctx: &Arc<Context>) -> Self {
//         Poly {
//             coeffs: vec![0; ctx.degree],
//             ctx: ctx.clone(),
//         }
//     }

//     pub fn new(ctx: &Arc<Context>, coeffs: Vec<u64>) -> Self {
//         coeffs.iter().for_each(|v| assert!(v < &ctx.moduli.q));
//         Poly {
//             coeffs,
//             ctx: ctx.clone(),
//         }
//     }
// }

// impl AddAssign<&Poly> for Poly {
//     fn add_assign(&mut self, rhs: &Self) {
//         debug_assert!(self.ctx == rhs.ctx);

//         self.coeffs
//             .iter_mut()
//             .zip(rhs.coeffs.iter())
//             .for_each(|(a, b)| *a = self.ctx.moduli.add_mod(*a, *b));
//     }
// }

// impl Add<&Poly> for &Poly {
//     type Output = Poly;
//     fn add(self, rhs: &Poly) -> Self::Output {
//         let mut res = self.clone();
//         res += rhs;
//         res
//     }
// }

// impl Add for Poly {
//     type Output = Poly;

//     fn add(self, mut rhs: Poly) -> Self::Output {
//         rhs += &self;
//         rhs
//     }
// }

// impl SubAssign<&Poly> for Poly {
//     fn sub_assign(&mut self, rhs: &Poly) {
//         debug_assert!(self.coeffs.len() == rhs.coeffs.len());
//         debug_assert!(self.ctx == rhs.ctx);

//         self.coeffs
//             .iter_mut()
//             .zip(rhs.coeffs.iter())
//             .for_each(|(a, b)| *a = self.ctx.moduli.sub_mod(*a, *b));
//     }
// }

// impl Sub for Poly {
//     type Output = Poly;
//     fn sub(self, mut rhs: Poly) -> Self::Output {
//         rhs -= &self;
//         rhs
//     }
// }

// impl MulAssign<&Poly> for Poly {
//     fn mul_assign(&mut self, rhs: &Poly) {
//         debug_assert!(self.ctx == rhs.ctx);

//         let mut res = Poly::zero(&self.ctx);

//         for i in 0..res.ctx.degree {
//             for j in 0..i + 1 {
//                 let tmp = res.ctx.moduli.mul_mod(self.coeffs[j], rhs.coeffs[i - j]);
//                 res.coeffs[i] = res.ctx.moduli.add_mod(res.coeffs[i], tmp);
//             }

//             for j in i + 1..res.ctx.degree {
//                 let tmp = res
//                     .ctx
//                     .moduli
//                     .mul_mod(self.coeffs[j], rhs.coeffs[res.ctx.degree + i - j]);
//                 res.coeffs[i] = res.ctx.moduli.sub_mod(res.coeffs[i], tmp);
//             }

//             res.coeffs[i] %= res.ctx.moduli.q;
//         }

//         *self = res;
//     }
// }

// impl Mul<&Poly> for &Poly {
//     type Output = Poly;
//     fn mul(self, rhs: &Poly) -> Self::Output {
//         let mut res = self.clone();
//         res *= rhs;
//         res
//     }
// }

// impl Mul<Poly> for Poly {
//     type Output = Poly;
//     fn mul(self, mut rhs: Poly) -> Self::Output {
//         rhs *= &self;
//         rhs
//     }
// }

// impl Mul<u64> for &Poly {
//     type Output = Poly;
//     fn mul(self, rhs: u64) -> Self::Output {
//         assert!(rhs < self.ctx.moduli.q);
//         let mut res = Poly::zero(&self.ctx);
//         res.coeffs
//             .iter_mut()
//             .for_each(|v| *v = res.ctx.moduli.mul_mod(*v, rhs));
//         res
//     }
// }

// impl Mul<u64> for Poly {
//     type Output = Poly;
//     fn mul(self, rhs: u64) -> Self::Output {
//         &self * rhs
//     }
// }

// impl Neg for &Poly {
//     type Output = Poly;
//     fn neg(self) -> Self::Output {
//         Poly {
//             coeffs: self
//                 .coeffs
//                 .iter()
//                 .map(|v| self.ctx.moduli.neg(*v))
//                 .collect(),
//             ctx: self.ctx.clone(),
//         }
//     }
// }

// impl Neg for Poly {
//     type Output = Poly;
//     fn neg(self) -> Self::Output {
//         Poly {
//             coeffs: self
//                 .coeffs
//                 .iter()
//                 .map(|v| self.ctx.moduli.neg(*v))
//                 .collect(),
//             ctx: self.ctx.clone(),
//         }
//     }
// }

// impl Poly {
//     pub fn random_uniform(ctx: &Arc<Context>) -> Self {
//         let mut polynomial = Poly::zero(ctx);
//         polynomial
//             .coeffs
//             .as_mut_slice()
//             .copy_from_slice(&sample_uniform_vec(
//                 polynomial.ctx.degree,
//                 polynomial.ctx.moduli.q,
//             ));
//         polynomial
//     }

//     pub fn random_gaussian(ctx: &Arc<Context>, std_dev: f64) -> Self {
//         let mut polynomial = Poly::zero(ctx);
//         polynomial
//             .coeffs
//             .as_mut_slice()
//             .copy_from_slice(&Vec::<u64>::from_iter(
//                 sample_gaussian_vec(polynomial.ctx.degree, std_dev)
//                     .into_iter()
//                     .map(|mut v| {
//                         if v < 0.0 {
//                             v += polynomial.ctx.moduli.q as f64;
//                         }

//                         (v as u64) % polynomial.ctx.moduli.q
//                     }),
//             ));
//         polynomial
//     }

//     pub fn switch_context(&mut self, ctx: &Arc<Context>) {
//         self.ctx = ctx.clone();
//     }

//     pub fn decompose(&self, base: u64, l: u64) -> Vec<Poly> {
//         // base should be power of 2
//         debug_assert!((base & (base - 1)) == 0);

//         let decomposed_coeffs: Vec<Vec<u64>> = self
//             .coeffs
//             .iter()
//             .map(|v| decompose_value(*v, self.ctx.moduli.q, base, l))
//             .collect();

//         // Change poly from Rq => Rb
//         let poly_ctx = Arc::new(Context::new(Modulus { q: base }, self.ctx.degree));

//         (0..l as usize)
//             .into_iter()
//             .map(|d_index| {
//                 let mut p = Poly::zero(&poly_ctx);
//                 p.coeffs
//                     .iter_mut()
//                     .enumerate()
//                     .for_each(|(index, c)| *c = decomposed_coeffs[index][d_index]);
//                 p
//             })
//             .collect()
//     }

//     /// Transforms a poly of form `Σ b_i • X^i`
//     /// to `Σ b_i • X^(k*i)`
//     pub fn shift_powers(&self, k: u64) -> Self {
//         let mut res = Poly::zero(&self.ctx);
//         self.coeffs.iter().enumerate().for_each(|(power, c)| {
//             let new_power = k * power as u64;
//             let rounds = new_power / self.ctx.degree as u64;
//             let new_reduced_power = new_power % self.ctx.degree as u64;

//             // if `-1^rounds` is `+`
//             if rounds % 2 == 0 {
//                 res.coeffs[new_reduced_power as usize] = res
//                     .ctx
//                     .moduli
//                     .add_mod(res.coeffs[new_reduced_power as usize], *c);
//             } else {
//                 res.coeffs[new_reduced_power as usize] = res
//                     .ctx
//                     .moduli
//                     .sub_mod(res.coeffs[new_reduced_power as usize], *c)
//             }
//         });
//         res
//     }
// }

// /// Decomposes `value` into `window_size` bits.
// /// For `value` it returns [a0, a1, ..., ak]
// /// where each value `a{x}` is `window_size` bit
// /// values and `value = a0 * k^0 + a1 * k + a1 * k^2`
// /// where `k = 2 ** window_size`
// ///
// pub fn decompose_value(mut value: u64, moduli: u64, base: u64, l: u64) -> Vec<u64> {
//     assert!(value < moduli);

//     let q_bits = (moduli as f64).log2().to_u64().unwrap();
//     let base_bits = (base as f64).log2().to_u64().unwrap();
//     let precision_bits = l * base_bits;

//     assert!(q_bits >= precision_bits);

//     let mut bitvec = Vec::new();
//     for _ in 0..q_bits {
//         bitvec.push(value & 1);
//         value >>= 1;
//     }

//     // trim lsb
//     bitvec = bitvec[precision_bits as usize..].to_vec();

//     let decomp: Vec<_> = bitvec
//         .chunks(base_bits as usize)
//         .map(|chunk| chunk.iter().rev().fold(0, |acc, v| (acc << 1) + *v))
//         .collect();
//     decomp
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_mul() {
//         let ctx = Arc::new(Context::new(Modulus { q: 4 }, 4));
//         let p1 = Poly::new(&ctx, vec![1, 2, 3, 3]);
//         let p2 = Poly::new(&ctx, vec![1, 3, 3, 3]);
//         let product = p1 * p2;
//         assert!(product.coeffs == vec![1, 3, 3, 1]);
//     }

//     #[test]
//     fn test_decompose_value() {
//         let q = 65536;
//         let value = 2532;
//         let base = 4;
//         let l = 4;
//         let decomp = decompose_value(value, q, base, l);

//         let value1 = decomp
//             .iter()
//             .rev()
//             .fold(0, |acc, value| (acc * base) + value);
//         assert!(value == value1);
//     }

//     #[test]
//     fn test_decompose_poly() {
//         let ctx = Arc::new(Context::new(Modulus { q: 65536 }, 8));
//         let poly = Poly::zero(&ctx);

//         let d = poly.decompose(4, 4);
//     }
// }
