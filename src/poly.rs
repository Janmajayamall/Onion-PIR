use super::utils::{sample_gaussian_vec, sample_uniform_vec};
use std::{
    clone,
    intrinsics::assert_inhabited,
    ops::{Add, AddAssign, Deref, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    sync::Arc,
};
#[derive(Clone, PartialEq, Debug)]
pub struct Modulus {
    pub q: u64,
}

impl Modulus {
    fn add_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a + b) % self.q
    }

    fn sub_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a + self.q - b) % self.q
    }

    fn mul_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a * b) % self.q
    }

    fn neg(&self, a: u64) -> u64 {
        debug_assert!(a <= self.q);
        (self.q - a) % self.q
    }

    pub fn convert(&self, a: u64) -> u64 {
        a % self.q
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

#[derive(Clone, Debug)]
/// Polynomial Ring: Rq = Zq[x] / (x^n + 1)
pub struct Poly {
    pub coeffs: Vec<u64>,
    pub ctx: Arc<Context>,
}

// define operations
impl Poly {
    pub fn zero(ctx: &Arc<Context>) -> Self {
        Poly {
            coeffs: vec![0; ctx.degree],
            ctx: ctx.clone(),
        }
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Self) {
        debug_assert!(self.ctx == rhs.ctx);

        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(a, b)| *a = self.ctx.moduli.add_mod(*a, *b));
    }
}

impl Add<&Poly> for &Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Self::Output {
        let mut res = self.clone();
        res += rhs;
        res
    }
}

impl Add for Poly {
    type Output = Poly;

    fn add(self, mut rhs: Poly) -> Self::Output {
        rhs += &self;
        rhs
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        debug_assert!(self.coeffs.len() == rhs.coeffs.len());
        debug_assert!(self.ctx == rhs.ctx);

        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(a, b)| *a = self.ctx.moduli.sub_mod(*a, *b));
    }
}

impl Sub for Poly {
    type Output = Poly;
    fn sub(self, mut rhs: Poly) -> Self::Output {
        rhs -= &self;
        rhs
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, rhs: &Poly) {
        debug_assert!(self.ctx == rhs.ctx);

        let mut res = Poly::zero(&self.ctx);

        for i in 0..res.ctx.degree {
            for j in 0..i + 1 {
                let tmp = res.ctx.moduli.mul_mod(self.coeffs[j], rhs.coeffs[i - j]);
                res.coeffs[i] = res.ctx.moduli.add_mod(res.coeffs[i], tmp);
            }

            for j in i + 1..res.ctx.degree {
                let tmp = res
                    .ctx
                    .moduli
                    .mul_mod(self.coeffs[j], rhs.coeffs[res.ctx.degree + i - j]);
                res.coeffs[i] = res.ctx.moduli.sub_mod(res.coeffs[i], tmp);
            }

            res.coeffs[i] %= res.ctx.moduli.q;
        }

        *self = res;
    }
}

impl Mul<&Poly> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Self::Output {
        let mut res = self.clone();
        res *= rhs;
        res
    }
}

impl Mul<Poly> for Poly {
    type Output = Poly;
    fn mul(self, mut rhs: Poly) -> Self::Output {
        rhs *= &self;
        rhs
    }
}

impl Mul<u64> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: u64) -> Self::Output {
        assert!(rhs < self.ctx.moduli.q);
        let mut res = Poly::zero(&self.ctx);
        res.coeffs
            .iter_mut()
            .for_each(|v| *v = res.ctx.moduli.mul_mod(*v, rhs));
        res
    }
}

impl Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Self::Output {
        Poly {
            coeffs: self
                .coeffs
                .iter()
                .map(|v| self.ctx.moduli.neg(*v))
                .collect(),
            ctx: self.ctx.clone(),
        }
    }
}

impl Neg for Poly {
    type Output = Poly;
    fn neg(self) -> Self::Output {
        Poly {
            coeffs: self
                .coeffs
                .iter()
                .map(|v| self.ctx.moduli.neg(*v))
                .collect(),
            ctx: self.ctx.clone(),
        }
    }
}

impl Poly {
    pub fn random_uniform(ctx: &Arc<Context>) -> Self {
        let mut polynomial = Poly::zero(ctx);
        polynomial
            .coeffs
            .as_mut_slice()
            .copy_from_slice(&sample_uniform_vec(
                polynomial.ctx.degree,
                polynomial.ctx.moduli.q,
            ));
        polynomial
    }

    pub fn random_gaussian(ctx: &Arc<Context>, std_dev: f64) -> Self {
        let mut polynomial = Poly::zero(ctx);
        polynomial
            .coeffs
            .as_mut_slice()
            .copy_from_slice(&Vec::<u64>::from_iter(
                sample_gaussian_vec(polynomial.ctx.degree, std_dev)
                    .into_iter()
                    .map(|mut v| {
                        if v < 0.0 {
                            v += polynomial.ctx.moduli.q as f64;
                        }

                        (v as u64) % polynomial.ctx.moduli.q
                    }),
            ));
        polynomial
    }

    pub fn switch_context(&mut self, ctx: &Arc<Context>) {
        debug_assert!(self.ctx.moduli.q < ctx.moduli.q);
        self.ctx = ctx.clone();
    }

    pub fn decompose(&self, base: u64) -> Vec<Poly> {
        // base should be power of 2
        debug_assert!((base & (base - 1)) == 0);

        let l = ((self.ctx.moduli.q as f64).log2() / (base as f64).log2()).floor() as u64;

        let decomposed_coeffs: Vec<Vec<u64>> = self
            .coeffs
            .iter()
            .map(|v| decompose_value(*v, base.try_into().unwrap()))
            .collect();

        let poly_ctx = Arc::new(Context::new(Modulus { q: base }, self.ctx.degree));

        (0..l as usize)
            .into_iter()
            .map(|d_index| {
                let mut p = Poly::zero(&poly_ctx);
                p.coeffs
                    .iter_mut()
                    .enumerate()
                    .for_each(|(index, c)| *c = decomposed_coeffs[index][d_index]);
                p
            })
            .collect()
    }
}

/// Decomposes `value` into `window_size` bits.
/// For `value` it returns [a0, a1, ..., ak]
/// where each value `a{x}` is `window_size` bit
/// value and `value = a0 * k^0 + a1 * k + a1 * k^2`
/// where `k = 2 ** window_size`
pub fn decompose_value(mut value: u64, window_size: usize) -> Vec<u64> {
    let mut bitvec = Vec::new();
    value = value.to_be();
    while value != 0 {
        bitvec.push(value & 1);
        value >>= 1;
    }
    bitvec
        .chunks(window_size)
        .map(|chunk| chunk.iter().rev().fold(0, |acc, v| (acc << 1) + *v))
        .collect()
}
