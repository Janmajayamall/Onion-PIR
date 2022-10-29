use itertools::izip;
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use rand::{distributions::Uniform, CryptoRng, Rng, RngCore};
#[derive(Clone, PartialEq, Debug)]
pub struct Modulus {
    pub p: u64,
    barret_hi: u64,
    barret_lo: u64,
}

/// Implements Modulus
///
/// TODO:
/// (1) Implement optimized lazy reduction (https://github.com/tlepoint/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L653-L654)
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
    /// (1) https://github.com/tlepoint/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L581-L582
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

    /// Ref - https://github.com/tlepoint/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/mod.rs#L426
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

    pub fn inv_vec(&self, a: &[u64]) -> Vec<u64> {
        a.iter().map(|ai| self.inv(*ai)).collect()
    }

    pub fn scalar_mul_vec(&self, a: &[u64], b: u64) -> Vec<u64> {
        let b_shoup = self.shoup(b);
        a.iter()
            .map(|ai| self.lazy_mul_shoup(*ai, b, b_shoup))
            .collect()
    }

    pub fn add(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p);
        debug_assert!(b < self.p);
        Self::reduce_ct(a + b, self.p)
    }

    pub fn sub(&self, a: u64, b: u64) -> u64 {
        // dbg!(a, b, self.p);
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

    pub fn neg_vec(&self, a: &[u64]) -> Vec<u64> {
        a.iter().map(|ai| self.neg(*ai)).collect()
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

    pub fn mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        Self::reduce_ct(self.lazy_mul_shoup(a, b, b_shoup), self.p)
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

    pub fn bits(&self) -> u64 {
        u64::MAX - self.p.leading_zeros().to_u64().unwrap()
    }
}
