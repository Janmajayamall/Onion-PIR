use crate::poly::Modulus;
use num_traits::cast::ToPrimitive;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

#[derive(Debug)]
struct NttOperator {
    p: Modulus,
    p_twice: u64,
    omegas: Vec<u64>,
    omegas_inv: Vec<u64>,
    omegas_shoup: Vec<u64>,
    zetas_inv: Vec<u64>,
    zetas_inv_shoup: Vec<u64>,
    size: usize,
    size_inv: u64,
    size_inv_shoup: u64,
}

impl NttOperator {
    pub fn new(p: &Modulus, size: usize) -> Self {
        // check whether params are valid for NttOperator
        assert!(supports_ntt(p.p, size));

        let omega = primitive_root(2 * size, p);
        let omega_inv = p.inv(omega);

        let mut powers = Vec::with_capacity(size + 1);
        let mut powers_inv = Vec::with_capacity(size + 1);
        let mut v = 1u64;
        let mut v_inv = 1u64;
        (0..(size + 1)).into_iter().for_each(|_| {
            powers.push(v);
            powers_inv.push(v_inv);
            v = p.mul(v, omega);
            v_inv = p.mul(v_inv, omega_inv);
        });

        let mut omegas = Vec::with_capacity(size);
        let mut omegas_inv = Vec::with_capacity(size);
        let mut zetas_inv = Vec::with_capacity(size);
        (0..size).into_iter().for_each(|i| {
            // this arranges powers of omegas in desired
            // sequence as we access them in ntt transformation
            let j = i.reverse_bits() >> (size.leading_zeros() + 1) as usize;
            omegas.push(powers[j]);
            omegas_inv.push(powers_inv[j]);
            zetas_inv.push(powers_inv[j + 1]);
        });

        let omegas_shoup = p.shoup_vec(&omegas);
        let zetas_inv_shoup = p.shoup_vec(&zetas_inv);

        let size_inv = p.inv(size as u64);
        let size_inv_shoup = p.shoup(size_inv);

        Self {
            p: p.clone(),
            p_twice: p.p * 2,
            omegas,
            omegas_inv,
            omegas_shoup,
            size,
            size_inv,
            size_inv_shoup,
            zetas_inv,
            zetas_inv_shoup,
        }
    }

    /// Forward ntt transformation
    ///
    /// Ref -
    /// (1) https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/ntt.rs#L89
    /// (2) Algorithm 2 - https://arxiv.org/pdf/2103.16400.pdf
    ///
    pub fn forward(&self, a: &mut [u64]) {
        assert!(a.len() == self.size);

        let n = self.size;
        let a_ptr = a.as_mut_ptr();
        let mut l = n >> 1;
        let mut m = 1;
        let mut k = 1;

        while l > 0 {
            for i in 0..m {
                unsafe {
                    let omega = *self.omegas.get_unchecked(k);
                    let omega_shoup = *self.omegas_shoup.get_unchecked(k);
                    k += 1;

                    let s = 2 * i * l;
                    match l {
                        1 => {
                            // last layer reduces residues from [0, 4p) to [0, p)
                            let uj = &mut *a_ptr.add(s);
                            let ujl = &mut *a_ptr.add(s + 1);
                            self.butterfly(uj, ujl, omega, omega_shoup);
                            *uj = reduce_4p(*uj, &self.p);
                            *ujl = reduce_4p(*ujl, &self.p);
                        }
                        _ => {
                            for j in s..(s + l) {
                                self.butterfly(
                                    &mut *a_ptr.add(j),
                                    &mut *a_ptr.add(j + l),
                                    omega,
                                    omega_shoup,
                                );
                            }
                        }
                    }
                }
            }

            l >>= 1;
            m <<= 1;
        }
    }

    fn backward(&self, a: &mut [u64]) {
        assert!(a.len() == self.size);

        let a_ptr = a.as_mut_ptr();

        let mut k = 0;
        let mut m = self.size >> 1;
        let mut l = 1;

        while m > 0 {
            for i in 0..m {
                let s = 2 * i * l;
                unsafe {
                    let zeta_inv = *self.zetas_inv.get_unchecked(k);
                    let zeta_inv_shoup = *self.zetas_inv_shoup.get_unchecked(k);
                    k += 1;

                    for j in s..(s + l) {
                        self.butterfly_inv(
                            &mut *a_ptr.add(j),
                            &mut *a_ptr.add(j + l),
                            zeta_inv,
                            zeta_inv_shoup,
                        )
                    }
                    // match l {
                    //     1 => self.butterfly_inv(
                    //         &mut *a_ptr.add(s),
                    //         &mut *a_ptr.add(s + l),
                    //         zeta_inv,
                    //         zeta_inv_shoup,
                    //     ),
                    //     _ => {
                    //         for j in s..(s + l) {
                    //             self.butterfly_inv(
                    //                 &mut *a_ptr.add(j),
                    //                 &mut *a_ptr.add(j + l),
                    //                 zeta_inv,
                    //                 zeta_inv_shoup,
                    //             )
                    //         }
                    //     }
                    // }
                }
            }
            l <<= 1;
            m >>= 1;
        }

        a.iter_mut().for_each(|ainv| {
            *ainv = self
                .p
                .lazy_mul_shoup(*ainv, self.size_inv, self.size_inv_shoup);
        });
    }

    /// Harvey NTT butterfly
    ///
    /// Ref -
    /// (1) https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/ntt.rs#L309
    /// (2) Algorithm 4 of https://arxiv.org/pdf/2103.16400.pdf
    fn butterfly(&self, x: &mut u64, y: &mut u64, w: u64, w_shoup: u64) {
        debug_assert!(*x < self.p.p * 4);
        debug_assert!(*y < self.p.p * 4);
        debug_assert!(w < self.p.p);
        debug_assert!(w_shoup == self.p.shoup(w));

        *x = Modulus::reduce_ct(*x, self.p_twice);
        let t = self.p.lazy_mul_shoup(*y, w, w_shoup);
        *y = *x + self.p_twice - t;
        *x += t;

        debug_assert!(*x < self.p.p * 4);
        debug_assert!(*y < self.p.p * 4);
    }

    fn butterfly_inv(&self, x: &mut u64, y: &mut u64, z_inv: u64, z_inv_shoup: u64) {
        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
        debug_assert!(z_inv < self.p.p);
        debug_assert!(z_inv_shoup == self.p.shoup(z_inv));

        let t = *x;
        *x = Modulus::reduce_ct(*y + t, self.p_twice);
        *y = self
            .p
            .lazy_mul_shoup(self.p_twice + t - *y, z_inv, z_inv_shoup);

        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
    }
}

/// Checks the following for ntt support
/// 1. p % 2n === 1 (to make sure that 2nth roots of unit exist and can be precalculated (https://eprint.iacr.org/2017/727))
/// 2. `p` is prime
/// 3. n is power of 2 and greater than 8
///
/// TODO: Add is_prim check
fn supports_ntt(p: u64, n: usize) -> bool {
    if !n.is_power_of_two() || n <= 8 {
        return false;
    }

    (p % (n << 1).to_u64().unwrap()) == 1
}

/// Finds and returns 2n-th primitive root modulo p
///
/// Ref - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/ntt.rs#L373
fn primitive_root(n: usize, p: &Modulus) -> u64 {
    let lambda = (p.p - 1) / (2 * n as u64);

    let mut rng: ChaCha8Rng = SeedableRng::seed_from_u64(0);

    for _ in 0..100 {
        let mut root = rng.gen_range(0..p.p);
        root = p.pow(root, lambda);

        if is_primitive_root(root, 2 * n, p) {
            return root;
        }
    }

    assert!(false);
    0
}

///
/// A primitive root of unity is such that x^n = 1 mod p, and x^(n/p) != 1 mod p
/// for all prime p dividing n.
///
/// Ref - https://github.com/Janmajayamall/fhe.rs/blob/8aafe4396d0b771e6aa25257c7daa61c109eb367/crates/fhe-math/src/zq/ntt.rs#L394
fn is_primitive_root(v: u64, n: usize, p: &Modulus) -> bool {
    debug_assert!(v < p.p);
    p.pow(v, n as u64) == 1 && p.pow(v, (n / 2) as u64) != 1
}

/// Reduces `a` modulo p
///
/// a < p * 4
fn reduce_4p(a: u64, p: &Modulus) -> u64 {
    debug_assert!(a < p.p * 4);
    let r = Modulus::reduce_ct(a, 2 * p.p);
    Modulus::reduce_ct(r, p.p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn construct_ntt_operator() {
        for n in [8, 1024] {
            for q in [1153, 4611686018326724609] {
                let moduli = Modulus::new(q);
                let _ = NttOperator::new(&moduli, n);
            }
        }
    }

    #[test]
    fn bijection() {
        let mut rng = thread_rng();
        for n in [1024] {
            for q in [4611686018326724609] {
                let moduli = Modulus::new(q);
                let operator = NttOperator::new(&moduli, n);

                let mut poly = moduli.random_vec(n, &mut rng);
                let poly_clone = poly.clone();
                operator.forward(&mut poly);
                operator.backward(&mut poly);
                assert_eq!(poly[0], poly_clone[0]);
            }
        }
    }
}
