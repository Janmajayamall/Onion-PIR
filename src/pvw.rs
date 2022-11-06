use std::vec;

use crate::{modulus::Modulus, utils::sample_vec_cbd};
use itertools::{izip, Itertools};
use ndarray::{Array2, Axis};
use rand::{thread_rng, Rng};
use rand_distr::Uniform;
use sha2::digest::typenum::Mod;

#[derive(Debug)]
pub struct PvwParams {
    pub n: usize,
    pub ell: usize,
    pub m: usize,
    pub q: u64,
    pub variance: usize,
}

impl Default for PvwParams {
    fn default() -> Self {
        Self {
            n: 450,
            ell: 4,
            m: 16000,
            q: 65537,
            variance: 2,
        }
    }
}

impl PvwParams {
    pub fn new() -> PvwParams {
        todo!()
    }
}

pub struct PvwSk(Vec<Vec<u64>>);

impl PvwSk {
    pub fn gen_sk(params: &PvwParams) -> PvwSk {
        let mut rng = thread_rng();

        // Sample from uniform distribution that is mod `q`
        // and return 2d matrix
        let mut sk = vec![];
        let q_mod = Modulus::new(params.q);
        for _ in 0..params.ell {
            sk.push(q_mod.random_vec(params.n, &mut rng));
        }

        PvwSk(sk)
    }

    /// Encrypts message vector `m`
    /// and returns Ciphertext
    pub fn gen_pk(&self, params: &PvwParams) -> PvwPk {
        let mut rng = thread_rng();
        let q_mod = Modulus::new(params.q);
        let sk = self.0.clone();

        let mut a = vec![];
        for _ in 0..params.m {
            a.push(q_mod.random_vec(params.n, &mut rng));
        }

        let mut p = vec![vec![0; params.ell]; params.m];
        for l_index in 0..params.ell {
            let err = q_mod
                .reduce_vec_i64(sample_vec_cbd(params.m, params.variance, &mut rng).as_slice());
            for m_index in 0..params.m {
                let mut sum = 0;
                for n_index in 0..params.n {
                    let product = q_mod.mul(sk[l_index][n_index], a[m_index][n_index]);
                    sum = q_mod.add(sum, product);
                }
                // add error value
                sum = q_mod.add(sum, err[m_index]);
                p[m_index][l_index] = sum;
            }
        }

        // P = sk * A + e
        let pk = (0..params.m)
            .into_iter()
            .map(|index| PvwCiphertext {
                b: p[index].clone(),
                a: a[index].clone(),
            })
            .collect();

        PvwPk(pk)
    }

    pub fn decrypt(&self, params: &PvwParams, ct: &PvwCiphertext) -> Vec<u64> {
        let sk = self.0.clone();
        let q_mod = Modulus::new(params.q);

        let mut b = ct.b.clone();

        // b - sk * a
        for l_index in 0..params.ell {
            let mut sum = 0;
            for n_index in 0..params.n {
                let product = q_mod.mul(sk[l_index][n_index], ct.a[n_index]);
                sum = q_mod.add(sum, product);
            }
            b[l_index] = q_mod.sub(b[l_index], sum);
        }

        b.iter_mut().for_each(|v| {
            *v = ((*v + (params.q / 4)) / (params.q / 2)) % 2;
        });
        b
    }
}

pub struct PvwPk(Vec<PvwCiphertext>);

impl PvwPk {
    pub fn encrypt(&self, params: &PvwParams, m: Vec<u64>) -> PvwCiphertext {
        let pk = self.0.clone();
        let q_mod = Modulus::new(params.q);
        let rng = thread_rng();

        let dist = Uniform::from(0..2u64);
        let e = rng.sample_iter(dist).take(params.m).collect_vec();
        let mut b = vec![0u64; params.ell];
        let mut a = vec![0u64; params.n];

        // Ae, Pe + t
        for m_index in 0..params.m {
            for l_index in 0..params.ell {
                let product = q_mod.mul(pk[m_index].b[l_index], e[m_index]);
                b[l_index] = q_mod.add(b[l_index], product);
            }

            for n_index in 0..params.n {
                let product = q_mod.mul(pk[m_index].a[n_index], e[m_index]);
                a[n_index] = q_mod.add(a[n_index], product);
            }
        }

        izip!(b.iter_mut(), m.iter()).for_each(|(v, m_v)| {
            *v = q_mod.add(*v, (params.q / 2) * m_v);
        });

        PvwCiphertext { b, a }
    }
}

#[derive(Clone, Debug)]
pub struct PvwCiphertext {
    /// a \in (Z_q)^n
    pub a: Vec<u64>,
    /// b \in (Z_q)^ell
    pub b: Vec<u64>,
}

mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let params = PvwParams::default();

        for _ in 0..5 {
            let sk = PvwSk::gen_sk(&params);
            let pk = sk.gen_pk(&params);
            let mut rng = thread_rng();
            let distr = Uniform::from(0..2);
            let m = rng.sample_iter(distr).take(params.ell).collect_vec();
            let ct = pk.encrypt(&params, m.clone());
            let d_m = sk.decrypt(&params, &ct);

            assert_eq!(d_m, m);
        }
    }
}
