use ndarray::Axis;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rand_distr::Open01;

use crate::rns::ScalingFactor;
use crate::rq::{Poly, Representation, RqContext, RqScaler};
use crate::{
    poly::{Context, Modulus},
    utils::{generate_prime, sample_vec_cbd},
};
use std::sync::Arc;

#[derive(Debug, PartialEq)]
pub struct BfvParameters {
    degree: usize,

    plaintext_modulus_u64: u64,
    plaintext_modulus: Modulus,
    plaintext_context: Arc<RqContext>,

    ciphertext_moduli: Vec<u64>,
    rq_context: Arc<RqContext>,
    scalar: RqScaler,
    delta: Poly,
    q_mod_t: u64,

    /// Error variance
    variance: isize,
}

impl BfvParameters {
    pub fn new(
        degree: usize,
        plaintext_modulus_u64: u64,
        ciphertext_moduli: Vec<u64>,
        variance: isize,
    ) -> Self {
        let pt_context = Arc::new(RqContext::new(vec![ciphertext_moduli[0]], degree));
        let rq_context = Arc::new(RqContext::new(ciphertext_moduli.clone(), degree));

        // scaler for scaling down
        // delta_m by (t/q) (i.e. from ct space to pt space)
        let scalar = RqScaler::new(
            &rq_context,
            &pt_context,
            ScalingFactor::new(
                BigUint::from_u64(plaintext_modulus_u64).unwrap(),
                rq_context.rns.product.clone(),
            ),
        );

        // delta = [-t^(-1)]_q
        let delta_rests: Vec<u64> = ciphertext_moduli
            .iter()
            .map(|qi| {
                let qi = Modulus::new(*qi);
                qi.inv(qi.neg(plaintext_modulus_u64))
            })
            .collect();
        let delta = rq_context.rns.lift((&delta_rests).into());
        let delta = Poly::try_from_bigint(&rq_context, &vec![delta]);

        let q_mod_t = (rq_context.rns.product.clone() % plaintext_modulus_u64)
            .to_u64()
            .unwrap();

        Self {
            degree,
            plaintext_modulus_u64,
            plaintext_modulus: Modulus::new(plaintext_modulus_u64),
            plaintext_context: pt_context,
            ciphertext_moduli,
            rq_context,
            scalar,
            variance,
            delta,
            q_mod_t,
        }
    }

    pub fn default(num_of_moduli: usize, degree: usize) -> BfvParameters {
        let moduli = BfvParameters::generate_moduli(&vec![62usize; num_of_moduli], degree).unwrap();
        BfvParameters::new(degree, 1153u64, moduli, 10)
    }

    /// Generate ciphertext moduli with the specified sizes
    ///
    /// Ref - https://github.com/tlepoint/fhe.rs/blob/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d/crates/fhe/src/bfv/parameters.rs#L281
    fn generate_moduli(moduli_sizes: &[usize], degree: usize) -> Option<Vec<u64>> {
        let mut moduli = vec![];
        for size in moduli_sizes {
            assert!(*size <= 62 && *size >= 10);

            let mut upper_bound = 1 << size;
            loop {
                if let Some(prime) = generate_prime(*size, 2 * degree as u64, upper_bound) {
                    if !moduli.contains(&prime) {
                        moduli.push(prime);
                        break;
                    } else {
                        upper_bound = prime;
                    }
                } else {
                    return None;
                }
            }
        }

        Some(moduli)
    }
}

#[derive(Debug)]
pub struct SecretKey {
    pub params: Arc<BfvParameters>,
    coeffs: Box<[i64]>,
}

impl SecretKey {
    pub fn new(params: &Arc<BfvParameters>) -> Self {
        let coeffs = sample_vec_cbd(params.degree, params.variance).into_boxed_slice();

        Self {
            params: params.clone(),
            coeffs,
        }
    }

    pub fn encrypt_poly(&self, pt: &Poly) -> BfvCipherText {
        let mut sk = Poly::try_from_vec_i64(&pt.context, &self.coeffs);
        sk.change_representation(Representation::Ntt);

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        thread_rng().fill(&mut seed);
        let mut a = Poly::random_from_seed(&pt.context, Representation::Ntt, seed);
        // a * sk
        a *= &sk;

        let mut b = Poly::random_small(&pt.context, Representation::Ntt, self.params.variance);
        // a * sk - e
        b -= &a;
        // b = a * sk - e + pt
        b += pt;

        BfvCipherText { cts: vec![b, a] }
    }

    pub fn encrypt(&self, pt: &BfvPlaintext) -> BfvCipherText {
        let poly = pt.to_poly();
        self.encrypt_poly(&poly)
    }

    pub fn decrypt(&self, bfvCt: &BfvCipherText) -> BfvPlaintext {
        let mut sk = Poly::try_from_vec_i64(&bfvCt.cts[0].context, &self.coeffs);

        // c[0] = b = -a * sk + e + delta_m
        // c[1] = a
        // e + delta_m = c[0] + c[1]sk
        let mut m = bfvCt.cts[0].clone();
        // sk * c[1]
        sk *= &bfvCt.cts[1];
        m += &sk;

        // We scale `m` by (t/q) scaling factor and switch
        // its rq context from `R_q` to `R_q{i}`. `R_q{i}`
        // is any of ct moduli.
        // i.e. [(t/q)* m ]_r where `r` is one of the `qi`s.
        m = self.params.scalar.scale(&m);
        let mut m: Vec<u64> = m
            .coefficients()
            .index_axis(Axis(0), 0)
            .iter()
            .map(|a| *a + self.params.plaintext_modulus.p)
            .collect();
        m = m[..self.params.degree].to_vec();
        let q = Modulus::new(self.params.ciphertext_moduli[0]);
        q.reduce_vec_u64(&mut m);
        self.params.plaintext_modulus.reduce_vec_u64(&mut m);

        BfvPlaintext {
            params: self.params.clone(),
            values: m.into_boxed_slice(),
        }
    }
}

#[derive(Clone)]
pub struct BfvCiphertext {
    cts: Vec<Poly>,
}

#[derive(Debug, Clone)]
pub struct BfvCipherText {
    pub cts: Vec<Poly>,
}

pub struct BfvPlaintext {
    params: Arc<BfvParameters>,
    values: Box<[u64]>,
}

impl BfvPlaintext {
    pub fn new(params: Arc<BfvParameters>, values: Vec<u64>) -> Self {
        assert!(values.len() <= params.degree);
        Self {
            params,
            values: values.into_boxed_slice(),
        }
    }

    pub fn to_poly(&self) -> Poly {
        // scale m values by q
        let values = self
            .params
            .plaintext_modulus
            .scalar_mul_vec(&self.values, self.params.q_mod_t);

        let mut poly = Poly::try_from_vec_u64(&self.params.rq_context, &values);
        poly.change_representation(Representation::Ntt);

        // Multiply poly by delta
        // delta = constant polynomial `(-t)^-1)` in rq.
        // so `delta * poly` results into a poly with
        // coefficients scaled by delta. Thus producing
        // scaled plaintext ((q * t)*m) in rq.
        poly *= &self.params.delta;

        poly
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let params = BfvParameters::default(6, 8);
    }
}

// impl BfvCipherText {
//     pub fn add_ciphertexts(ct0: &BfvCipherText, ct1: &BfvCipherText) -> Self {
//         debug_assert!(ct0.params == ct1.params);
//         let c0 = &ct0.c[0] + &ct1.c[0];
//         let c1 = &ct0.c[1] + &ct1.c[1];
//         BfvCipherText::new(&ct0.params, vec![c0, c1])
//     }

//     pub fn sub_ciphertexts(ct0: &BfvCipherText, ct1: &BfvCipherText) -> Self {
//         debug_assert!(ct0.params == ct1.params);
//         let ct1 = BfvCipherText::negate(ct1);
//         BfvCipherText::add_ciphertexts(ct0, &ct1)
//     }

//     pub fn negate(ct: &BfvCipherText) -> BfvCipherText {
//         // to negate: -1 * Ct
//         // -1 mod t == t-1 mod t
//         // TODO: Handle such ops in a better way
//         BfvCipherText::multiply_constant(ct, ct.params.t - 1)
//     }

//     pub fn multiply_pt_poly(ct: &BfvCipherText, pt_poly: &Poly) -> BfvCipherText {
//         // plaintext polynomial must be Rt
//         assert!(pt_poly.ctx.moduli.q <= ct.params.t);
//         let mut pt_poly = pt_poly.clone();
//         pt_poly.switch_context(&ct.params.poly_ctx);

//         let c0 = &ct.c[0] * &pt_poly;
//         let c1 = &ct.c[1] * &pt_poly;
//         BfvCipherText::new(&ct.params, vec![c0, c1])
//     }

//     pub fn multiply_constant(ct: &BfvCipherText, value: u64) -> BfvCipherText {
//         // constant must be in Z_t
//         assert!(value < ct.params.t);

//         let c0 = &ct.c[0] * value;
//         let c1 = &ct.c[1] * value;
//         BfvCipherText::new(&ct.params, vec![c0, c1])
//     }
// }
