use crate::ksk::Ksk;
use crate::rns::ScalingFactor;
use crate::rq::{BitDecomposition, Poly, Representation, RqContext, RqScaler, Substitution};
use crate::{
    poly::Modulus,
    utils::{generate_prime, sample_vec_cbd},
};
use ndarray::Axis;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use sha2::digest::typenum::bit;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
    sync::Arc,
};

#[derive(PartialEq)]
pub struct BfvParameters {
    pub degree: usize,

    pub plaintext_modulus_u64: u64,
    pub plaintext_modulus: Modulus,
    pub plaintext_context: Arc<RqContext>,

    pub ciphertext_moduli: Vec<u64>,
    pub rq_context: Arc<RqContext>,
    pub scalar: RqScaler,
    pub delta: Poly,
    pub q_mod_t: u64,

    /// Error variance
    pub variance: usize,
}

impl BfvParameters {
    pub fn new(
        degree: usize,
        plaintext_modulus_u64: u64,
        ciphertext_moduli: Vec<u64>,
        variance: usize,
        bit_decomp: BitDecomposition,
    ) -> Self {
        let pt_context = Arc::new(RqContext::new(
            vec![ciphertext_moduli[0]],
            degree,
            bit_decomp.clone(),
        ));
        let rq_context = Arc::new(RqContext::new(
            ciphertext_moduli.clone(),
            degree,
            bit_decomp,
        ));

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
        let mut delta = Poly::try_from_bigint(&rq_context, &[delta]);
        delta.change_representation(Representation::Ntt);

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

    pub fn default(
        num_of_moduli: usize,
        degree: usize,
        bit_decomp: BitDecomposition,
    ) -> BfvParameters {
        let moduli = BfvParameters::generate_moduli(&vec![62usize; num_of_moduli], degree).unwrap();
        BfvParameters::new(degree, 1153u64, moduli, 10, bit_decomp)
    }

    /// Generate ciphertext moduli with the specified sizes
    ///
    /// Ref - https://github.com/tlepoint/fhe.rs/blob/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d/crates/fhe/src/bfv/parameters.rs#L281
    pub fn generate_moduli(moduli_sizes: &[usize], degree: usize) -> Option<Vec<u64>> {
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

impl Debug for BfvParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BfvParameters")
            .field("polynomial_degree", &self.degree)
            .field("plaintext_modulus", &self.plaintext_modulus)
            .field("moduli", &self.ciphertext_moduli)
            // .field("moduli_sizes", &self.moduli_sizes)
            // .field("variance", &self.variance)
            // .field("ctx", &self.ctx)
            // .field("op", &self.op)
            // .field("delta", &self.delta)
            // .field("q_mod_t", &self.q_mod_t)
            // .field("scaler", &self.scaler)
            // .field("plaintext", &self.plaintext)
            // .field("mul_params", &self.mul_params)
            // .field("matrix_reps_index_map", &self.matrix_reps_index_map)
            .finish()
    }
}

#[derive(Debug)]
pub struct SecretKey {
    pub params: Arc<BfvParameters>,
    pub coeffs: Box<[i64]>,
}

impl SecretKey {
    pub fn generate(params: &Arc<BfvParameters>) -> Self {
        let mut rng = thread_rng();
        let coeffs = sample_vec_cbd(params.degree, params.variance, &mut rng).into_boxed_slice();
        // let coeffs = (0..(params.degree)).into_iter().map(|_| 1i64).collect();
        // let mut coeffs: Box<[i64]> = (0..(params.degree)).into_iter().map(|_| 0i64).collect();
        // coeffs[0] = 1;

        Self {
            params: params.clone(),
            coeffs,
        }
    }

    pub fn new(coeffs: Vec<i64>, params: &Arc<BfvParameters>) -> Self {
        Self {
            params: params.clone(),
            coeffs: coeffs.into_boxed_slice(),
        }
    }

    pub fn encrypt_poly(&self, pt: &Poly) -> BfvCipherText {
        let mut sk = Poly::try_from_vec_i64(&pt.context, &self.coeffs);
        sk.change_representation(Representation::Ntt);

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        thread_rng().fill(&mut seed);
        let a = Poly::random_from_seed(&pt.context, Representation::Ntt, seed);

        let mut a_s = a.clone();
        // a * sk
        a_s *= &sk;

        let mut rng = thread_rng();
        // e
        let mut b = Poly::random_small(
            &pt.context,
            Representation::Ntt,
            self.params.variance,
            &mut rng,
        );
        // -(a * sk) + e
        b -= &a_s;
        // b = -(a * sk) + e + pt
        b += pt;

        BfvCipherText {
            cts: vec![b, a],
            params: self.params.clone(),
        }
    }

    pub fn encrypt(&self, pt: &Plaintext) -> BfvCipherText {
        let poly = pt.to_poly();
        self.encrypt_poly(&poly)
    }

    pub fn decrypt_trial(&self, ct0: &Poly, ct1: &Poly) -> Poly {
        let mut sk = Poly::try_from_vec_i64(&self.params.rq_context, &self.coeffs);
        sk.change_representation(Representation::Ntt);

        let mut v = &sk * ct1;
        v += ct0;

        v
    }

    pub fn decrypt(&self, ct: &BfvCipherText) -> Plaintext {
        let mut sk = Poly::try_from_vec_i64(&ct.cts[0].context, &self.coeffs);
        sk.change_representation(Representation::Ntt);

        // c[0] = b = -(a * sk) + e + delta_m
        // c[1] = a
        // e + delta_m = c[0] + c[1]sk
        let mut m = ct.cts[0].clone();
        // sk * c[1]
        sk *= &ct.cts[1];
        m += &sk;
        m.change_representation(Representation::PowerBasis);

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
        m = q.reduce_vec_u64(&m);
        m = self.params.plaintext_modulus.reduce_vec_u64(&m);

        Plaintext {
            params: self.params.clone(),
            values: m.into_boxed_slice(),
        }
    }

    /// Measures difference between ideal message value in ct space
    /// and the real message value in ct space
    pub fn measure_noise(&self, ct: &BfvCipherText) -> usize {
        let pt = self.decrypt(ct);
        let ideal_m = pt.to_poly();

        let mut sk = Poly::try_from_vec_i64(&self.params.rq_context, &self.coeffs);
        sk.change_representation(Representation::Ntt);
        let mut m = &ct.cts[1] * &sk;
        m += &ct.cts[0];

        m -= &ideal_m;
        m.change_representation(Representation::PowerBasis);

        let mut noise = 0usize;
        let ct_moduli = &self.params.rq_context.rns.product.clone();

        Vec::<BigUint>::from(&m).iter().for_each(|coeff| {
            noise = std::cmp::max(
                noise,
                std::cmp::min(coeff.bits(), (ct_moduli - coeff).bits()) as usize,
            );
        });
        noise
    }
}

#[derive(Clone)]
pub struct BfvCiphertext {
    cts: Vec<Poly>,
}

#[derive(Debug, Clone)]
pub struct BfvCipherText {
    pub params: Arc<BfvParameters>,
    pub cts: Vec<Poly>,
}

impl Add<&BfvCipherText> for &BfvCipherText {
    type Output = BfvCipherText;
    fn add(self, rhs: &BfvCipherText) -> Self::Output {
        assert!(self.params == rhs.params);

        let c0 = &self.cts[0] + &rhs.cts[0];
        let c1 = &self.cts[1] + &rhs.cts[1];

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

impl Add<&Poly> for &BfvCipherText {
    type Output = BfvCipherText;
    fn add(self, rhs: &Poly) -> Self::Output {
        assert!(self.params.rq_context == rhs.context);

        let c0 = &self.cts[0] + rhs;
        let c1 = &self.cts[1] + rhs;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

impl Sub<&BfvCipherText> for &BfvCipherText {
    type Output = BfvCipherText;
    fn sub(self, rhs: &BfvCipherText) -> Self::Output {
        assert!(self.params == rhs.params);

        let c0 = &self.cts[0] - &rhs.cts[0];
        let c1 = &self.cts[1] - &rhs.cts[1];

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

impl Mul<&Plaintext> for &BfvCipherText {
    type Output = BfvCipherText;
    fn mul(self, rhs: &Plaintext) -> Self::Output {
        assert!(rhs.params == self.params);

        // TODO: store this in BfvPlaintext
        let mut pt_poly = Poly::try_from_vec_u64(&self.params.rq_context, &rhs.values);
        pt_poly.change_representation(Representation::Ntt);

        let c0 = &self.cts[0] * &pt_poly;
        let c1 = &self.cts[1] * &pt_poly;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

impl Mul<&Poly> for &BfvCipherText {
    type Output = BfvCipherText;
    fn mul(self, rhs: &Poly) -> Self::Output {
        assert!(rhs.context == self.params.rq_context);
        assert!(rhs.representation == Representation::Ntt);

        let c0 = &self.cts[0] * rhs;
        let c1 = &self.cts[1] * rhs;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

impl Mul<&BigUint> for &BfvCipherText {
    type Output = BfvCipherText;
    fn mul(self, rhs: &BigUint) -> Self::Output {
        let c0 = &self.cts[0] * rhs;
        let c1 = &self.cts[1] * rhs;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Plaintext {
    pub params: Arc<BfvParameters>,
    pub values: Box<[u64]>,
}

impl Plaintext {
    pub fn new(params: &Arc<BfvParameters>, values: &Vec<u64>) -> Self {
        assert!(values.len() <= params.degree);
        Self {
            params: params.clone(),
            values: values.clone().into_boxed_slice(),
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

/// Special key that perform key switching operation
/// from `s^i` to `s`, where `i` is the substitution
/// exponent
pub struct GaliosKey {
    gk: Ksk,
    substitution: Substitution,
}

impl GaliosKey {
    pub fn new(sk: &SecretKey, i: &Substitution) -> Self {
        // TODO: Q. why is this the case ?
        assert!(sk.params.ciphertext_moduli.len() != 1);

        let mut sk_poly = Poly::try_from_vec_i64(&sk.params.rq_context, &sk.coeffs);
        let mut sk_i_poly = sk_poly.substitute(i);
        sk_i_poly.change_representation(Representation::Ntt);

        let gk = Ksk::new(sk, &sk_i_poly);

        GaliosKey {
            gk,
            substitution: i.clone(),
        }
    }

    pub fn relinearize(&self, ct: &BfvCipherText) -> BfvCipherText {
        assert!(self.gk.params.rq_context == ct.cts[0].context);

        debug_assert!(ct.cts[0].representation == Representation::Ntt);

        let mut c1_r = ct.cts[1].substitute(&self.substitution);
        c1_r.change_representation(Representation::PowerBasis);
        let (mut c0, c1) = self.gk.key_switch(&c1_r);
        c0 += &ct.cts[0].substitute(&self.substitution);

        BfvCipherText {
            params: ct.params.clone(),
            cts: vec![c0, c1],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(
            6,
            8,
            BitDecomposition { base: 4, l: 8 },
        ));

        for _ in 0..100 {
            let sk = SecretKey::generate(&params);
            let pt = Plaintext {
                params: params.clone(),
                values: params
                    .plaintext_modulus
                    .random_vec(params.degree, &mut rng)
                    .into_boxed_slice(),
            };
            let ct = sk.encrypt(&pt);
            let pt_after = sk.decrypt(&ct);

            assert_eq!(pt.values, pt_after.values);
        }
    }

    #[test]
    fn galois_key() {
        for _ in 0..100 {
            let mut rng = thread_rng();
            let params = Arc::new(BfvParameters::default(
                2,
                8,
                BitDecomposition { base: 4, l: 8 },
            ));
            let sk = SecretKey::generate(&params);

            let subs = Substitution::new(3);
            let gk = GaliosKey::new(&sk, &subs);

            let v = params.plaintext_modulus.random_vec(params.degree, &mut rng);
            let pt = Plaintext::new(&params, &v);
            let ct = sk.encrypt(&pt);

            let ct2 = gk.relinearize(&ct);
            let vr = sk.decrypt(&ct2);

            let pt_ctx = Arc::new(RqContext::new(
                vec![params.plaintext_modulus_u64],
                params.degree,
                BitDecomposition { base: 4, l: 8 },
            ));
            let expected_vr = Poly::try_from_vec_u64(&pt_ctx, &v);

            assert_eq!(
                vr.values.to_vec(),
                Vec::<u64>::from(&expected_vr.substitute(&subs))
            );
        }
    }
}
