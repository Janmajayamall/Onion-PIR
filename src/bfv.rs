use crate::ksk::Ksk;
use crate::rns::ScalingFactor;
use crate::rq::{BitDecomposition, Poly, Representation, RqContext, RqScaler, Substitution};
use crate::utils::div_ceil;
use crate::{
    modulus::Modulus,
    utils::{generate_prime, sample_vec_cbd},
};

use itertools::Itertools;
use ndarray::Axis;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
    sync::Arc,
};

#[derive(PartialEq)]
pub struct BfvParameters {
    pub degree: usize,

    pub plaintext_modulus_u64: u64,
    pub pt_ctx: Arc<RqContext>,

    pub ciphertext_moduli: Vec<u64>,

    pub q_ctxs: Vec<Arc<RqContext>>,
    pub e_ctxs: Vec<Arc<RqContext>>,
    pub mul_params: Vec<MultiplicationParameters>,
    pub delta_ts: Vec<Poly>,
    pub q_mod_ts: Vec<u64>,
    pub scalars: Vec<RqScaler>,

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
        // We need n+1 moduli of 62 bits for multiplication
        let mut extended_basis_moduli = vec![];
        let mut upper_bound = 1 << 62;
        while extended_basis_moduli.len() != (ciphertext_moduli.len() + 1) {
            upper_bound = generate_prime(62, (2 * degree) as u64, upper_bound).unwrap();
            if !extended_basis_moduli.contains(&upper_bound)
                && !ciphertext_moduli.contains(&upper_bound)
            {
                extended_basis_moduli.push(upper_bound);
            }
        }

        let t_invs_rests = ciphertext_moduli
            .iter()
            .map(|q| {
                let q = Modulus::new(*q);
                q.inv(q.neg(plaintext_modulus_u64))
            })
            .collect_vec();

        // Setting parameters for every level afforded
        // by the parent modulus. This is needed so that
        // plaintexts can be encoded at eny level during
        // operations.
        let mut q_ctxs = vec![];
        let mut e_ctxs = vec![];
        let mut mul_params = vec![];
        let mut q_mod_ts = vec![];
        let mut t_inv_polys = vec![];
        let mut scalars = vec![];

        let pt_ctx = Arc::new(RqContext::new(
            vec![ciphertext_moduli[0]],
            degree,
            bit_decomp.clone(),
        ));
        for i in 0..ciphertext_moduli.len() {
            let moduli = ciphertext_moduli[..(ciphertext_moduli.len() - i)].to_vec();
            let q_ctx_i = Arc::new(RqContext::new(moduli.clone(), degree, bit_decomp.clone()));
            let qi_mod_t_i = (q_ctx_i.rns.modulus() % plaintext_modulus_u64)
                .to_u64()
                .unwrap();
            q_ctxs.push(q_ctx_i.clone());
            q_mod_ts.push(qi_mod_t_i);

            // (-t)^-1
            let mut t_inv_poly =
                Poly::try_from_bigint(&q_ctx_i, &[q_ctx_i.rns.lift((&t_invs_rests).into())]);
            t_inv_poly.change_representation(Representation::Ntt);
            t_inv_polys.push(t_inv_poly);

            // scaler
            let scalar_i = RqScaler::new(
                &q_ctx_i,
                &pt_ctx,
                ScalingFactor::new(
                    BigUint::from_u64(plaintext_modulus_u64).unwrap(),
                    q_ctx_i.rns.modulus().clone(),
                ),
            );
            scalars.push(scalar_i);

            // parameters for multiplication
            let moduli_size_sum: usize = moduli
                .iter()
                .map(|q| 64 - (q.leading_zeros()) as usize)
                .sum();
            // extended basis requires: `moduli size` + 60 bits
            let extended_basis_count = div_ceil(moduli_size_sum + 60, 62);
            // this is needed for layer 1 multiplication
            let extended_basis_i = extended_basis_moduli[..extended_basis_count].to_vec();
            let e_ctx_i = Arc::new(RqContext::new(extended_basis_i, degree, bit_decomp.clone()));
            let mul_i = MultiplicationParameters::new(
                &q_ctx_i,
                &e_ctx_i,
                ScalingFactor::new(BigUint::one(), BigUint::one()),
                // Going from PQ to Q, we need to scale down
                // ct by t/Q since after multiplication ct contains
                // delta^2.
                ScalingFactor::new(
                    BigUint::from_u64(plaintext_modulus_u64).unwrap(),
                    q_ctx_i.rns.modulus().clone(),
                ),
            );
            e_ctxs.push(e_ctx_i);
            mul_params.push(mul_i);
        }

        Self {
            degree,
            plaintext_modulus_u64,
            pt_ctx,

            ciphertext_moduli,

            q_ctxs,
            e_ctxs,
            q_mod_ts,
            scalars,
            delta_ts: t_inv_polys,
            mul_params,

            variance,
        }
    }

    pub fn ctx_level(&self, of_ctx: &Arc<RqContext>) -> Option<usize> {
        for i in 0..(self.q_ctxs.len()) {
            if *of_ctx == self.q_ctxs[i] {
                return Some(i);
            }
        }
        None
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

#[derive(PartialEq)]
pub struct MultiplicationParameters {
    from: Arc<RqContext>,
    to: Arc<RqContext>,
    up_scaler: RqScaler,
    down_scalar: RqScaler,
}

impl MultiplicationParameters {
    pub fn new(
        from: &Arc<RqContext>,
        to: &Arc<RqContext>,
        up_factor: ScalingFactor,
        down_factor: ScalingFactor,
    ) -> Self {
        Self {
            from: from.clone(),
            to: to.clone(),
            up_scaler: RqScaler::new(from, to, up_factor),
            down_scalar: RqScaler::new(from, to, down_factor),
        }
    }
}

impl Debug for BfvParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BfvParameters")
            .field("polynomial_degree", &self.degree)
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
            level: 0, // FIXME: change this
        }
    }

    pub fn encrypt(&self, pt: &Plaintext) -> BfvCipherText {
        let poly = pt.to_poly();
        self.encrypt_poly(&poly)
    }

    pub fn decrypt_trial(&self, ct0: &Poly, ct1: &Poly) -> Poly {
        assert_eq!(ct0.context, ct1.context);
        let mut sk = Poly::try_from_vec_i64(&ct0.context, &self.coeffs);
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
        m = self.params.scalars[ct.level].scale(&m);
        let mut m: Vec<u64> = m
            .coefficients()
            .index_axis(Axis(0), 0)
            .iter()
            .map(|a| *a + self.params.plaintext_modulus_u64)
            .collect();
        m = m[..self.params.degree].to_vec();
        let q = self.params.q_ctxs[ct.level].moduli[0].clone();
        m = q.reduce_vec_u64(&m);
        m = Modulus::new(self.params.plaintext_modulus_u64).reduce_vec_u64(&m);

        Plaintext {
            params: self.params.clone(),
            values: m.into_boxed_slice(),
            level: ct.level,
        }
    }

    /// Measures difference between ideal message value in ct space
    /// and the real message value in ct space
    pub fn measure_noise(&self, ct: &BfvCipherText) -> usize {
        let pt = self.decrypt(ct);
        let ideal_m = pt.to_poly();
        let ctx = self.params.q_ctxs[ct.level].clone();

        let mut sk = Poly::try_from_vec_i64(&ctx, &self.coeffs);
        sk.change_representation(Representation::Ntt);
        let mut m = &ct.cts[1] * &sk;
        m += &ct.cts[0];

        m -= &ideal_m;
        m.change_representation(Representation::PowerBasis);

        let mut noise = 0usize;
        let ct_moduli = &ctx.rns.product.clone();

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
    pub level: usize,
}

impl Add<&BfvCipherText> for &BfvCipherText {
    type Output = BfvCipherText;
    fn add(self, rhs: &BfvCipherText) -> Self::Output {
        assert!(self.params == rhs.params);
        assert!(self.level == rhs.level);

        let c0 = &self.cts[0] + &rhs.cts[0];
        let c1 = &self.cts[1] + &rhs.cts[1];

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
            level: self.level,
        }
    }
}

impl Add<&Poly> for &BfvCipherText {
    type Output = BfvCipherText;
    fn add(self, rhs: &Poly) -> Self::Output {
        assert!(self.cts[0].context == rhs.context);

        let c0 = &self.cts[0] + rhs;
        let c1 = &self.cts[1] + rhs;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
            level: self.level,
        }
    }
}

impl Sub<&BfvCipherText> for &BfvCipherText {
    type Output = BfvCipherText;
    fn sub(self, rhs: &BfvCipherText) -> Self::Output {
        assert!(self.params == rhs.params);
        assert!(self.level == rhs.level);

        let c0 = &self.cts[0] - &rhs.cts[0];
        let c1 = &self.cts[1] - &rhs.cts[1];

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
            level: self.level,
        }
    }
}

impl Mul<&Plaintext> for &BfvCipherText {
    type Output = BfvCipherText;
    fn mul(self, rhs: &Plaintext) -> Self::Output {
        assert!(self.level == rhs.level);

        // TODO: store this in BfvPlaintext
        let mut pt_poly = Poly::try_from_vec_u64(&self.params.q_ctxs[self.level], &rhs.values);
        pt_poly.change_representation(Representation::Ntt);

        let c0 = &self.cts[0] * &pt_poly;
        let c1 = &self.cts[1] * &pt_poly;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
            level: self.level,
        }
    }
}

impl Mul<&Poly> for &BfvCipherText {
    type Output = BfvCipherText;
    fn mul(self, rhs: &Poly) -> Self::Output {
        assert!(rhs.context == self.cts[0].context);
        assert!(rhs.representation == Representation::Ntt);

        let c0 = &self.cts[0] * rhs;
        let c1 = &self.cts[1] * rhs;

        BfvCipherText {
            params: self.params.clone(),
            cts: vec![c0, c1],
            level: self.level,
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
            level: self.level,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Plaintext {
    pub params: Arc<BfvParameters>,
    pub values: Box<[u64]>,
    pub level: usize,
}

impl Plaintext {
    pub fn new(params: &Arc<BfvParameters>, values: &Vec<u64>, level: usize) -> Self {
        assert!(values.len() <= params.degree);
        Self {
            params: params.clone(),
            values: values.clone().into_boxed_slice(),
            level,
        }
    }

    pub fn to_poly(&self) -> Poly {
        let q_mod_t = self.params.q_mod_ts[self.level];
        let q_ctx = self.params.q_ctxs[self.level].clone();
        let delta = self.params.delta_ts[self.level].clone();

        // scale m values by q
        let values =
            Modulus::new(self.params.plaintext_modulus_u64).scalar_mul_vec(&self.values, q_mod_t);

        let mut poly = Poly::try_from_vec_u64(&q_ctx, &values);
        poly.change_representation(Representation::Ntt);

        // Multiply poly by delta
        // delta = constant polynomial `(-t)^-1)` in rq.
        // so `delta * poly` results into a poly with
        // coefficients scaled by delta. Thus producing
        // scaled plaintext ((q * t)*m) in rq.
        poly *= &delta;

        poly
    }
}

/// Special key that perform key switching operation
/// from `s^i` to `s`, where `i` is the substitution
/// exponent
pub struct GaliosKey {
    gk: Ksk,
    substitution: Substitution,
    galios_level: usize,
}

impl GaliosKey {
    pub fn new(sk: &SecretKey, i: &Substitution, galios_level: usize) -> Self {
        assert!(sk.params.ciphertext_moduli.len() != 1);

        let ctx = sk.params.q_ctxs[galios_level].clone();

        let mut sk_poly = Poly::try_from_vec_i64(&ctx, &sk.coeffs);
        let mut sk_i_poly = sk_poly.substitute(i);
        sk_i_poly.change_representation(Representation::Ntt);

        let gk = Ksk::new(sk, &sk_i_poly, galios_level, galios_level);

        GaliosKey {
            gk,
            substitution: i.clone(),
            galios_level,
        }
    }

    pub fn relinearize(&self, ct: &BfvCipherText) -> BfvCipherText {
        assert!(ct.level == self.galios_level);

        debug_assert!(ct.cts[0].representation == Representation::Ntt);

        let mut c1_r = ct.cts[1].substitute(&self.substitution);
        c1_r.change_representation(Representation::PowerBasis);
        let (mut c0, c1) = self.gk.key_switch(&c1_r);
        c0 += &ct.cts[0].substitute(&self.substitution);

        BfvCipherText {
            params: ct.params.clone(),
            cts: vec![c0, c1],
            level: self.galios_level,
        }
    }
}

#[cfg(test)]
mod tests {
    use sha2::digest::typenum::Mod;

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
                values: Modulus::new(params.plaintext_modulus_u64)
                    .random_vec(params.degree, &mut rng)
                    .into_boxed_slice(),
                level: 0,
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
            let gk = GaliosKey::new(&sk, &subs, 0);

            let v = Modulus::new(params.plaintext_modulus_u64).random_vec(params.degree, &mut rng);
            let pt = Plaintext::new(&params, &v, 0);
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
