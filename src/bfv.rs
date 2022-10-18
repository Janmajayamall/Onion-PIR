use itertools::ProcessResults;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::rq::{Poly, Representation, RqContext};
use crate::{
    poly::{Context, Modulus},
    utils::sample_vec_cbd,
};
use std::sync::Arc;

#[derive(Debug, PartialEq)]
pub struct BfvParameters {
    polynomial_degree: usize,

    plaintext_modulus_64: u64,
    plaintext_modulus: Modulus,

    ciphertext_moduli: Vec<u64>,
    rq_context: Arc<RqContext>,

    /// Error variance
    variance: isize,
}

#[derive(Debug)]
pub struct SecretKey {
    pub params: Arc<BfvParameters>,
    coeffs: Box<[i64]>,
}

impl SecretKey {
    pub fn new(params: &Arc<BfvParameters>) -> Self {
        let coeffs = sample_vec_cbd(params.polynomial_degree, params.variance).into_boxed_slice();

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
        // TODO: Scale plaintext poly by q / t
        // return encrypt_poly
        todo!()
    }

    pub fn decrypt(&self, bfvCt: &BfvCipherText) {
        let mut sk = Poly::try_from_vec_i64(&bfvCt.cts[0].context, &self.coeffs);

        // c[0] = b = -a * sk + e + delta_m
        // c[1] = a
        // e + delta_m = c[0] + c[1]sk
        let mut m = bfvCt.cts[0].clone();
        // sk * c[1]
        sk *= &bfvCt.cts[1];
        m += &sk;

        // descale
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

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn encrypt_decrypt() {
    //     let params = Arc::new(BfvParameters::default());
    //     let sk = BfvSecretKey::new(&params);
    //     let pk = BfvPublicKey::new(&sk);
    //     let pt = BfvPlaintext::new(&params, vec![3, 3, 3, 2]);
    //     let ct = pk.encrypt(&pt.poly.coeffs.into());
    //     let pt2 = sk.decrypt(&ct);
    //     assert_eq!(pt.values.into_vec(), pt2);
    // }
}
