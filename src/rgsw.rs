use std::{sync::Arc, task::Poll};

use super::bfv::BfvCipherText;
use crate::{
    bfv::{BfvParameters, BfvPlaintext, SecretKey},
    ksk::Ksk,
    rq::{Poly, Representation},
};

/// Key switching key
///
/// Params ß & l
/// Note: Least significant first
/// [
///     RLWE_s'( q/ß^l * s)
///     RLWE_s'( q/ß^{l-1} * s)
///     .
///     .   
///     .
///     RLWE_s'( q/ß^1 * s)
/// ]
/// `RLWE_s'` means encrypted under `s'`,
/// where `s'` is the new key
// pub struct Ksk {
//     pub cts: Vec<BfvCipherText>,
//     pub new_pk: BfvPublicKey,
//     pub beta: u64,
//     pub l: u64,
//     pub is_subs: bool,
//     pub subs_k: Option<u64>,
// }

// impl Ksk {
//     pub fn new(
//         curr_sk: &BfvSecretKey,
//         new_sk: &BfvSecretKey,
//         beta: u64,
//         is_subs: bool,
//         subs_k: Option<u64>,
//     ) -> Self {
//         assert!(is_subs && subs_k.is_some() || subs_k.is_none());

//         let q_ref = curr_sk.params.poly_ctx.moduli.q;
//         let l = ((q_ref as f64).log2() / (beta as f64).log2()) as u64;

//         let new_pk = BfvPublicKey::new(new_sk);
//         let cts = (l..0).into_iter().map(|index| {
//             // encrypt under new_sk
//             let pt = &curr_sk.poly * (q_ref / beta.pow(index as u32));
//             let pt = BfvPlaintext::new(&new_sk.params, pt.coeffs.into());
//             new_pk.encrypt(&pt.poly.coeffs)
//         });

//         Self {
//             cts: cts.collect(),
//             new_pk,
//             beta,
//             l,
//             is_subs,
//             subs_k,
//         }
//     }

//     pub fn key_switch(ksk: &Ksk, ct: BfvCipherText) -> BfvCipherText {
//         let a_decomposed = ct.c[1].decompose(ksk.beta, ksk.l);

//         debug_assert!(a_decomposed.len() == ksk.cts.len());

//         // RLWE encryption of `a * curr_sk` under new_sk
//         let mut a_curr_s = a_decomposed
//             .iter()
//             .zip(ksk.cts.iter())
//             .map(|(a_pt, rlwe_ct)| BfvCipherText::multiply_pt_poly(rlwe_ct, a_pt))
//             .fold(
//                 BfvCipherText::new(
//                     &ksk.new_pk.params,
//                     vec![
//                         Poly::zero(&ksk.new_pk.params.poly_ctx),
//                         Poly::zero(&ksk.new_pk.params.poly_ctx),
//                     ],
//                 ),
//                 |acc, ct| BfvCipherText::add_ciphertexts(&acc, &ct),
//             );

//         // RLWE encryption of `b` under new_sk
//         let b_under_new_sk = ksk.new_pk.encrypt(&ct.c[0].coeffs);

//         // (/delta * M) + E_old = B - (A * S_curr)
//         // Thus homomorphic computation of `B - (A * S_curr)` under new_sk
//         // returns `(/delta * M) + E_old` under new_sk.
//         BfvCipherText::sub_ciphertexts(&b_under_new_sk, &a_curr_s)
//     }

//     /// Subs(•,k) converts `RLWE(Σ b_i • X^i)`
//     /// to `RLWE(Σ b_i • (X^i)^k)`.
//     /// The operation uses key switching
//     /// along with the fact that:
//     /// If `ct(X) = (c0(X), c1(X))` for `M(X)`, where
//     /// `M(X) = Σ b_i • X^i` under S(X) as secret key
//     /// then
//     /// `ct(X^k) = (c0(X^k), c1(X^k))` for `M(X^k)`
//     /// (i.e. Σ b_i • (X^i)^k) under S(X^k) as secret key.
//     ///
//     /// So the trick is as follows:
//     /// For `Subs(•,k)`, take ct(X) transform it to ct(X^k).
//     /// Create new `Ksk` with params `curr_sk=S(X^k)` and
//     /// `new_key=S(X)` and then perform key switch operation
//     /// so that you obtain M(X^k) encrypt under `new_key` that is
//     /// `S(X)`
//     ///
//     /// Also note that since we are operating on polynomial
//     /// integers modulus some cyclotomic polynomial `X^N + 1`,
//     /// following holds whenever k = N + 1
//     /// Subs(•, N + 1) transforms
//     /// `RLWE(Σ b_i • X^i)`
//     /// => `RLWE(Σ b_i • (X^i)^(N+1))`
//     /// => `RLWE(Σ b_i • (X^(Ni+i))` ---- (1)
//     /// Since X^N + 1 = 0 mod (X^N + 1)
//     /// => X^N = -1 mod (X^N + 1)
//     /// => (X^N)^i = X^(N*i) = -1^i mod (X^N + 1) ---- (2)
//     /// Therefore following from (1) and (2)
//     /// we can write (1) as following
//     /// => `RLWE(Σ b_i • -1^(i) • X^(i))`
//     ///
//     /// `ksk` is key switching key corresponding to
//     /// Subs(•,k).
//     pub fn substitute(ksk: &Ksk, ct: &BfvCipherText) -> BfvCipherText {
//         assert!(ksk.is_subs && ksk.subs_k.is_some());

//         // ct(X) -> ct(X^k)
//         let ct_k = BfvCipherText::new(
//             &ct.params,
//             ct.c.iter()
//                 .map(|poly| poly.shift_powers(ksk.subs_k.unwrap()))
//                 .collect(),
//         );

//         Ksk::key_switch(ksk, ct_k)
//     }

//     /// Generates key switching key for
//     /// Subs(•, k) operation
//     ///
//     /// Takes in secret key `sk` and returns
//     /// key switching key for necessary ops.
//     pub fn gen_substitution_key(sk: &BfvSecretKey, k: u64, beta: u64) -> Self {
//         // S(X^k)
//         let sk_k_poly = sk.poly.shift_powers(k);
//         let sk_k = BfvSecretKey::new_with_poly(&sk.params, sk_k_poly);
//         Ksk::new(&sk_k, sk, beta, true, Some(k))
//     }
// }

/// Structure:
/// Params: \ß & l
/// Note: Least significant first
/// [
///     [
///         RLWE( q/ß^l * -s * m)
///         RLWE( q/ß^{l-1} * -s * m)
///         .
///         .
///         .
///         RLWE( q/ß^1 * -s * m)
///     ]
///     [
///         RLWE( q/ß^l * m)
///         RLWE( q/ß^{l-1} * m)
///         .
///         .
///         .
///         RLWE( q/ß^{1} * m)
///     ]
/// ]
/// E (2l X 2)
#[derive(Debug)]
pub struct RgswCt {
    pub ksk0: Ksk,
    pub ksk1: Ksk,
}

impl RgswCt {
    pub fn encrypt_poly(sk: &SecretKey, m: &Poly) -> Self {
        assert!(sk.params.rq_context == m.context);
        assert!(m.representation == Representation::Ntt);

        let mut sk_poly = Poly::try_from_vec_i64(&sk.params.rq_context, &sk.coeffs);
        sk_poly.change_representation(Representation::Ntt);

        let s_m = &sk_poly * m;

        let ksk0 = Ksk::new(sk, &s_m);
        let ksk1 = Ksk::new(sk, m);

        RgswCt { ksk0, ksk1 }
    }

    pub fn encrypt(sk: &SecretKey, m: &BfvPlaintext) -> Self {
        assert!(sk.params == m.params);

        let mut m = Poly::try_from_vec_u64(&sk.params.rq_context, &m.values);
        m.change_representation(Representation::Ntt);

        RgswCt::encrypt_poly(sk, &m)
    }

    pub fn external_product(rgsw: &RgswCt, bfv: &BfvCipherText) -> (Poly, Poly) {
        let mut c1 = bfv.cts[1].clone();
        let mut c0 = bfv.cts[0].clone();
        c1.change_representation(Representation::PowerBasis);
        c0.change_representation(Representation::PowerBasis);

        // c1 * s * m
        let (c1_0, c1_1) = rgsw.ksk0.key_switch(&c1);
        // c0 * m
        let (mut c0_0, mut c0_1) = rgsw.ksk1.key_switch(&c0);

        c0_0 += &c1_0;
        c0_1 += &c1_1;

        (c0_0, c0_1)
    }

    pub fn try_from_ksks(ksk0: &Ksk, ksk1: &Ksk) -> Self {
        Self {
            ksk0: ksk0.clone(),
            ksk1: ksk1.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bfv::{BfvParameters, BfvPlaintext},
        poly::Modulus,
        rns::RnsContext,
        rq::RqContext,
        utils::generate_prime,
    };
    use itertools::izip;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use rand::thread_rng;
    use std::sync::Arc;

    #[test]
    fn encrypt_poly() {
        // let pt_moduli: u64 = 4611686018326724609;
        let pt_moduli = generate_prime(60, 2 * 1048576, (1 << 60) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[61, 61, 61], 64).unwrap();

        let params = Arc::new(BfvParameters::new(
            64, pt_moduli,
            // vec![1125899906840833, 36028797018963841, 36028797018963457],
            ct_moduli, 10,
        ));
        params.rq_context.rns.garner.iter().for_each(|gi| {
            dbg!(gi);
            dbg!(gi % pt_moduli);
        });

        let sk = SecretKey::generate(&params);
        let mut sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        sk_poly.change_representation(Representation::Ntt);

        let mut m = Poly::try_from_vec_i64(&sk.params.rq_context, &sk.coeffs);
        m.change_representation(Representation::Ntt);
        let sk_rgsw = RgswCt::encrypt_poly(&sk, &m);

        let pt_one = BfvPlaintext::new(&params, &vec![1u64]);
        let ex_rgsw = RgswCt::encrypt(&sk, &pt_one);

        let bis = [
            1084615932887104652u64,
            88464820691527941,
            1132762255576341106,
        ]
        .iter()
        .map(|gi| sk.encrypt(&BfvPlaintext::new(&params, &vec![*gi])));
        // let skis: Vec<BfvCipherText> = bis
        //     .map(|bi| {
        //         let (c0, c1) = RgswCt::external_product(&sk_rgsw, &bi);
        //         BfvCipherText {
        //             params: params.clone(),
        //             cts: vec![c0.clone(), c1.clone()],
        //         }
        //     })
        //     .collect();
        let skis: Vec<BfvCipherText> = izip!(ex_rgsw.ksk1.c0.iter(), ex_rgsw.ksk1.c1.iter())
            .map(|(c0, c1)| {
                let (c02, c12) = RgswCt::external_product(
                    &sk_rgsw,
                    &BfvCipherText {
                        params: params.clone(),
                        cts: vec![c0.clone(), c1.clone()],
                    },
                );
                BfvCipherText {
                    params: params.clone(),
                    cts: vec![c02.clone(), c12.clone()],
                }
            })
            .collect();

        // dbg!(ex_rgsw.ksk0.c0.len());
        // dbg!(ex_rgsw.ksk0.c1.len());
        // dbg!(skis.len());

        izip!(ex_rgsw.ksk0.c0.iter(), ex_rgsw.ksk0.c1.iter(), skis.iter()).for_each(
            |(ec0, ec1, cti)| {
                println!(
                    "Decrypted ex {:?}",
                    sk.decrypt(&BfvCipherText {
                        params: params.clone(),
                        cts: vec![ec0.clone(), ec1.clone()]
                    })
                    .values
                );

                println!(
                    "Decrypted ski {:?}",
                    sk.decrypt(&BfvCipherText {
                        params: params.clone(),
                        cts: vec![cti.cts[0].clone(), cti.cts[1].clone()]
                    })
                    .values
                );

                let mut ex_poly = ec1 * &sk_poly;
                ex_poly += ec0;

                let mut poly = &cti.cts[1] * &sk_poly;
                poly += &cti.cts[0];

                let mut sub = &ex_poly - &poly;

                ex_poly.change_representation(Representation::PowerBasis);
                poly.change_representation(Representation::PowerBasis);
                sub.change_representation(Representation::PowerBasis);

                let rns = RnsContext::new(&params.ciphertext_moduli.clone());
                izip!(
                    Vec::<BigUint>::from(&ex_poly).iter(),
                    Vec::<BigUint>::from(&poly).iter(),
                    Vec::<BigUint>::from(&sub).iter()
                )
                .for_each(|(ex, p, s)| {
                    println!(
                        "{} {} {}",
                        (rns.modulus() - ex).bits(),
                        (rns.modulus() - p).bits(),
                        (rns.modulus() - s).bits()
                    );
                });
            },
        );
    }

    #[test]
    fn external_product() {
        for _ in 0..100 {
            let mut rng = thread_rng();
            let params = Arc::new(BfvParameters::default(6, 8));
            let sk = SecretKey::generate(&params);

            let v1 = params.plaintext_modulus.random_vec(1, &mut rng);
            let v2 = params.plaintext_modulus.random_vec(1, &mut rng);

            let m1 = BfvPlaintext::new(&params, &v1);
            let m2 = BfvPlaintext::new(&params, &v2);

            let bfv_ct = sk.encrypt(&m1);
            let rgsw_ct = RgswCt::encrypt(&sk, &m2);

            let (ec0, ec1) = RgswCt::external_product(&rgsw_ct, &bfv_ct);
            let product = sk.decrypt(&BfvCipherText {
                cts: vec![ec0, ec1],
                params: params.clone(),
            });

            let rp = Arc::new(RqContext::new(
                vec![params.plaintext_modulus.p],
                params.degree,
            ));
            let mut v1_poly = Poly::try_from_vec_u64(&rp, &v1);
            let mut v2_poly = Poly::try_from_vec_u64(&rp, &v2);
            v1_poly.change_representation(Representation::Ntt);
            v2_poly.change_representation(Representation::Ntt);

            let mut expected_product = &v1_poly * &v2_poly;
            expected_product.change_representation(Representation::PowerBasis);
            let coeffs = expected_product.coefficients();
            let values = params
                .plaintext_modulus
                .reduce_vec_u64(coeffs.as_slice().unwrap());
            assert_eq!(product.values, values.into());
        }
    }
}
