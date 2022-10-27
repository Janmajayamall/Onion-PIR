use std::{sync::Arc, task::Poll};

use super::bfv::BfvCipherText;
use crate::{
    bfv::{BfvParameters, BfvPlaintext, SecretKey},
    ksk::Ksk,
    poly::Modulus,
    rq::{Poly, Representation},
};

pub struct RGSW {
    pub sm: Vec<BfvCipherText>,
}

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

    // pub fn encrypt_pt(sk: &SecretKey, pt: &BfvPlaintext) {
    //     let sk_m = vec![];
    //     let pt_mod = Modulus::new(sk.params.plaintext_modulus_u64);
    //     let sk_coeffs = sk.coeffs;
    //     let m = pt.values;
    //     s
    // }

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
        rq::{BitDecomposition, RqContext},
        utils::generate_prime,
    };
    use itertools::izip;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use rand::thread_rng;
    use std::{sync::Arc, vec};

    #[test]
    fn encrypt_poly() {
        // let pt_moduli: u64 = 4611686018326724609;
        let pt_moduli = generate_prime(60, 2 * 1048576, (1 << 60) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[61, 61, 61], 64).unwrap();

        let params = Arc::new(BfvParameters::new(
            64,
            pt_moduli,
            // vec![1125899906840833, 36028797018963841, 36028797018963457],
            ct_moduli,
            10,
            BitDecomposition { base: 4, l: 8 },
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
    fn contruct_ksk_from_rgsw() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(
            6,
            8,
            BitDecomposition { base: 4, l: 8 },
        ));
        let sk = SecretKey::generate(&params);

        let m = BfvPlaintext::new(&params, &vec![2]);
        let rgsw_ct = RgswCt::encrypt(&sk, &m);

        let one_pt = BfvPlaintext::new(&params, &vec![1]);
        let ksk = Ksk::new_with_pt(&sk, &one_pt);

        let two_by_gi = izip!(ksk.c0.iter(), ksk.c1.iter()).map(|(c0, c1)| {
            RgswCt::external_product(
                &rgsw_ct,
                &BfvCipherText {
                    params: params.clone(),
                    cts: vec![c0.clone(), c1.clone()],
                },
            )
        });

        let two_ksk = Ksk {
            params: params.clone(),
            c0: two_by_gi.clone().map(|(c0, _)| c0).collect(),
            c1: two_by_gi.map(|(_, c1)| c1).collect(),
            ct_context: params.rq_context.clone(),
        };

        let three = Poly::try_from_vec_u64(&params.rq_context, &vec![3]);
        let product = two_ksk.key_switch(&three);

        dbg!(sk.decrypt(&BfvCipherText {
            params: params.clone(),
            cts: vec![product.0, product.1]
        }));
    }

    #[test]
    fn external_product() {
        for _ in 0..100 {
            let mut rng = thread_rng();
            let params = Arc::new(BfvParameters::default(
                6,
                8,
                BitDecomposition { base: 4, l: 8 },
            ));
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
                BitDecomposition { base: 4, l: 8 },
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
