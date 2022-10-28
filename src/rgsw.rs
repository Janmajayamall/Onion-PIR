use std::{sync::Arc, task::Poll};

use itertools::izip;
use num_bigint::BigUint;

use super::bfv::BfvCipherText;
use crate::{
    bfv::{BfvParameters, Plaintext, SecretKey},
    ksk::Ksk,
    poly::Modulus,
    rq::{BitDecomposition, Poly, Representation, RqContext},
};

use itertools::Itertools;

pub struct RGSW {
    pub sm: Vec<BfvCipherText>,
}

#[derive(Debug)]
pub struct RgswCt {
    pub ksk0: Ksk,
    pub ksk1: Ksk,
}

/// IMPORTANT SWITCH BFV TO RLWE (i.e. scaling of M isn't default)
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

    pub fn encrypt(sk: &SecretKey, m: &Plaintext) -> Self {
        assert!(sk.params == m.params);
        let params = sk.params.clone();

        let mut m_poly = Poly::try_from_vec_u64(&sk.params.rq_context, &m.values);
        m_poly.change_representation(Representation::Ntt);

        let mut s_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        s_poly.change_representation(Representation::Ntt);

        let sm_poly = &s_poly * &m_poly;

        RgswCt {
            ksk0: Ksk::new(sk, &sm_poly),
            ksk1: Ksk::new(sk, &m_poly),
        }
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
        bfv::{BfvParameters, Plaintext},
        poly::Modulus,
        rns::RnsContext,
        rq::{BitDecomposition, RqContext},
        utils::generate_prime,
    };

    use itertools::izip;
    use num_bigint::BigUint;
    use num_traits::ToPrimitive;
    use rand::thread_rng;
    use sha2::digest::typenum::Mod;
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

        let pt_one = Plaintext::new(&params, &vec![1u64]);
        let ex_rgsw = RgswCt::encrypt(&sk, &pt_one);

        let bis = [
            1084615932887104652u64,
            88464820691527941,
            1132762255576341106,
        ]
        .iter()
        .map(|gi| sk.encrypt(&Plaintext::new(&params, &vec![*gi])));
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
    fn construct_full_from_half() {
        let mut rng = thread_rng();
        let params = Arc::new(BfvParameters::default(
            6,
            8,
            BitDecomposition { base: 4, l: 8 },
        ));
        let sk = SecretKey::generate(&params);
        let mut sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        sk_poly.change_representation(Representation::Ntt);

        let sk_rgsw = RgswCt::encrypt_poly(&sk, &sk_poly);

        let b = Plaintext {
            params: params.clone(),
            values: vec![1].to_vec().into_boxed_slice(),
        };
        let b_ksk = Ksk::new_with_pt(&sk, &b);

        let upper_half = izip!(b_ksk.c0.iter(), b_ksk.c1.iter())
            .map(|(c0, c1)| {
                RgswCt::external_product(
                    &sk_rgsw,
                    &BfvCipherText {
                        params: params.clone(),
                        cts: vec![c0.clone(), c1.clone()],
                    },
                )
            })
            .map(|(c0, c1)| BfvCipherText {
                params: params.clone(),
                cts: vec![c0, c1],
            })
            .collect_vec();
        let upper_half = Ksk::try_from_decomposed_rlwes(&params, &upper_half);

        let m = Modulus::new(params.plaintext_modulus_u64).random_vec(params.degree, &mut rng);
        let m_ct = sk.encrypt(&Plaintext::new(&params, &m));
        let b_rgsw_constructed = RgswCt::try_from_ksks(&upper_half, &b_ksk);
        let product = RgswCt::external_product(&b_rgsw_constructed, &m_ct);
        let product = sk.decrypt(&BfvCipherText {
            params: params.clone(),
            cts: vec![product.0.clone(), product.1.clone()],
        });
        assert_eq!(product.values, m.into());
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
        let sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);

        // let m = BfvPlaintext::new(&params, &vec![2]);
        // let rgsw_ct = RgswCt::encrypt(&sk, &m);
        let mut m = &Plaintext::new(&params, &vec![2]);
        let rgsw_ct = RgswCt::encrypt(&sk, m);

        let ksk1 = Ksk::new_with_pt(&sk, &Plaintext::new(&params, &vec![2]));

        // construct ksk2 from ksk1 using rgsw
        let ksk2 = izip!(ksk1.c0.iter(), ksk1.c1.iter()).map(|(c0, c1)| {
            RgswCt::external_product(
                &rgsw_ct,
                &BfvCipherText {
                    params: params.clone(),
                    cts: vec![c0.clone(), c1.clone()],
                },
            )
        });

        let ksk2_s = izip!(ksk1.c0.iter(), ksk1.c1.iter()).map(|(c0, c1)| {
            let f = &BfvCipherText {
                params: params.clone(),
                cts: vec![c0.clone(), c1.clone()],
            } * &sk_poly;
            RgswCt::external_product(&rgsw_ct, &f)
        });
        let ksk2 = Ksk {
            params: params.clone(),
            c0: ksk2.clone().map(|(c0, _)| c0).collect(),
            c1: ksk2.map(|(_, c1)| c1).collect(),
            ct_context: params.rq_context.clone(),
        };

        let p1 = Poly::try_from_vec_u64(&params.rq_context, &vec![3]);
        let product = ksk2.key_switch(&p1);

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

            let v1 = params.plaintext_modulus.random_vec(params.degree, &mut rng);
            let v2 = params.plaintext_modulus.random_vec(params.degree, &mut rng);

            let m1 = Plaintext::new(&params, &v1);
            let m2 = Plaintext::new(&params, &v2);

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

fn record_difference(real: &Poly, ideal: &Poly) {
    let mut diff = real - ideal;
    diff.change_representation(Representation::PowerBasis);

    Vec::<BigUint>::from(&diff).iter().for_each(|coeff| {
        dbg!(std::cmp::min(
            coeff.bits(),
            real.context.rns.product.bits() - coeff.bits()
        ));
    });
}

fn measure_noise_tmp(
    params: &Arc<BfvParameters>,
    ideal_m: Poly,
    sk: &SecretKey,
    ct: &BfvCipherText,
) -> usize {
    let pt = Poly::try_from_vec_u64(&ct.params.rq_context, &[1u64]);

    // TODO: REMOVE
    let constant = sk.params.rq_context.rns.garner[0].clone();

    let mut bg = Poly::try_from_bigint(&sk.params.rq_context, &[constant.clone()]);
    bg.change_representation(Representation::Ntt);

    let mut sk = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
    sk.change_representation(Representation::Ntt);
    let mut m = &ct.cts[1] * &sk;
    m += &ct.cts[0];

    bg -= &ideal_m;
    bg.change_representation(Representation::PowerBasis);

    let mut noise = 0usize;
    let ct_moduli = &params.rq_context.rns.product.clone();

    Vec::<BigUint>::from(&bg).iter().for_each(|coeff| {
        noise = std::cmp::max(
            noise,
            std::cmp::min(coeff.bits(), (ct_moduli - coeff).bits()) as usize,
        );
    });
    noise
}
