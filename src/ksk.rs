use std::sync::Arc;

use crate::{bfv::BfvCipherText, rq::Representation};
use itertools::{izip, Itertools};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::{
    bfv::{BfvParameters, SecretKey},
    rq::{Poly, RqContext},
};

#[derive(Clone, Debug)]
pub struct Ksk {
    pub params: Arc<BfvParameters>,

    c0: Vec<Poly>,
    c1: Vec<Poly>,

    /// Context of poly that will be key switched
    ct_context: Arc<RqContext>,
}

impl Ksk {
    pub fn new(sk: &SecretKey, from: &Poly) -> Self {
        let params = sk.params.clone();

        let c1s: Vec<Poly> = (0..params.rq_context.rns.moduli.len())
            .map(|_| {
                let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
                thread_rng().fill(&mut seed);
                Poly::random_from_seed(&params.rq_context, Representation::Ntt, seed)
            })
            .collect();

        let mut sk = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        sk.change_representation(Representation::Ntt);

        debug_assert!(c1s.len() == params.rq_context.rns.garner.len());
        let mut c0s = Vec::<Poly>::with_capacity(c1s.len());
        izip!(c1s.iter(), params.rq_context.rns.garner.iter()).for_each(|(c1, garner)| {
            let mut a_s = c1.clone();
            a_s *= &sk;

            let mut rng = thread_rng();
            // e
            let mut b = Poly::random_small(
                &params.rq_context,
                Representation::Ntt,
                params.variance,
                &mut rng,
            );
            // -(a*s) + e
            b -= &a_s;

            // m = garner * from
            let mut m = from.clone();
            m *= garner;

            // -(a*s) + e + m
            b += &m;

            c0s.push(b);
        });

        Ksk {
            params: params.clone(),
            c0: c0s,
            c1: c1s,
            ct_context: params.rq_context.clone(),
        }
    }

    pub fn key_switch(&self, poly: &Poly) -> (Poly, Poly) {
        assert!(poly.context == self.ct_context);
        assert!(poly.representation == Representation::PowerBasis);

        let mut c0 = Poly::zero(&self.ct_context, Representation::Ntt);
        let mut c1 = Poly::zero(&self.ct_context, Representation::Ntt);
        izip!(
            poly.coefficients().outer_iter(),
            self.c0.iter(),
            self.c1.iter()
        )
        .for_each(|(decomposed_poly_i, c0_i, c1_i)| {
            let mut decomposed_poly_i =
                Poly::try_from_vec_u64(&self.ct_context, decomposed_poly_i.as_slice().unwrap());
            decomposed_poly_i.change_representation(Representation::Ntt);

            let b = &decomposed_poly_i * c0_i;
            let a = &decomposed_poly_i * c1_i;
            c0 += &b;
            c1 += &a;
        });

        (c0, c1)
    }

    pub fn try_from_decomposed_rlwes(
        params: &Arc<BfvParameters>,
        cts: &Vec<BfvCipherText>,
    ) -> Self {
        assert!(params.rq_context.moduli.len() == cts.len());

        let c0 = cts
            .iter()
            .map(|c| {
                assert!(&c.params == params);
                c.cts[0].clone()
            })
            .collect();
        let c1 = cts
            .iter()
            .map(|c| {
                assert!(&c.params == params);
                c.cts[1].clone()
            })
            .collect();

        Ksk {
            params: params.clone(),
            c0,
            c1,
            ct_context: params.rq_context.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::*;
    use crate::{
        bfv::{BfvParameters, BfvPlaintext},
        rns::RnsContext,
    };

    #[test]
    fn trial() {
        let mut rng = thread_rng();
        let degree = 64;
        // let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        // dbg!(&ct_moduli);
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        // let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        let params = Arc::new(BfvParameters::new(
            degree,
            pt_moduli,
            vec![1125899906840833, 36028797018963841, 36028797018963457],
            10,
        ));
        let sk = &SecretKey::generate(&params);

        let mut one_poly = Poly::try_from_vec_u64(&params.rq_context, &[0u64]);
        let mut sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        one_poly.change_representation(Representation::Ntt);
        sk_poly.change_representation(Representation::Ntt);
        sk_poly *= &one_poly;
        let ksk = Ksk::new(&sk, &sk_poly);
        izip!(ksk.c0.iter(), ksk.c1.iter()).for_each(|(c0, c1)| {
            let p = sk.decrypt(&BfvCipherText {
                params: params.clone(),
                cts: vec![c0.clone(), c1.clone()],
            });
            println!("{:?}", p.values);
        });

        let mut two = BfvPlaintext::new(&params, &[2u64].to_vec());
        let mut two_poly = sk.encrypt(&two);
        let mut c0 = two_poly.cts[0].clone();
        let mut c1 = two_poly.cts[1].clone();
        c1.change_representation(Representation::PowerBasis);
        let (c0_2, c1_2) = ksk.key_switch(&c1);
        println!(
            "two * sk: {:?}",
            sk.decrypt(&BfvCipherText {
                params: params.clone(),
                cts: vec![&c0 + &c0_2, c1_2],
            })
            .values
        );
    }

    #[test]
    fn ksk() {
        for _ in 0..100 {
            let params = Arc::new(BfvParameters::default(2, 8));
            let sk = SecretKey::generate(&params);
            let mut sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
            sk_poly.change_representation(Representation::Ntt);

            let mut rng = thread_rng();
            let mut from_poly = Poly::random(&params.rq_context, &mut rng, Representation::Ntt);
            let ksk = Ksk::new(&sk, &from_poly);

            let mut input = Poly::random(&params.rq_context, &mut rng, Representation::PowerBasis);
            let (c0, c1) = ksk.key_switch(&input);

            let mut product = &c1 * &sk_poly;
            product += &c0;

            from_poly.change_representation(Representation::Ntt);
            input.change_representation(Representation::Ntt);
            let mut product_ex = &from_poly * &input;

            product_ex -= &product;

            product_ex.change_representation(Representation::PowerBasis);

            // dbg!(product_ex.coefficients());

            let rns = RnsContext::new(&params.ciphertext_moduli.clone());
            Vec::<BigUint>::from(&product_ex).iter().for_each(|b| {
                // TODO: Why does this need to be less than 70 ?
                // Ref - https://github.com/tlepoint/fhe.rs/blob/f7cddb358f2ce28483944f99e223c07ae41b0c1c/crates/fhe/src/bfv/keys/key_switching_key.rs#L335
                assert!((std::cmp::min(b.bits(), (rns.modulus() - b).bits())) <= 70);
            });
        }
    }
}
