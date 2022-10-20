use std::sync::Arc;

use crate::rq::Representation;
use itertools::izip;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::{
    bfv::{BfvParameters, SecretKey},
    rq::{Poly, RqContext},
};

struct Ksk {
    params: Arc<BfvParameters>,

    c0: Vec<Poly>,
    c1: Vec<Poly>,

    /// Context of poly that will be key switched
    ct_context: Arc<RqContext>,
}

impl Ksk {
    pub fn new(sk: SecretKey, from: &Poly) -> Self {
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
}
