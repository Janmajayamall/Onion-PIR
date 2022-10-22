use std::{collections::HashMap, f32::consts::E, ops::Sub, sync::Arc, vec};

use crate::{
    bfv::{BfvCipherText, BfvParameters, BfvPlaintext, GaliosKey, SecretKey},
    rq::{Poly, RqContext, Substitution},
    utils::ilog2,
};

struct Evaluation {
    gk: HashMap<usize, GaliosKey>,
}

impl Evaluation {
    pub fn new(sk: &SecretKey, params: &Arc<BfvParameters>) -> Self {
        let mut gk_map = HashMap::new();
        let n = params.degree;
        for i in 0..ilog2(n) {
            let k = (n / 2usize.pow(i as u32)) + 1;
            let gk_i = GaliosKey::new(sk, &Substitution::new(k));
            gk_map.insert(k, gk_i);
        }

        Evaluation { gk: gk_map }
    }

    pub fn unpack(&self, ct: &BfvCipherText) -> Vec<BfvCipherText> {
        let mut bits_ct = vec![ct.clone()];

        let x_pows = Evaluation::gen_x_pows(ct.params.degree, &ct.params.rq_context);
        debug_assert!(x_pows.len() == ilog2(ct.params.degree));

        let n = ct.params.degree;
        for i in 0..ilog2(n) {
            let mut branch_cts = Vec::<BfvCipherText>::new();
            for j in 0..bits_ct.len() {
                let c = &bits_ct[j];

                let exponent = (n / 2usize.pow(i as u32)) + 1;
                let gk = self.gk.get(&exponent).unwrap();

                let c_r = gk.relinearize(&bits_ct[j]);

                let c_even = c + &c_r;
                let mut c_odd = c - &c_r;
                c_odd.cts[0] *= &x_pows[i];
                c_odd.cts[1] *= &x_pows[i];

                branch_cts.push(c_even);
                branch_cts.push(c_odd);
            }
            bits_ct = branch_cts;
        }

        dbg!(bits_ct.len());
        assert!(bits_ct.len() == ct.params.degree);

        bits_ct
    }

    fn gen_x_pows(n: usize, ctx: &Arc<RqContext>) -> Vec<Poly> {
        let mut polys = Vec::<Poly>::new();
        // X
        let x = Poly::try_from_vec_u64(&ctx, &[0u64, 1]);
        for i in 0..ilog2(n) {
            let mut xi = x.clone();
            // X^-(2^i) = -X^(N - 2^i), since X^-k = -(X^(N-k))
            let ex = n - 2usize.pow(i as u32);
            xi = -&xi.substitute(&Substitution::new(ex));
            xi.change_representation(crate::rq::Representation::Ntt);
            polys.push(xi);
        }
        polys
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{BfvParameters, BfvPlaintext};

    #[test]
    fn unpack() {
        let params = Arc::new(BfvParameters::new(
            8,
            1153,
            [
                4611686018427387761,
                4611686018427387617,
                4611686018427387409,
                2305843009213693921,
                1152921504606846577,
                2017,
            ]
            .to_vec(),
            10,
        ));
        let sk = SecretKey::generate(&params);

        let evaluation = Evaluation::new(&sk, &params);

        let v = vec![1u64, 0, 1, 1, 1, 1, 1, 0];
        let pt = BfvPlaintext::new(&params, &v);

        let ct = sk.encrypt(&pt);

        let unpacked_cts = evaluation.unpack(&ct);
        unpacked_cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            dbg!(p);
        });
    }
}
