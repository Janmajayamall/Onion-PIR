use itertools::{izip, Itertools};
use num_traits::ToPrimitive;
use std::{
    collections::HashMap, default::Default, f32::consts::E, future::pending, ops::Sub,
    os::unix::process::parent_id, sync::Arc, vec,
};

use crate::{
    bfv::{BfvCipherText, BfvParameters, BfvPlaintext, GaliosKey, SecretKey},
    poly::Modulus,
    rq::{Poly, Representation, RqContext, Substitution},
    utils::ilog2,
};

struct Client {
    first_dim: u64,
    dim: u64,
}

impl Default for Client {
    fn default() -> Self {
        Client {
            first_dim: 4,
            dim: 2,
        }
    }
}

impl Client {
    pub fn new(first_dim: u64, dim: u64) -> Client {
        Client { first_dim, dim }
    }

    pub fn encode(
        &self,
        query_dims: Vec<u64>,
        db_rows: u64,
        bfv_params: &Arc<BfvParameters>,
    ) -> Vec<u64> {
        assert!(query_dims[0] < self.first_dim);
        let mut first_dim_bfv: Vec<u64> = (0..self.first_dim)
            .map(|i| (i == query_dims[0]).into())
            .collect();

        assert!(((db_rows / self.first_dim) / 4u64.pow((query_dims.len() - 1) as u32)) == 1);

        // final length of `rgsws` = garner.len() * self.dim * (query_dims - 1)
        let mut rgsws: Vec<u64> = vec![];
        for i in 1..query_dims.len() {
            assert!(query_dims[i] < self.dim);
            for j in 0..self.dim {
                if j == query_dims[i] {
                    izip!(bfv_params.rq_context.rns.garner.iter()).for_each(|g| {
                        rgsws.push(bfv_params.plaintext_modulus.reduce_biguint(g));
                    });
                } else {
                    (0..(bfv_params.rq_context.rns.garner.len()))
                        .into_iter()
                        .for_each(|_| rgsws.push(0u64));
                }
            }
        }

        first_dim_bfv.extend(rgsws.iter());
        first_dim_bfv
    }
}

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
        let n = ct.params.degree;
        let mut bits_ct = vec![ct.clone(); n];

        let x_pows = Evaluation::gen_x_pows(ct.params.degree, &ct.params.rq_context);
        debug_assert!(x_pows.len() == ilog2(ct.params.degree));

        for i in 0..ilog2(n) {
            let exponent = (n / (1 << i)) + 1;

            for j in 0..(1 << i) {
                let c = bits_ct[j].clone();

                let gk = self.gk.get(&exponent).unwrap();

                let c_r = gk.relinearize(&bits_ct[j]);

                let c_even = &c + &c_r;

                let mut x_p = x_pows[i].clone();
                x_p.change_representation(Representation::Ntt);
                let mut c_odd = &c - &c_r;
                let c0_odd = &c_odd.cts[0] * &x_p;
                let c1_odd = &c_odd.cts[1] * &x_p;
                c_odd.cts = vec![c0_odd, c1_odd];

                bits_ct[j] = c_even;
                bits_ct[j + (1 << i)] = c_odd;
            }
        }

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
            polys.push(xi);
        }
        polys
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bfv::{BfvParameters, BfvPlaintext},
        rq::Representation,
    };
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use rand::{thread_rng, Rng};
    use rand_distr::Uniform;

    #[test]
    fn trial() {
        let params = Arc::new(BfvParameters::default(6, 8));
        let mut v = Poly::try_from_vec_u64(&params.rq_context, &[1]);
        v.change_representation(Representation::Ntt);
        let f = BigUint::from_usize(10usize).unwrap();
        v *= &f;
        v.change_representation(Representation::PowerBasis);
        dbg!(v);
    }

    #[test]
    fn client_encode() {
        let params = Arc::new(BfvParameters::default(2, 64));

        let query_dims = [2u64];
        let client = Client::default();
        let query = client.encode(query_dims.to_vec(), 64, &params);
        dbg!(&query);
        let sk = SecretKey::generate(&params);
        let evaluation = Evaluation::new(&sk, &params);

        let pt = BfvPlaintext::new(&params, &query);
        let ct = sk.encrypt(&pt);

        let cts = evaluation.unpack(&ct);
        let mut query_r = vec![];
        cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            query_r.push(p.values[0]);
        });
        dbg!(query);
        dbg!(query_r);
    }

    #[test]
    fn unpack() {
        let mut rng = thread_rng();
        let degree = 2048;
        let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        // let params = Arc::new(BfvParameters::new(
        //     64,
        //     1153,
        //     [
        //         4611686018427387761,
        //         4611686018427387617,
        //         4611686018427387409,
        //         2305843009213693921,
        //         1152921504606846577,
        //         2017,
        //     ]
        //     .to_vec(),
        //     10,
        // ));
        // let params = Arc::new(BfvParameters::default(2, 2048));
        let sk = SecretKey::generate(&params);

        let evaluation = Evaluation::new(&sk, &params);

        let binary_dist = Uniform::from(0..pt_moduli);
        let values = rng.sample_iter(binary_dist).take(2048).collect();
        let pt = BfvPlaintext::new(&params, &values);

        let ct = sk.encrypt(&pt);
        let d_values = sk.decrypt(&ct);

        let unpacked_cts = evaluation.unpack(&ct);
        let mut r_values = vec![];
        unpacked_cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            let v = p.values[0] / (degree as u64);
            r_values.push(v);
        });

        // error counts
        let mut count = 0;
        let mut d_count = 0;
        izip!(values.iter(), r_values.iter(), d_values.values.iter()).for_each(|(v, r_v, d_v)| {
            if v != r_v {
                count += 1;
            }
            if v != d_v {
                d_count += 1;
            }
        });

        // println!("{:?}", values);
        // println!("{:?}", r_values);
        dbg!(count);
        dbg!(d_count);
    }
}
