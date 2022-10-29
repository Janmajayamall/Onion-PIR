use crate::{
    bfv::{BfvCipherText, BfvParameters, GaliosKey, Plaintext, SecretKey},
    ksk::Ksk,
    modulus::Modulus,
    rgsw::RgswCt,
    rq::{Poly, Representation, RqContext, Substitution},
    utils::ilog2,
};
use itertools::{enumerate, izip, Itertools};
use ndarray::Array2;
use num_bigint::BigUint;
use num_bigint_dig::{BigUint as BigUintDig, ModInverse};
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use sha2::digest::typenum::bit;
use std::{collections::HashMap, default::Default, sync::Arc, vec};

struct Client {
    first_dim: usize,
    dim: usize,
}

impl Default for Client {
    fn default() -> Self {
        Client {
            first_dim: 4,
            dim: 4,
        }
    }
}

impl Client {
    pub fn new(first_dim: usize, dim: usize) -> Client {
        Client { first_dim, dim }
    }

    pub fn encode(
        &self,
        query_dims: Vec<usize>,
        db_rows: usize,
        bfv_params: &Arc<BfvParameters>,
        sk: &SecretKey,
    ) -> Vec<u64> {
        assert!(query_dims[0] < self.first_dim);
        let mut first_dim_bfv: Vec<u64> = (0..self.first_dim)
            .map(|i| (i == query_dims[0]).into())
            .collect();

        assert!(((db_rows / self.first_dim) / self.dim.pow((query_dims.len() - 1) as u32)) == 1);

        // final length of `rgsws` = garner.len() * self.dim * (query_dims - 1)
        let mut rgsws: Vec<u64> = vec![];
        for i in 1..query_dims.len() {
            assert!(query_dims[i] < self.dim);
            for j in 0..self.dim {
                rgsws.push((j == query_dims[i]) as u64)
            }
        }

        first_dim_bfv.extend(rgsws.iter());
        // dbg!(&first_dim_bfv);
        first_dim_bfv
    }
}

struct Server {
    first_dim: usize,
    dim: usize,
    db: Vec<Plaintext>,
    db_scaled: Vec<Poly>,
}

impl Server {
    pub fn new(
        first_dim: usize,
        dim: usize,
        db: &Vec<Plaintext>,
        params: &Arc<BfvParameters>,
    ) -> Self {
        let db_scaled = db
            .iter()
            .map(|p| {
                let mut pt = Poly::try_from_vec_u64(&params.rq_context, &p.values);
                pt.change_representation(Representation::Ntt);
                pt
            })
            .collect();
        Self {
            first_dim,
            dim,
            db: db.clone(),
            db_scaled,
        }
    }

    fn no_of_cts(&self) -> usize {
        let mut len = self.db.len();
        let mut count = 0usize;

        count += self.first_dim as usize;
        len /= self.first_dim as usize;

        while len > 1 {
            len /= self.dim as usize;
            count += self.dim as usize;
        }
        count
    }

    pub fn decode(
        &self,
        query_cts: Vec<BfvCipherText>,
        eval_key: &Evaluation,
        params: &Arc<BfvParameters>,
    ) -> Vec<BfvCipherText> {
        assert!(params.rq_context.moduli.len() + 1 == query_cts.len());

        let mut first_dim_cts = eval_key.unpack(query_cts[0].clone());
        first_dim_cts = first_dim_cts.as_slice()[..self.first_dim].to_vec();

        let dimensions = ilog2(self.db.len() / self.first_dim) / ilog2(self.dim) + 1;

        let rgsws_count = (dimensions - 1) * self.dim;
        let mut rgsws = vec![];
        let mut vals = vec![];
        for i in 1..query_cts.len() {
            let cts = eval_key.unpack(query_cts[i].clone()).as_slice()[..rgsws_count].to_vec();
            vals.push(cts);
        }
        for i in 0..rgsws_count {
            let mut lower_half = vec![];
            for j in 0..vals.len() {
                lower_half.push(vals[j][i].clone());
            }
            let upper_half = lower_half
                .iter()
                .map(|ct| {
                    let (c0, c1) = RgswCt::external_product(&eval_key.sk_rgsw, ct);
                    BfvCipherText {
                        params: params.clone(),
                        cts: vec![c0, c1],
                    }
                })
                .collect();
            let ksk1 = Ksk::try_from_decomposed_rlwes(params, &lower_half);
            let ksk0 = Ksk::try_from_decomposed_rlwes(params, &upper_half);
            rgsws.push(RgswCt::try_from_ksks(&ksk0, &ksk1));
        }

        let rgsws_by_dim = rgsws
            .iter()
            .chunks(self.dim)
            .into_iter()
            .map(|cts| cts.into_iter().collect_vec())
            .collect_vec();

        // dbg!(dimensions_queries.len());
        // dbg!(dimensions_queries[0].len());

        // First dimension
        let mut db: Vec<BfvCipherText> = self
            .db_scaled
            .chunks(self.first_dim)
            .into_iter()
            .map(|pt_polys| {
                let rows = izip!(first_dim_cts.iter(), pt_polys.iter()).map(|(q, p)| q * p);

                let mut sum = BfvCipherText {
                    params: params.clone(),
                    cts: vec![],
                };
                rows.enumerate().for_each(|(index, product)| {
                    if index == 0 {
                        sum = product;
                    } else {
                        sum = &sum + &product
                    }
                });
                sum
            })
            .collect();

        // dbg!(&db.len());
        // dbg!(sk.decrypt(&db[0].clone()).values);
        // dbg!(sk.decrypt(&db[1].clone()).values);

        // Rest of dimensions
        rgsws_by_dim
            .iter()
            .enumerate()
            .for_each(|(dim_index, dim_rgsws)| {
                let tmp = db
                    .chunks(self.dim)
                    .into_iter()
                    .enumerate()
                    .map(|(row_index, db_chunk_cts)| {
                        let new_chunk_cts =
                            izip!(dim_rgsws.iter(), db_chunk_cts.iter()).map(|(rgsw, rlwe)| {
                                let (c0, c1) = RgswCt::external_product(rgsw, rlwe);
                                BfvCipherText {
                                    params: params.clone(),
                                    cts: vec![c0, c1],
                                }
                            });
                        let mut sum = BfvCipherText {
                            params: params.clone(),
                            cts: vec![],
                        };
                        new_chunk_cts.enumerate().for_each(|(index, product)| {
                            if index == 0 {
                                sum = product;
                            } else {
                                sum = &sum + &product
                            }
                        });
                        sum
                    })
                    .collect();
                db = tmp;
            });

        db
    }
}

struct Evaluation {
    gk: HashMap<usize, GaliosKey>,
    sk_rgsw: RgswCt,
}

impl Evaluation {
    pub fn build(sk: &SecretKey, params: &Arc<BfvParameters>) -> Self {
        let mut gk_map = HashMap::new();
        let n = params.degree;
        for i in 0..ilog2(n) {
            let k = (n / 2usize.pow(i as u32)) + 1;
            let gk_i = GaliosKey::new(sk, &Substitution::new(k));
            gk_map.insert(k, gk_i);
        }

        let mut m = Poly::try_from_vec_i64(&sk.params.rq_context, &sk.coeffs);
        m.change_representation(Representation::Ntt);
        let sk_rgsw = RgswCt::encrypt_poly(sk, &m);

        Evaluation {
            gk: gk_map,
            sk_rgsw,
        }
    }

    pub fn unpack(&self, mut ct: BfvCipherText) -> Vec<BfvCipherText> {
        let x_pows = Evaluation::gen_x_pows(ct.params.degree, &ct.params.rq_context);
        debug_assert!(x_pows.len() == ilog2(ct.params.degree));

        let n = ct.params.degree;

        // multiply ct by N^-1
        let n_biguint = BigUintDig::from_usize(n).unwrap();
        let n_inv = n_biguint
            .mod_inverse(ct.params.rq_context.rns.product_dig.clone())
            .unwrap()
            .to_biguint()
            .unwrap()
            .to_bytes_le();
        let n_inv = BigUint::from_bytes_le(&n_inv);

        ct.cts.iter_mut().for_each(|c| *c *= &n_inv);
        let mut bits_ct = vec![ct.clone(); n];

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
        bfv::{BfvParameters, Plaintext},
        rq::{BitDecomposition, Representation},
        utils::generate_prime,
    };
    use ndarray::Axis;
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use rand::{thread_rng, Rng};
    use rand_distr::Uniform;

    #[test]
    fn pack_and_construct_rgsw() {
        let mut rng = thread_rng();
        let degree = 64;
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62], degree).unwrap();
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));
        let sk = SecretKey::generate(&params);
        let mut sk_poly = Poly::try_from_vec_i64(&params.rq_context, &sk.coeffs);
        sk_poly.change_representation(Representation::Ntt);

        let query = vec![1u64, 0, 1, 0];
        let mut query_poly = Poly::try_from_vec_u64(&params.rq_context, &query);
        query_poly.change_representation(Representation::Ntt);
        let decomp_query_polys = params
            .rq_context
            .rns
            .garner
            .iter()
            .map(|gi| &query_poly * gi);

        let evaluation_key = Evaluation::build(&sk, &params);

        let decomp_query_cts = decomp_query_polys.map(|qp| sk.encrypt_poly(&qp));

        let unpacked_query_cts = decomp_query_cts
            .map(|ct| evaluation_key.unpack(ct))
            .collect_vec();

        // lower half
        let lower_half = params
            .rq_context
            .rns
            .garner
            .iter()
            .enumerate()
            .map(|(index, gi)| unpacked_query_cts[index][0].clone())
            .collect_vec();
        // upper half
        let upper_half = lower_half
            .iter()
            .map(|gi_b| {
                let (c0, c1) = RgswCt::external_product(&evaluation_key.sk_rgsw, gi_b);
                BfvCipherText {
                    params: params.clone(),
                    cts: vec![c0, c1],
                }
            })
            .collect();

        let lower_half = Ksk::try_from_decomposed_rlwes(&params, &lower_half);
        let upper_half = Ksk::try_from_decomposed_rlwes(&params, &upper_half);

        let rgsw = RgswCt::try_from_ksks(&upper_half, &lower_half);

        let mut p_1 = Poly::try_from_vec_u64(&params.rq_context, &[1]);
        p_1.change_representation(Representation::Ntt);
        let rgsw_ideal = RgswCt::encrypt_poly(&sk, &p_1);

        let ct1 = sk.encrypt(&Plaintext::new(
            &params,
            &params.plaintext_modulus.random_vec(degree, &mut rng),
        ));

        let product = RgswCt::external_product(&rgsw, &ct1);
        let product_ideal = RgswCt::external_product(&rgsw_ideal, &ct1);

        assert_eq!(
            sk.decrypt(&BfvCipherText {
                params: params.clone(),
                cts: vec![product.0.clone(), product.1.clone()],
            })
            .values,
            sk.decrypt(&BfvCipherText {
                params: params.clone(),
                cts: vec![product_ideal.0.clone(), product_ideal.1.clone()],
            })
            .values
        );

        // dbg!(sk.measure_noise(&BfvCipherText {
        //     params: params.clone(),
        //     cts: vec![product.0.clone(), product.1.clone()],
        // }));
    }

    #[test]
    fn query() {
        let mut rng = thread_rng();
        let degree = 64;
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62], degree).unwrap();
        // let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        // Pre-processing
        let db: Vec<Plaintext> = (0..16)
            .into_iter()
            .map(|i| Plaintext::new(&params, &vec![i, i, i, i]))
            .collect();

        // Client side
        let sk = &SecretKey::generate(&params);
        let query = vec![3usize, 1, 1];
        let client = Client::new(4, 2);
        let query_encoded = client.encode(query, db.len(), &params, &sk);
        dbg!(&query_encoded);
        let query_first_dim = Plaintext::new(&params, &query_encoded);
        let query_first_dim = sk.encrypt(&query_first_dim);
        let mut query_cts = vec![query_first_dim];
        params.rq_context.rns.garner.iter().for_each(|gi| {
            let mut p = Poly::try_from_vec_u64(&params.rq_context, &query_encoded.as_slice()[4..]);
            p.change_representation(Representation::Ntt);
            p *= gi;
            let ct = sk.encrypt_poly(&p);
            query_cts.push(ct);
        });

        // let query_ct = sk.encrypt_poly(&query_poly);
        let eval_key = Evaluation::build(&sk, &params);

        // Server side
        let server = Server::new(4, 2, &db, &params);
        let res_ct = server.decode(query_cts, &eval_key, &params);

        // Client
        dbg!(res_ct.len());
        let res_ct = res_ct[0].clone();
        dbg!(sk.decrypt(&res_ct).values);
    }

    #[test]
    fn client_encode() {
        let mut rng = thread_rng();
        let degree = 2048;
        let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let query_dims = [123usize, 2, 3];
        let client = Client::default();
        let query = client.encode(query_dims.to_vec(), 2048, &params, &sk);
        dbg!(&query);
        let sk = SecretKey::generate(&params);
        let evaluation = Evaluation::build(&sk, &params);

        let pt = Plaintext::new(&params, &query);
        let ct = sk.encrypt(&pt);

        let cts = evaluation.unpack(ct);
        let mut query_r = vec![];
        cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            query_r.push(p.values[0]);
        });
        dbg!(&query);
        dbg!(&query_r[..query.len()]);
    }

    #[test]
    fn unpack_decomp() {
        let mut rng = thread_rng();
        let degree = 8;
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62], degree).unwrap();
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let evaluation = Evaluation::build(&sk, &params);

        let binary_dist: Uniform<u64> = Uniform::from(0..100000000);
        let values = rng
            .clone()
            .sample_iter(binary_dist)
            .take(degree)
            .collect_vec();
        let values0 = rng.sample_iter(binary_dist).take(degree).collect_vec();
        let g0 = BigUint::from(5usize);
        let mut p1 = Poly::try_from_vec_u64(&params.rq_context, &values);
        p1.change_representation(Representation::Ntt);
        p1 *= &g0;

        let ct = sk.encrypt_poly(&p1);

        let unpacked_cts = evaluation.unpack(ct);
        izip!(unpacked_cts.iter(), values0.iter()).for_each(|(ct, v)| {
            let mut ex_p = Poly::try_from_vec_u64(&params.rq_context, &[*v]);
            ex_p.change_representation(Representation::Ntt);
            ex_p *= &g0;

            let p = sk.decrypt_trial(&ct.cts[0], &ct.cts[1]);

            let mut diff = &p - &ex_p;
            diff.change_representation(Representation::Ntt);

            let modulus = params.rq_context.rns.product.clone();

            Vec::<BigUint>::from(&diff).iter().for_each(|coeff| {
                dbg!(std::cmp::min(coeff.bits(), modulus.bits() - coeff.bits()));
            });
        });
    }

    #[test]
    fn unpack() {
        let mut rng = thread_rng();
        let degree = 64;
        let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let evaluation = Evaluation::build(&sk, &params);

        let binary_dist = Uniform::from(0..100);
        let values = rng.sample_iter(binary_dist).take(degree).collect_vec();
        let ct = sk.encrypt(&Plaintext::new(&params, &values));

        let unpacked_cts = evaluation.unpack(ct);

        let mut r_values = vec![];
        unpacked_cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            r_values.push(p.values[0]);
        });

        // error counts
        let mut count = 0;
        izip!(values.iter(), r_values.iter()).for_each(|(v, r_v)| {
            if v != r_v {
                count += 1;
            }
        });

        println!("{:?}", values);
        println!("{:?}", r_values);
        dbg!(count);
    }
}
