use crate::{
    bfv::{BfvCipherText, BfvParameters, GaliosKey, Plaintext, SecretKey},
    ksk::Ksk,
    poly::Modulus,
    rgsw::RgswCt,
    rq::{Poly, Representation, RqContext, Substitution},
    utils::ilog2,
};
use itertools::{enumerate, izip, Itertools};
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
}

impl Server {
    pub fn new(first_dim: usize, dim: usize, db: &Vec<Plaintext>) -> Self {
        Self {
            first_dim,
            dim,
            db: db.clone(),
        }
    }

    fn no_of_cts(&self, decomposition_base: usize) -> usize {
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
        decomposition_base: usize,
        query_ct: &BfvCipherText,
        eval_key: &Evaluation,
        //FIXEME: REMOVE
        sk: &SecretKey,
    ) -> Vec<BfvCipherText> {
        let mut cts = eval_key.unpack(query_ct.clone());
        let params = sk.params.clone();

        // truncate cts
        let no_cts = self.no_of_cts(decomposition_base);
        cts = cts[0..no_cts].to_vec();

        // cts.iter().for_each(|c| {
        //     println!("{:?} cs", sk.decrypt(c).values);
        // });

        let first_dim_bfvs = &cts[..(self.first_dim as usize)];

        let dimensions_queries: Vec<Vec<RgswCt>> = cts[(self.first_dim as usize)..]
            .chunks(self.dim as usize)
            .into_iter()
            .enumerate()
            .map(|(dim_index, dim_cts)| {
                dim_cts
                    .into_iter()
                    .enumerate()
                    .map(|(row_index, row_ct)| {
                        dbg!(row_index);
                        println!("{:?}", sk.measure_noise(&row_ct));

                        let rgsw_lower_half = params
                            .rq_context
                            .rns
                            .garner
                            .iter()
                            .map(|gi| row_ct * gi)
                            .collect_vec();
                        let rgsw_upper_half = rgsw_lower_half
                            .iter()
                            .map(|m_gi| {
                                let (c0, c1) = RgswCt::external_product(&eval_key.sk_rgsw, &m_gi);
                                BfvCipherText {
                                    params: params.clone(),
                                    cts: vec![c0, c1],
                                }
                            })
                            .collect();

                        let rgsw_lower_half_ksk =
                            Ksk::try_from_decomposed_rlwes(&params, &rgsw_lower_half);
                        let rgsw_upper_half_ksk =
                            Ksk::try_from_decomposed_rlwes(&params, &rgsw_upper_half);

                        RgswCt::try_from_ksks(&rgsw_upper_half_ksk, &rgsw_lower_half_ksk)
                    })
                    .collect()
            })
            .collect();

        // dbg!(dimensions_queries.len());
        // dbg!(dimensions_queries[0].len());

        // First dimension
        let mut db: Vec<BfvCipherText> = self
            .db
            .chunks(self.first_dim)
            .into_iter()
            .map(|pts| {
                let rows = izip!(first_dim_bfvs.iter(), pts.iter()).map(|(q, p)| q * p);

                let mut sum = BfvCipherText {
                    params: pts[0].params.clone(),
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
        dimensions_queries
            .iter()
            .enumerate()
            .for_each(|(dim_index, query_rgsws)| {
                let tmp = db
                    .chunks(self.dim)
                    .into_iter()
                    .enumerate()
                    .map(|(row_index, cts)| {
                        let rows = izip!(query_rgsws.iter(), cts.iter()).map(|(rgsw, bfv)| {
                            let (c0, c1) = RgswCt::external_product(rgsw, bfv);
                            let f = BfvCipherText {
                                params: bfv.params.clone(),
                                cts: vec![c0, c1],
                            };
                            // println!(
                            //     "dim: {}, row: {}, f: {:?}, bfv: {:?}",
                            //     dim_index,
                            //     row_index,
                            //     sk.decrypt(&f).values,
                            //     sk.decrypt(bfv).values
                            // );
                            f
                        });
                        let mut sum = BfvCipherText {
                            params: cts[0].params.clone(),
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
        ksk,
        rns::RnsContext,
        rq::{BitDecomposition, Representation},
        utils::generate_prime,
    };
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use rand::{thread_rng, Rng};
    use rand_distr::Uniform;

    #[test]
    fn query() {
        let mut rng = thread_rng();
        let degree = 64;
        let pt_moduli = generate_prime(50, 2 * 1048576, (1 << 50) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62, 62], degree).unwrap();
        // let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let degree = 64;

        // Pre-processing
        let db: Vec<Plaintext> = (0..8)
            .into_iter()
            .map(|i| {
                let v = [i, i, i, i];
                Plaintext::new(&params, &v.to_vec())
            })
            .collect();

        // Client side
        let sk = &SecretKey::generate(&params);
        let query = vec![3usize, 0];
        let client = Client::new(4, 2);
        let query_pt = client.encode(query, db.len(), &params, &sk);
        let query_pt = Plaintext::new(&params, &query_pt);
        let query_ct = sk.encrypt(&query_pt);
        let eval_key = Evaluation::build(&sk, &params);

        // Server side
        let server = Server::new(4, 2, &db);
        let res_ct = server.decode(
            params.rq_context.rns.moduli.len(),
            &query_ct,
            &eval_key,
            &sk,
        );

        // Client
        dbg!(res_ct.len());
        let res_ct = res_ct[0].clone();
        let res_pt = sk.decrypt(&res_ct);
        println!("{:?}", res_pt.values);
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
    fn unpack() {
        let mut rng = thread_rng();
        let degree = 64;
        let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        let pt_moduli: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));
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

        let evaluation = Evaluation::build(&sk, &params);

        let binary_dist = Uniform::from(0..2);
        let values = rng.sample_iter(binary_dist).take(60).collect();
        let pt = Plaintext::new(&params, &values);

        let ct = sk.encrypt(&pt);
        let d_values = sk.decrypt(&ct);

        let unpacked_cts = evaluation.unpack(ct);
        let mut r_values = vec![];
        unpacked_cts.iter().for_each(|c| {
            let p = sk.decrypt(c);
            r_values.push(p.values[0]);
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

    #[test]
    fn poly_test() {
        let mut rng = thread_rng();
        let degree = 64;
        // let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        // dbg!(&ct_moduli);
        let pt_moduli = generate_prime(60, 2 * 1048576, (1 << 60) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62, 62], degree).unwrap();
        // let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let mut one_poly = Poly::try_from_vec_u64(&params.rq_context, &[1u64]);
        one_poly.change_representation(Representation::Ntt);

        let b_ksk = Ksk::new(&sk, &one_poly);
    }

    #[test]
    fn trial() {
        let mut rng = thread_rng();
        let degree = 8;
        // let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        // dbg!(&ct_moduli);
        let pt_moduli = generate_prime(60, 2 * 1048576, (1 << 60) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62], degree).unwrap();
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let mut one_poly = Poly::try_from_vec_u64(&params.rq_context, &[1u64]);
        one_poly.change_representation(Representation::Ntt);
        let b_ksk = Ksk::new(&sk, &one_poly);
    }

    #[test]
    fn decomp() {
        let mut rng = thread_rng();
        let degree = 64;
        // let ct_moduli = BfvParameters::generate_moduli(&[50, 55, 55], degree).unwrap();
        // dbg!(&ct_moduli);
        let pt_moduli = generate_prime(60, 2 * 1048576, (1 << 60) - 1).unwrap();
        let ct_moduli = BfvParameters::generate_moduli(&[62, 62, 62, 62, 62], degree).unwrap();
        // let params = Arc::new(BfvParameters::new(degree, pt_moduli, ct_moduli, 10));
        let bit_decomp = BitDecomposition { base: 4, l: 8 };
        let params = Arc::new(BfvParameters::new(
            degree, pt_moduli, ct_moduli, 10, bit_decomp,
        ));

        let sk = SecretKey::generate(&params);

        let mut one_poly = Poly::try_from_vec_u64(&params.rq_context, &[1u64]);
        one_poly.change_representation(Representation::Ntt);

        let b_ksk = Ksk::new(&sk, &one_poly);

        let e_decomp = izip!(b_ksk.c0.iter(), b_ksk.c1.iter())
            .map(|(c0, c1)| BfvCipherText {
                params: params.clone(),
                cts: vec![c0.clone(), c1.clone()],
            })
            .collect_vec();
        let decomp = e_decomp
            .iter()
            .map(|f| sk.decrypt(f).values[0])
            .map(|bi| {
                sk.encrypt(&Plaintext {
                    params: params.clone(),
                    values: [bi].to_vec().into_boxed_slice(),
                })
            })
            .collect_vec();

        izip!(e_decomp.iter(), decomp.iter()).for_each(|(e, r)| {
            println!("{:?}", sk.decrypt(e).values);
            println!("{:?}", sk.decrypt(r).values);
            println!();
        });

        let mut sk_poly = Poly::try_from_vec_i64(&sk.params.rq_context, &sk.coeffs);
        sk_poly.change_representation(Representation::Ntt);
        let sk_rgsw = RgswCt::encrypt_poly(&sk, &sk_poly);

        izip!(
            e_decomp.iter(),
            decomp.iter(),
            sk.params.rq_context.rns.garner.iter()
        )
        .for_each(|(e, r, gi)| {
            let mut p = Poly::try_from_vec_u64(&sk.params.rq_context, &[1u64]);
            p.change_representation(Representation::Ntt);
            dbg!(&p.representation);
            // let zx = BigUint::from_u64(5).unwrap();
            p *= gi;
            println!(
                "ex: {:?}, re: {:?}",
                measure_noise_tmp(&sk.params, p.clone(), &sk, e),
                measure_noise_tmp(&sk.params, p, &sk, r)
            );
        });

        // {
        //     let mut p = BfvPlaintext::new(&sk.params, &vec![5u64].to_vec()).to_poly();
        //     p.change_representation(Representation::Ntt);
        //     dbg!(&p.representation);
        //     let zx = BigUint::from_u64(913091389012891028).unwrap();
        //     p *= &zx;
        //     // p.change_representation(Representation::PowerBasis);
        //     // println!("{:?}", p);

        //     let ct = sk.encrypt_poly(&p);
        //     let noise = measure_noise_tmp(&sk.params, p.clone(), &sk, &ct);
        //     dbg!(noise);
        //     // dbg!(sk.decrypt(&ct).values);
        // }

        // izip!(
        //     e_decomp.iter(),
        //     decomp.iter(),
        //     sk.params.rq_context.rns.garner.iter()
        // )
        // .for_each(|(e, r, gi)| {
        //     let e = RgswCt::external_product(&sk_rgsw, e);
        //     // let e = (e.cts[0].clone(), e.cts[1].clone());
        //     let mut em = &e.1 * &sk_poly;
        //     em += &e.0;

        //     // println!(
        //     //     "e extern {:?}",
        //     //     sk.decrypt(&BfvCipherText {
        //     //         params: params.clone(),
        //     //         cts: vec![e.0.clone(), e.1.clone()]
        //     //     })
        //     //     .values
        //     // );

        //     let r = RgswCt::external_product(&sk_rgsw, r);
        //     // let r = (r.cts[0].clone(), r.cts[1].clone());
        //     let mut rm = &r.1 * &sk_poly;
        //     rm += &r.0;

        //     // println!(
        //     //     "r extern {:?}",
        //     //     sk.decrypt(&BfvCipherText {
        //     //         params: params.clone(),
        //     //         cts: vec![r.0.clone(), r.1.clone()]
        //     //     })
        //     //     .values
        //     // );

        //     let mut sub = &em - &rm;
        //     sub.change_representation(Representation::PowerBasis);

        //     let rns = RnsContext::new(&params.ciphertext_moduli.clone());
        //     Vec::<BigUint>::from(&sub).iter().for_each(|coeff| {
        //         dbg!(std::cmp::min(coeff.bits(), (rns.modulus() - coeff).bits()));
        //     });

        //     println!();
        // });
    }

    #[test]
    fn test_bit_vector() {
        let m = BigUint::from_u128(u128::MAX - 1).unwrap();
        let bits = bit_vector(128, &m);

        let value = bits
            .iter()
            .rev()
            .fold(BigUint::zero(), |acc, bit| (acc * 2usize) + (*bit as u8));

        assert_eq!(value, m);
    }

    #[test]
    fn test_decompose_bits() {
        let value = BigUint::from_u128(u128::MAX - 1).unwrap();
        let decomposed = decompose_bits(&value, 128, 4);
        let base = BigUint::one() << 4;
        let recovered = decomposed
            .iter()
            .rev()
            .fold(BigUint::zero(), |acc, b| (acc * &base) + *b);

        assert_eq!(value, recovered);
    }
}

pub fn bit_vector(bit_size: usize, m: &BigUint) -> Vec<bool> {
    let mut bit_vec = vec![false; bit_size];

    let mut value = m.clone();
    let mask = BigUint::from_usize(1).unwrap();
    let mut index = 0;
    while !value.is_zero() {
        let bit = &value & &mask;

        if bit.is_one() {
            bit_vec[index] = true;
        } else {
            bit_vec[index] = false;
        }

        value >>= 1;
        index += 1;
    }

    bit_vec
}

pub fn decompose_bits(m: &BigUint, parent_bits: usize, decomp_bits: usize) -> Vec<u64> {
    let bit_vector = bit_vector(parent_bits, m);

    let decomp_m = bit_vector
        .chunks(decomp_bits)
        .into_iter()
        .map(|chunk| {
            chunk
                .iter()
                .rev()
                .fold(0u64, |acc, b| acc * 2 + (*b as u64))
        })
        .collect();

    decomp_m
}

pub fn measure_noise_tmp(
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
