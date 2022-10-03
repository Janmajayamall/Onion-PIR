use rand_distr::num_traits::Pow;
use std::{sync::Arc, vec};

use crate::{
    bfv::{BfvCipherText, BfvParameters, BfvPlaintext, BfvPublicKey},
    poly::{Context, Modulus, Poly},
    rgsw::Ksk,
    utils::{ilog2, mod_inverse},
};

pub struct QueryParams {
    // 128 by default
    pub first_dim: u64,
    // 4 by default
    pub normal_dim: u64,
    pub gadget_beta: u64,
}

fn build_query(
    query_params: QueryParams,
    mut db_dim: u64,
    mut query_index: u64,
    bfv_pk: &BfvPublicKey,
) -> BfvCipherText {
    // `gadget_beta` is for RLWE to RGSW trick. Make sure
    // it is a power of 2 and smaller than `t`.
    assert!(query_params.gadget_beta.is_power_of_two());

    let mut no_of_dims = 0;
    let mut first = true;
    while db_dim != 0 {
        if first {
            db_dim /= query_params.first_dim;
            first = false;
        } else {
            db_dim /= query_params.normal_dim;
        }
        no_of_dims += 1;
    }

    let mut dims = Vec::new();
    (1..no_of_dims + 1).into_iter().for_each(|dim| {
        if dim == 1 {
            dims.push(query_index % query_params.first_dim);
            query_index /= query_params.first_dim;
        } else {
            dims.push(query_index % query_params.normal_dim);
            query_index /= query_params.normal_dim
        }
    });

    // Imagine the db vector as a hypercube of `d` dimensions:
    // [d0, d1, d2, d3, d4...d{d-1}].
    // To reach a specific query index you iteratively
    // divide current db dimension with next hypercube dimension.
    // At each iteration (i.e. for current dimension) your dimension specific query
    // index is equal to last query index `mod` current db dimension.
    // You next query index is equal to last query index / current db dimension.
    //
    // Once you have obtained vector containing dimension specific query
    // indexes. You need to expand them into a vec to 0s and 1s for bfv ciphertext
    // encoding. For example, consider following as dimension specific query indexes
    // [8, 2, 3].
    // You expand them to
    // [
    //     [
    //         // first_dim len
    //         0,0,0,0,0,0,0,1,0,0,...,0
    //     ],
    //     [
    //         // normal_dim len
    //         0,0,1,0
    //     ],
    //     [
    //         // normal_dim len
    //         0,0,0,1
    //     ]
    // ]
    // and encode them into Bfv ciphertext. This means simply flat
    // map all values as coefficients of the Bfv plaintext
    // polynomial.
    let mut expanded_dims = Vec::new();
    dims.iter().enumerate().for_each(|(index, value)| {
        let mut expanded_dim = Vec::<u64>::new();
        (0..(if index == 0 {
            query_params.first_dim
        } else {
            query_params.normal_dim
        }))
            .into_iter()
            .for_each(|i| {
                if i == *value {
                    expanded_dim.push(1);
                } else {
                    expanded_dim.push(0);
                }
            });
        expanded_dims.push(expanded_dim);
    });

    let bfv_params = bfv_pk.params.clone();

    // generate coefficients
    let l = ((bfv_params.q as f64).log2() / (query_params.gadget_beta as f64).log2()) as u64;
    let pt_space = Modulus { q: bfv_params.t };
    let beta_inv = mod_inverse(query_params.gadget_beta, bfv_params.t).unwrap();
    let mut q_beta_inv = pt_space.convert(bfv_params.q);
    let coeffs: Vec<_> = (0..l)
        .map(|_| {
            q_beta_inv = pt_space.mul_mod(q_beta_inv, beta_inv);
            q_beta_inv
        })
        .collect();

    // We need to further expand each value to `l` values in dimension vector of every dimension, except first
    // to perform RLWE to RGSW conversion trick
    // TODO: Explain what the trick is
    let mut normal_dims_3d: Vec<Vec<Vec<u64>>> = vec![];

    // handle first dim separately
    let d: Vec<_> = expanded_dims[0].iter().map(|bit| vec![*bit]).collect();
    normal_dims_3d.push(d);

    expanded_dims.iter().skip(1).for_each(|dim_vec| {
        let dim_vec_2d: Vec<Vec<_>> = dim_vec
            .iter()
            .map(|bit_value| {
                let f: Vec<_> = coeffs.iter().map(|coeff| bit_value * coeff).collect();
                f
            })
            .collect();
        normal_dims_3d.push(dim_vec_2d)
    });

    // encode dimension specific query vector bits into bfv ciphertext
    let pt = BfvPlaintext::new(
        &bfv_pk.params,
        normal_dims_3d
            .iter()
            .flatten()
            .into_iter()
            .flatten()
            .into_iter()
            .copied()
            .collect(),
    );

    // So now we have a Bfv ciphertext that encrypts
    // a polynomial that consists of bits as its
    // coefficients.
    // In order to determine indexes the client is
    // interested in hypercube, the server needs to
    // convert ciphertext encrypting polynomial with
    // bits encoded as coefficients to ciphertexts
    // encrypting individual bits. This is achieved
    // using `Subs(•,k)` operation.
    bfv_pk.encrypt(&pt)
}

///
/// Let's consider the following bit vector
/// [0, 1, 0, 1, 1, 0...]
///
/// We encode the bit vector into a `BfvPlaintext` as following
/// `0•X^0 + 1•X^1 + 0•X^2 + 1•X^3 + ... + Nth_bit•X^{N-1}
/// and then encrypt it under secret key `SK`. Let's call the
/// encrypted ciphertext `C`.
///
/// We need to recover individual bits encrypted under `SK`
/// from the ciphertext `C`. We achieve this using the following
/// way:
///
/// Assume that our poly operations are modulo cyclotomic
/// polynomial `X^16 + 1`. (i.e. N = 16)
///
/// Note: that `X^N = -1` and `(X^i)^N+1 = X^(Ni+i) = X^Ni * X^i`
/// Therefore, `(X^i)^N+1 = -1^i * X^i`
///
/// We will use the statement above to separate encryption of odd
/// terms and even terms in `C` into two different encryptions.
///
/// `C_even = C + Subs(C, N+1)` --> Encryption of only even terms
/// `C_odd = C - Subs(C, N+1)` --> Encryption of only odd terms
///
/// Note: that `Subs(C, N+1)` transforms `Enc(Σ b_i • X^i)` to
/// `Enc(Σ b_i • (X^i)^(N+1))`. Also, `Σ b_i • (X^i)^(N+1)` is equal
/// to `Σ b_i • (X^i)` except that all its odd position are -ve.
/// Therefore,
/// `Enc(Σ b_i • X^i)` + `Enc(Σ b_i • (X^i)^(N+1))` = Ct_even
/// AND
/// `Enc(Σ b_i • X^i)` - `Enc(Σ b_i • (X^i)^(N+1))` = Ct_odd
///
/// Also notice that in case of `Ct_odd` the powers are odd, which
/// wouldn't isn't useful in subsequent steps. Thus we change Ct_odd
/// `Ct_odd = Ct_odd * X^-1`
///
/// Therefore now we have
/// `Ct_even = C + Subs(C, N+1) = Enc(Σ 2b_2i • X^2i)`
/// `Ct_odd = (C - Subs(C, N+1)) • X^-1 = Enc(Σ 2b_2i+1 • X^2i)`
///
/// Now to obtain `Enc(b_i)`s of bits we perform the same operation
/// recursively. After `log(N)` branches we will have `Enc(N * b_i)`
/// as N leaf nodes.
/// `k = n/2^(i) + 1`
///                                    [x0, x1, x2, x3, x4, x5, x6, x7] depth of tree: logN
///                                      / i = 0                                       \
///                                     / C + Subs(•, k)                                \ C - Subs(•, k)          
///                                  2^(i+1)*[x0, x2, x4, x6]                     2^(i+1)*[x1, x3, x5, x7] * x^-(2^i) => 2^i*[x0, x2, x4, x6, x8]
///                                    / i = 1              \
///                                   / C + Subs(•, k)       \ C - Subs(•, k)          
///                             2^(i+1)*[x0, x4]             2^(i+1)*[x2, x6] * x^-(2^i)
///                             / i = 2        \
///                            / C + Subs(•, k) \ C - Subs(•, k)          
///                     2^(i+1)*[x0]         2^(i+1)*[x4] * x^-(2^i)
///
/// Ref - Algorithm 3 & 4 of https://eprint.iacr.org/2019/736.pdf
fn resolve_query(query_ct: &BfvCipherText, ksks: Vec<Ksk>) -> Vec<BfvCipherText> {
    let mut enc_bits: Vec<BfvCipherText> = vec![query_ct.clone()];
    let n = query_ct.params.poly_ctx.degree;
    let logn = ilog2(n);

    // TODO: Figure out way to check the order of `ksks`
    // such that their `subs_k` match: n, n/2, n/4...n/2^logn - 1
    assert!(ksks.len() == logn - 1);

    // Should be in `Rt` since all X^-(2^i) polys are plaintext
    let rt_ctx = Arc::new(Context::new(
        Modulus {
            q: query_ct.params.t,
        },
        query_ct.params.n,
    ));
    let x_polys = gen_x_pow_polys(&rt_ctx, n);

    for i in 0..logn {
        let k = (n as usize / 2.pow(i - 1) as usize) + 1;
        // make sure that `Ksk` is for subs ops at current branch `i`
        assert!(k as u64 == ksks[i].subs_k.unwrap());

        let mut curr_branch: Vec<BfvCipherText> = vec![];
        for j in 0..enc_bits.len() {
            // c_even = c + Subs(•, k)
            curr_branch.push(BfvCipherText::add_ciphertexts(
                &enc_bits[j],
                &Ksk::substitute(&ksks[i], &enc_bits[j]),
            ));

            // c_odd = c - Subs(•, k)
            // c_odd = c_odd * X^-(2^i)
            let c_odd = BfvCipherText::multiply_pt_poly(
                &BfvCipherText::sub_ciphertexts(
                    &enc_bits[j],
                    &Ksk::substitute(&ksks[i], &enc_bits[j]),
                ),
                &x_polys[i],
            );
            curr_branch.push(c_odd);
        }

        enc_bits = curr_branch;
    }

    assert!(enc_bits.len() == query_ct.params.n);

    // Note that enc_bits = [RLWE(n * b0), RLWE(n * b1), ...RLWE(n * b{n-1})].
    // Thus we multiply all RLWEs with inverse of `n` in pt modulus `t` to get
    // RLWE(bi)
    let n_inv = mod_inverse(query_ct.params.n as u64, query_ct.params.t).unwrap();
    enc_bits
        .iter()
        // change this to inverse
        .map(|ct| BfvCipherText::multiply_constant(ct, n_inv))
        .collect()

    // Now we structure cts of bits into encrypted query vector for every dimension.
    // Note that query vector of `first_dim` consists of RLWE cts.
    // For the rest of the dimensions we have vector of RGSW cts, thus we need
    // to convert RLWE cts to RGSW cts.
}

pub fn gen_x_pow_polys(
    // x^-(2^0) ... x^-(2^(logn))
    poly_ctx: &Arc<Context>,
    n: usize,
) -> Vec<Poly> {
    // X^-k = -X^(N-k)
    (0..ilog2(n))
        .into_iter()
        .map(|i| {
            let mut poly = Poly::zero(poly_ctx);
            poly.coeffs[1] = 1; // X
            poly.shift_powers((poly_ctx.degree - 2.pow(i) as usize) as u64); // X^(N-k)
            poly = -poly; // -X^(N-k)
            poly
        })
        .collect()
}
