use std::{sync::Arc, vec};

use rand_distr::num_traits::Pow;

use crate::{
    bfv::{BfvCipherText, BfvParameters, BfvPlaintext, BfvPublicKey},
    poly::Poly,
    rgsw::Ksk,
    utils::ilog2,
};

struct QueryParams {
    first_dim: u64,
    normal_dim: u64,
}

fn build_query(
    query_params: QueryParams,
    mut db_dim: u64,
    mut query_index: u64,
    bfv_pk: &BfvPublicKey,
) -> BfvCipherText {
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

    // encode dimension specific query vector bits into bfv ciphertext
    let pt = BfvPlaintext::new(
        &bfv_pk.params,
        expanded_dims.iter().flatten().copied().collect(),
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
/// `k = n/2^(i-1) + 1`
///                                    [x0, x1, x2, x3, x4, x5, x6, x7] depth of tree: logN
///                                      / i = 1                                       \
///                                     / C + Subs(•, k)                                \ C - Subs(•, k)          
///                         2^k*[x0, x2, x4, x6], n/2^(k-1) + 1         2^k*[x1, x3, x5, x7] * x^-k => 2*[x0, x2, x4, x6, x8]
///                                    / i = 2              \
///                                   / C + Subs(•, k)       \ C - Subs(•, k)          
///            Subs(2^k*[x0, x4], n/2^(k-1) + 1)          2^k*[x2, x6] * -k
///                                  / i = 3                    \
///                                 / C + Subs(•, k)             \ C - Subs(•, k)          
///                             2^k*[x0]                        2^k*[x4] * x^-k
///
/// Ref - Algorithm 3 & 4 of https://eprint.iacr.org/2019/736.pdf
fn resolve_query(query_ct: &BfvCipherText, ksks: Vec<Ksk>) {
    let mut enc_bits: Vec<BfvCipherText> = vec![query_ct.clone()];
    let n = query_ct.params.poly_ctx.degree;
    let logn = ilog2(n);

    // TODO: Figure out way to check the order of `ksks`
    // such that their `subs_k` match: n, n/2, n/4...n/2^logn - 1
    assert!(ksks.len() == logn - 1);

    for i in 1..logn {
        let k = (n as usize / 2.pow(i - 1) as usize) + 1;
        let mut curr_branch: Vec<BfvCipherText> = vec![];
        for j in 0..enc_bits.len() {
            // c_even = c + Subs(•, k)
            curr_branch.push(BfvCipherText::add_ciphertexts(
                &enc_bits[j],
                &Ksk::substitute(&ksks[i], &enc_bits[j]),
            ));

            // c_odd = c - Subs(•, k)
            curr_branch.push(BfvCipherText::add_ciphertexts(
                &enc_bits[j],
                &Ksk::substitute(&ksks[i], &enc_bits[j]),
            ));
        }
    }
}
