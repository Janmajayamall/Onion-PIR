use super::bfv::BfvCipherText;
use crate::{
    bfv::{BfvPlaintext, BfvPublicKey, BfvSecretKey},
    poly::Poly,
};

/// Key switching key
///
/// Params \beta & l
/// Note: Least significant first
/// [
///     RLWE_s'( q/beta^l * s)
///     RLWE_s'( q/beta^{l-1} * s)
///     .
///     .   
///     .
///     RLWE_s'( q/beta^1 * s)
/// ]
/// `RLWE_s'` means encrypted under `s'`,
/// where `s'` is the new key
struct Ksk {
    cts: Vec<BfvCipherText>,
    new_pk: BfvPublicKey,
    beta: u64,
    l: u64,
}

impl Ksk {
    pub fn new(curr_sk: &BfvSecretKey, new_sk: &BfvSecretKey, beta: u64) -> Self {
        let q_ref = curr_sk.params.poly_ctx.moduli.q;
        let l = ((q_ref as f64).log2() / (beta as f64).log2()).floor() as u64;

        let new_pk = BfvPublicKey::new(new_sk);
        let cts = (l..0).into_iter().map(|index| {
            // encrypt under new_sk
            let pt = &curr_sk * (q_ref / beta.pow(index as u32));
            let pt = BfvPlaintext::new(&new_sk.params, pt.coeffs.into());
            new_pk.encrypt(&pt)
        });

        Self {
            cts: cts.collect(),
            new_pk,
            beta,
            l,
        }
    }

    pub fn key_switch(ksk: Ksk, ct: BfvCipherText) -> BfvCipherText {
        let a_decomposed = ct.c[1].decompose(ksk.beta);

        debug_assert!(a_decomposed.len() == ksk.cts.len());

        // RLWE encryption of `a * curr_sk` under new_sk
        let a_curr_s = a_decomposed
            .iter()
            .zip(ksk.cts.iter())
            .map(|(a_pt, rlwe_ct)| rlwe_ct.multiply_pt_poly(a_pt))
            .fold(
                BfvCipherText::new(
                    &ksk.new_pk.params,
                    vec![
                        Poly::zero(&ksk.new_pk.params.poly_ctx),
                        Poly::zero(&ksk.new_pk.params.poly_ctx),
                    ],
                ),
                |acc, ct| BfvCipherText::add_ciphertexts(&acc, &ct),
            );

        // RLWE encryption of `b` under new_sk
        let b_under_new_sk = ksk.new_pk.encrypt(&BfvPlaintext::new(
            &ksk.new_pk.params,
            ct.c[0].coeffs.into(),
        ));

        // (/delta * M) + E_old = B - (A * S_curr)
        // Thus homomorphic computation of `B - (A * S_curr)` under new_sk
        // returns `(/delta * M) + E_old` under new_sk.
        //
        // Note that
        BfvCipherText::add_ciphertexts(&b_under_new_sk, &(&a_curr_s.negate()))
    }
}

/// Structure:
/// Params: \beta & l
/// Note: Least significant first
/// [
///     [
///         RLWE( q/beta^l * -s * m)
///         RLWE( q/beta^{l-1} * -s * m)
///         .
///         .
///         .
///         RLWE( q/beta^1 * -s * m)
///     ]
///     [
///         RLWE( q/beta^l * m)
///         RLWE( q/beta^{l-1} * m)
///         .
///         .
///         .
///         RLWE( q/beta^{1} * m)
///     ]
/// ]
/// E (2l X 2)
struct Rgsw {
    cts: Vec<Vec<BfvCipherText>>,
    beta: u64,
    l: usize,
}

impl Rgsw {
    fn external_product(rgsw: Rgsw, rlwe: BfvCipherText) -> BfvCipherText {
        // Steps (C[0] is B, C[1] is A)
        // (1) Decompose B and A of rlwe into B` and A`
        // (2) Perform <B`, RGWS[1]> and <A`, RGSW[0]>

        let b_decomposed = rlwe.c[0].decompose(rgsw.beta);
        let a_decomposed = rlwe.c[1].decompose(rgsw.beta);

        debug_assert!(b_decomposed.len() == rgsw.cts[1].len());
        debug_assert!(a_decomposed.len() == rgsw.cts[0].len());

        let bfv_ref = rgsw.cts[0][0].params.clone();

        let bm = b_decomposed
            .iter()
            .zip(rgsw.cts[1].iter())
            .map(|(pt_poly, rlwe)| rlwe.multiply_pt_poly(pt_poly))
            .fold(
                BfvCipherText::new(
                    &bfv_ref,
                    vec![Poly::zero(&bfv_ref.poly_ctx), Poly::zero(&bfv_ref.poly_ctx)],
                ),
                |acc, poly| BfvCipherText::add_ciphertexts(&acc, &poly),
            );

        let n_asm = a_decomposed
            .iter()
            .zip(rgsw.cts[0].iter())
            .map(|(pt_poly, rlwe)| rlwe.multiply_pt_poly(pt_poly))
            .fold(
                BfvCipherText::new(
                    &bfv_ref,
                    vec![Poly::zero(&bfv_ref.poly_ctx), Poly::zero(&bfv_ref.poly_ctx)],
                ),
                |acc, poly| BfvCipherText::add_ciphertexts(&acc, &poly),
            );

        BfvCipherText::add_ciphertexts(&bm, &n_asm)
    }
}
