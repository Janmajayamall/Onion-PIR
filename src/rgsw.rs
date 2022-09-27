use std::task::Poll;

use crate::poly::Poly;

use super::bfv::BfvCipherText;

/// Structure:
/// Params: \beta & l
/// [
///     [
///         RLWE( q/beta^1 * -s * m)
///         RLWE( q/beta^2 * -s * m)
///         .
///         .
///         .
///         RLWE( q/beta^l * -s * m)
///     ]
///     [
///         RLWE( q/beta^1 * m)
///         RLWE( q/beta^2 * m)
///         .
///         .
///         .
///         RLWE( q/beta^l * m)
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
