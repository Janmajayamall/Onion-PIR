use super::bfv::BfvCipherText;
use crate::{
    bfv::{BfvPlaintext, BfvPublicKey, BfvSecretKey},
    poly::Poly,
};

/// Key switching key
///
/// Params ß & l
/// Note: Least significant first
/// [
///     RLWE_s'( q/ß^l * s)
///     RLWE_s'( q/ß^{l-1} * s)
///     .
///     .   
///     .
///     RLWE_s'( q/ß^1 * s)
/// ]
/// `RLWE_s'` means encrypted under `s'`,
/// where `s'` is the new key
pub struct Ksk {
    pub cts: Vec<BfvCipherText>,
    new_pk: BfvPublicKey,
    beta: u64,
    l: u64,
    is_subs: bool,
    subs_k: Option<u64>,
}

impl Ksk {
    pub fn new(
        curr_sk: &BfvSecretKey,
        new_sk: &BfvSecretKey,
        beta: u64,
        is_subs: bool,
        subs_k: Option<u64>,
    ) -> Self {
        assert!(is_subs && subs_k.is_some() || subs_k.is_none());

        let q_ref = curr_sk.params.poly_ctx.moduli.q;
        let l = ((q_ref as f64).log2() / (beta as f64).log2()).floor() as u64;

        let new_pk = BfvPublicKey::new(new_sk);
        let cts = (l..0).into_iter().map(|index| {
            // encrypt under new_sk
            let pt = &curr_sk.poly * (q_ref / beta.pow(index as u32));
            let pt = BfvPlaintext::new(&new_sk.params, pt.coeffs.into());
            new_pk.encrypt(&pt)
        });

        Self {
            cts: cts.collect(),
            new_pk,
            beta,
            l,
            is_subs,
            subs_k,
        }
    }

    pub fn key_switch(ksk: &Ksk, ct: BfvCipherText) -> BfvCipherText {
        let a_decomposed = ct.c[1].decompose(ksk.beta);

        debug_assert!(a_decomposed.len() == ksk.cts.len());

        // RLWE encryption of `a * curr_sk` under new_sk
        let mut a_curr_s = a_decomposed
            .iter()
            .zip(ksk.cts.iter())
            .map(|(a_pt, rlwe_ct)| BfvCipherText::multiply_pt_poly(rlwe_ct, a_pt))
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
            ct.c[0].coeffs.clone(),
        ));

        // (/delta * M) + E_old = B - (A * S_curr)
        // Thus homomorphic computation of `B - (A * S_curr)` under new_sk
        // returns `(/delta * M) + E_old` under new_sk.
        BfvCipherText::sub_ciphertexts(&b_under_new_sk, &a_curr_s)
    }

    /// Subs(•,k) converts `RLWE(Σ b_i • X^i)`
    /// to `RLWE(Σ b_i • (X^i)^k)`.
    /// The operation uses key switching
    /// along with the fact that:
    /// If `ct(X) = (c0(X), c1(X))` for `M(X)`, where
    /// `M(X) = Σ b_i • X^i` under S(X) as secret key
    /// then
    /// `ct(X^k) = (c0(X^k), c1(X^k))` for `M(X^k)`
    /// (i.e. Σ b_i • (X^i)^k) under S(X^k) as secret key.
    ///
    /// So the trick is as follows:
    /// For `Subs(•,k)`, take ct(X) transform it to ct(X^k).
    /// Create new `Ksk` with params `curr_sk=S(X^k)` and
    /// `new_key=S(X)` and then perform key switch operation
    /// so that you obtain M(X^k) encrypt under `new_key` that is
    /// `S(X)`
    ///
    /// Also note that since we are operating on polynomial
    /// integers modulus some cyclotomic polynomial `X^N + 1`,
    /// following holds whenever k = N + 1
    /// Subs(•, N + 1) transforms
    /// `RLWE(Σ b_i • X^i)`
    /// => `RLWE(Σ b_i • (X^i)^(N+1))`
    /// => `RLWE(Σ b_i • (X^(Ni+i))` ---- (1)
    /// Since X^N + 1 = 0 mod (X^N + 1)
    /// => X^N = -1 mod (X^N + 1)
    /// => (X^N)^i = X^(N*i) = -1^i mod (X^N + 1) ---- (2)
    /// Therefore following from (1) and (2)
    /// we can write (1) as following
    /// => `RLWE(Σ b_i • -1^(i) • X^(i))`
    ///
    /// `ksk` is key switching key corresponding to
    /// Subs(•,k).  
    pub fn substitute(ksk: &Ksk, ct: &BfvCipherText) -> BfvCipherText {
        assert!(ksk.is_subs && ksk.subs_k.is_some());

        // ct(X) -> ct(X^k)
        let ct_k = BfvCipherText::new(
            &ct.params,
            ct.c.iter()
                .map(|poly| poly.shift_powers(ksk.subs_k.unwrap()))
                .collect(),
        );

        Ksk::key_switch(ksk, ct_k)
    }

    /// Generates key switching key for
    /// Subs(•, k) operation
    ///
    /// Takes in secret key `sk` and returns
    /// key switching key for necessary ops.
    pub fn gen_substitution_key(sk: &BfvSecretKey, k: u64, beta: u64) -> Self {
        // S(X^k)
        let sk_k_poly = sk.poly.shift_powers(k);
        let sk_k = BfvSecretKey::new_with_poly(&sk.params, sk_k_poly);
        Ksk::new(&sk_k, sk, beta, true, Some(k))
    }
}

/// Structure:
/// Params: \ß & l
/// Note: Least significant first
/// [
///     [
///         RLWE( q/ß^l * -s * m)
///         RLWE( q/ß^{l-1} * -s * m)
///         .
///         .
///         .
///         RLWE( q/ß^1 * -s * m)
///     ]
///     [
///         RLWE( q/ß^l * m)
///         RLWE( q/ß^{l-1} * m)
///         .
///         .
///         .
///         RLWE( q/ß^{1} * m)
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
            .map(|(pt_poly, rlwe)| BfvCipherText::multiply_pt_poly(&rlwe, pt_poly))
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
            .map(|(pt_poly, rlwe)| BfvCipherText::multiply_pt_poly(&rlwe, pt_poly))
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
