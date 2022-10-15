use std::sync::Arc;

use crate::poly::{Context, Modulus};

use super::poly::Poly;

#[derive(Debug, PartialEq)]
pub struct BfvParameters {
    pub t: u64,
    pub q: u64,
    pub n: usize,
    pub st_dev: f64,
    pub poly_ctx: Arc<Context>,
    pub pt_poly_ctx: Arc<Context>,
}

impl Default for BfvParameters {
    fn default() -> Self {
        let t = 4;
        let q: u64 = 65536;
        let degree: usize = 4;
        let ctx = Context::new(Modulus { q }, degree);
        let pt_ctx = Context::new(Modulus { q: t }, degree);
        Self {
            t,
            q,
            n: degree,
            st_dev: 3.2,
            poly_ctx: Arc::new(ctx),
            pt_poly_ctx: Arc::new(pt_ctx),
        }
    }
}

#[derive(Debug)]
pub struct BfvSecretKey {
    pub params: Arc<BfvParameters>,
    pub poly: Poly,
}

impl BfvSecretKey {
    pub fn new(params: &Arc<BfvParameters>) -> Self {
        let mut sk = Poly::random_gaussian(
            &Arc::new(Context::new(Modulus { q: 2 }, params.n)),
            params.st_dev,
        );
        sk.switch_context(&params.poly_ctx);
        Self {
            params: params.clone(),
            poly: sk,
        }
    }

    pub fn new_with_poly(params: &Arc<BfvParameters>, mut poly: Poly) -> Self {
        poly.switch_context(&params.poly_ctx);
        Self {
            params: params.clone(),
            poly,
        }
    }

    pub fn decrypt(&self, ct: &BfvCipherText) -> Vec<u64> {
        let mut poly = &ct.c[1] * &self.poly + ct.c[0].clone();
        dbg!(&poly.coeffs);
        poly.coeffs
            .iter_mut()
            .for_each(|c| *c = ((*c * self.params.t) as f64 / self.params.q as f64).round() as u64);

        poly.coeffs
    }
}

#[derive(Debug)]
pub struct BfvPublicKey {
    pub params: Arc<BfvParameters>,
    pub ciphertext: BfvCipherText,
}

impl BfvPublicKey {
    pub fn new(sk: &BfvSecretKey) -> Self {
        let a = Poly::random_uniform(&sk.params.poly_ctx);
        let e = Poly::random_gaussian(&sk.params.poly_ctx, sk.params.st_dev);

        Self {
            params: sk.params.clone(),
            ciphertext: BfvCipherText::new(&sk.params, vec![-(&a * &sk.poly + e), a]),
        }
    }

    pub fn encrypt(&self, plaintext: &Vec<u64>) -> BfvCipherText {
        let mut u = Poly::random_gaussian(
            &Arc::new(Context::new(Modulus { q: 2 }, self.params.poly_ctx.degree)),
            self.params.st_dev,
        );
        u.switch_context(&self.params.poly_ctx);

        let e1 = Poly::random_gaussian(&self.params.poly_ctx, self.params.st_dev);
        let e2 = Poly::random_gaussian(&self.params.poly_ctx, self.params.st_dev);

        let delta = (self.params.q as f64 / self.params.t as f64).floor() as u64;
        let mut m = Poly::zero(&self.params.poly_ctx);
        m.coeffs
            .iter_mut()
            .zip(plaintext.iter())
            .for_each(|(sv, v)| *sv = (v * delta) % self.params.q);

        // p0 * u + e1 + delta_m
        let c0 = &self.ciphertext.c[0] + &u + e1 + m;
        // p1 * u + e2
        let c1 = &self.ciphertext.c[1] + &u + e2;

        BfvCipherText::new(&self.params, vec![c0, c1])
    }
}

#[derive(Debug)]
pub struct BfvPlaintext {
    pub params: Arc<BfvParameters>,
    pub values: Box<[u64]>,
    pub poly: Poly,
}

impl BfvPlaintext {
    pub fn new(params: &Arc<BfvParameters>, pt: Vec<u64>) -> Self {
        let mut poly = Poly::zero(&params.pt_poly_ctx);
        poly.coeffs
            .iter_mut()
            .zip(pt.iter().map(|v| *v % params.t))
            .for_each(|(c, v)| *c = v);

        BfvPlaintext {
            params: params.clone(),
            values: pt.into_boxed_slice(),
            poly,
        }
    }
}

#[derive(Debug, Clone)]
pub struct BfvCipherText {
    pub params: Arc<BfvParameters>,
    pub c: Vec<Poly>,
}

impl BfvCipherText {
    pub fn new(params: &Arc<BfvParameters>, c: Vec<Poly>) -> Self {
        Self {
            params: params.clone(),
            c,
        }
    }

    pub fn add_ciphertexts(ct0: &BfvCipherText, ct1: &BfvCipherText) -> Self {
        debug_assert!(ct0.params == ct1.params);
        let c0 = &ct0.c[0] + &ct1.c[0];
        let c1 = &ct0.c[1] + &ct1.c[1];
        BfvCipherText::new(&ct0.params, vec![c0, c1])
    }

    pub fn sub_ciphertexts(ct0: &BfvCipherText, ct1: &BfvCipherText) -> Self {
        debug_assert!(ct0.params == ct1.params);
        let ct1 = BfvCipherText::negate(ct1);
        BfvCipherText::add_ciphertexts(ct0, &ct1)
    }

    pub fn negate(ct: &BfvCipherText) -> BfvCipherText {
        // to negate: -1 * Ct
        // -1 mod t == t-1 mod t
        // TODO: Handle such ops in a better way
        BfvCipherText::multiply_constant(ct, ct.params.t - 1)
    }

    pub fn multiply_pt_poly(ct: &BfvCipherText, pt_poly: &Poly) -> BfvCipherText {
        // plaintext polynomial must be Rt
        assert!(pt_poly.ctx.moduli.q <= ct.params.t);
        let mut pt_poly = pt_poly.clone();
        pt_poly.switch_context(&ct.params.poly_ctx);

        let c0 = &ct.c[0] * &pt_poly;
        let c1 = &ct.c[1] * &pt_poly;
        BfvCipherText::new(&ct.params, vec![c0, c1])
    }

    pub fn multiply_constant(ct: &BfvCipherText, value: u64) -> BfvCipherText {
        // constant must be in Z_t
        assert!(value < ct.params.t);

        let c0 = &ct.c[0] * value;
        let c1 = &ct.c[1] * value;
        BfvCipherText::new(&ct.params, vec![c0, c1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let params = Arc::new(BfvParameters::default());
        let sk = BfvSecretKey::new(&params);
        let pk = BfvPublicKey::new(&sk);
        let pt = BfvPlaintext::new(&params, vec![3, 3, 3, 2]);
        let ct = pk.encrypt(&pt.poly.coeffs.into());
        let pt2 = sk.decrypt(&ct);
        assert_eq!(pt.values.into_vec(), pt2);
    }
}
