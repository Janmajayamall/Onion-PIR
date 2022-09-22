use std::sync::Arc;

use crate::poly::{Context, Modulus};

use super::poly::Poly;

#[derive(Debug)]
pub struct BfvParameters {
    t: u64,
    q: u64,
    n: usize,
    st_dev: f64,
    poly_ctx: Arc<Context>,
}

impl Default for BfvParameters {
    fn default() -> Self {
        let q: u64 = 65536;
        let degree: usize = 4;
        let ctx = Context::new(Modulus { q }, degree);
        Self {
            t: 4,
            q,
            n: degree,
            st_dev: 3.2,
            poly_ctx: Arc::new(ctx),
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

    pub fn decrypt(&self, ct: &BfvCipherText) -> Vec<u64> {
        let mut poly = &ct.c[1] * &self.poly + ct.c[0].clone();
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

    pub fn encrypt(&self, plaintext: &BfvPlaintext) -> BfvCipherText {
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
            .zip(plaintext.poly.coeffs.iter())
            .for_each(|(sv, v)| *sv = v * delta);

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
        let mut poly = Poly::zero(&params.poly_ctx);
        poly.coeffs
            .iter_mut()
            .zip(pt.iter().map(|v| *v % params.t))
            .for_each(|(c, v)| *c = params.poly_ctx.moduli.convert(v));

        BfvPlaintext {
            params: params.clone(),
            values: pt.into_boxed_slice(),
            poly,
        }
    }
}

#[derive(Debug)]
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encrypt_decrypt() {
        let params = Arc::new(BfvParameters::default());
        let sk = BfvSecretKey::new(&params);
        let pk = BfvPublicKey::new(&sk);
        let pt = BfvPlaintext::new(&params, vec![1, 2, 3, 2]);
        let ct = pk.encrypt(&pt);
        let pt2 = sk.decrypt(&ct);
        assert_eq!(pt.values.into_vec(), pt2);
    }
}
