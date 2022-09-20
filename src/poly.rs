use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

#[derive(Clone)]
pub struct Modulus {
    q: u64,
}

impl Modulus {
    fn add_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a + b) % self.q
    }

    fn sub_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a + self.q - b) % self.q
    }

    fn mul_mod(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a <= self.q);
        debug_assert!(b <= self.q);

        (a * b) % self.q
    }
}

#[derive(Clone)]
pub struct Context {
    moduli: Modulus,
    degree: usize,
}

pub struct Poly {
    pub coeffs: Vec<u64>,
    pub ctx: Context,
}

// define operations
impl Poly {
    pub fn zero(ctx: Context) -> Self {
        Poly {
            coeffs: vec![0; ctx.degree],
            ctx,
        }
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Self) {
        debug_assert!(self.coeffs.len() == rhs.coeffs.len());

        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(a, b)| *a = self.ctx.moduli.add_mod(*a, *b));
    }
}

impl Add for Poly {
    type Output = Poly;

    fn add(self, mut rhs: Poly) -> Self::Output {
        rhs += &self;
        rhs
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        debug_assert!(self.coeffs.len() == rhs.coeffs.len());

        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(a, b)| *a = self.ctx.moduli.sub_mod(*a, *b));
    }
}

impl Sub for Poly {
    type Output = Poly;
    fn sub(self, mut rhs: Poly) -> Self::Output {
        rhs -= &self;
        rhs
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, rhs: &Poly) {
        debug_assert!(self.coeffs.len() == rhs.coeffs.len());

        let mut ctx_res = self.ctx.clone();
        ctx_res.degree += rhs.ctx.degree;
        let mut res = Poly::zero(ctx_res);
        self.coeffs.iter().enumerate().for_each(|(i, a)| {
            rhs.coeffs.iter().enumerate().for_each(|(j, b)| {
                res.coeffs[i + j] += res.ctx.moduli.mul_mod(*a, *b);
            })
        });

        *self = res;
    }
}

impl Mul for Poly {
    type Output = Poly;
    fn mul(self, mut rhs: Poly) -> Self::Output {
        rhs *= &self;
        rhs
    }
}
