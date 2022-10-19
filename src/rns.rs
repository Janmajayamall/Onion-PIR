use crate::poly::Modulus;
use crate::utils::{div_ceil, ilog2, U256};
use crypto_bigint::U192;
use itertools::{izip, Itertools};
use ndarray::{ArrayView1, ArrayViewMut1};
use num_bigint::BigUint;
use num_bigint_dig::{BigUint as BigUintDig, ModInverse};
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use std::{cmp::min, sync::Arc};

#[derive(Debug, PartialEq, Clone)]
pub struct RnsContext {
    pub moduli_u64: Vec<u64>,
    pub moduli: Vec<Modulus>,
    pub q_star: Vec<BigUint>,
    pub q_tilde: Vec<u64>,
    /// garner = q_star * q_tilde
    pub garner: Vec<BigUint>,
    pub product: BigUint,
    pub product_dig: BigUintDig,
}

impl RnsContext {
    pub fn new(moduli_u64: Vec<u64>) -> Self {
        //TODO: Check that moduli are coprime
        let moduli = moduli_u64.iter().map(|m| Modulus::new(*m)).collect();

        let product = moduli_u64.iter().fold(BigUint::one(), |acc, q| acc * *q);
        let product_dig = moduli_u64.iter().fold(BigUintDig::one(), |acc, q| acc * *q);

        let q_star: Vec<BigUint> = moduli_u64.iter().map(|qi| &product / qi).collect();
        let q_tilde: Vec<u64> = moduli_u64
            .iter()
            .map(|qi| {
                (&product_dig / qi)
                    .mod_inverse(BigUintDig::from_u64(*qi).unwrap())
                    .unwrap()
                    .to_u64()
                    .unwrap()
            })
            .collect();

        let garner = izip!(q_star.iter(), q_tilde.iter())
            .map(|(q_s, q_t)| q_s * *q_t)
            .collect();

        Self {
            moduli_u64,
            moduli,
            q_star,
            q_tilde,
            garner,
            product,
            product_dig,
        }
    }

    pub fn project(&self, a: &BigUint) -> Vec<u64> {
        self.moduli_u64
            .iter()
            .map(|q| (a % q).to_u64().unwrap())
            .collect()
    }

    pub fn lift(&self, rests: ArrayView1<u64>) -> BigUint {
        let mut res = BigUint::zero();
        izip!(rests.iter(), self.garner.iter()).for_each(|(x, garner_i)| {
            res += garner_i * *x;
        });
        res
    }

    pub fn modulus(&self) -> &BigUint {
        &self.product
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ScalingFactor {
    numerator: BigUint,
    denominator: BigUint,
    is_one: bool,
}
impl ScalingFactor {
    pub fn new(numerator: BigUint, denominator: BigUint) -> Self {
        assert!(denominator != BigUint::zero());
        ScalingFactor {
            is_one: numerator == denominator,
            numerator,
            denominator,
        }
    }
    pub fn is_one(&self) -> bool {
        self.is_one
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RnsScaler {
    from: Arc<RnsContext>,
    to: Arc<RnsContext>,
    scaling_factor: ScalingFactor,

    gamma: Box<[u64]>,
    gamma_shoup: Box<[u64]>,
    theta_gamma_lo: u64,
    theta_gamma_hi: u64,
    theta_gamma_sign: bool,

    omega: Box<[Box<[u64]>]>,
    omega_shoup: Box<[Box<[u64]>]>,
    theta_omega_lo: Box<[u64]>,
    theta_omega_hi: Box<[u64]>,
    theta_omega_sign: Box<[bool]>,

    theta_garner_lo: Box<[u64]>,
    theta_garner_hi: Box<[u64]>,
    theta_garner_shift: usize,
}

impl RnsScaler {
    /// Create a RNS scaler by numerator / denominator.
    ///
    /// Aborts if denominator is equal to 0.
    pub fn new(
        from: &Arc<RnsContext>,
        to: &Arc<RnsContext>,
        scaling_factor: ScalingFactor,
    ) -> Self {
        // Let's define gamma = round(numerator * from.product / denominator)
        let (gamma, theta_gamma_lo, theta_gamma_hi, theta_gamma_sign) =
            Self::extract_projection_and_theta(
                to,
                &from.product,
                &scaling_factor.numerator,
                &scaling_factor.denominator,
                false,
            );
        let gamma_shoup = izip!(&gamma, &to.moduli)
            .map(|(wi, q)| q.shoup(*wi))
            .collect_vec();

        // Let's define omega_i = round(from.garner_i * numerator / denominator)
        let mut omega = Vec::with_capacity(to.moduli.len());
        let mut omega_shoup = Vec::with_capacity(to.moduli.len());
        for _ in &to.moduli {
            omega.push(vec![0u64; from.moduli.len()].into_boxed_slice());
            omega_shoup.push(vec![0u64; from.moduli.len()].into_boxed_slice());
        }
        let mut theta_omega_lo = Vec::with_capacity(from.garner.len());
        let mut theta_omega_hi = Vec::with_capacity(from.garner.len());
        let mut theta_omega_sign = Vec::with_capacity(from.garner.len());
        for i in 0..from.garner.len() {
            let (omega_i, theta_omega_i_lo, theta_omega_i_hi, theta_omega_i_sign) =
                Self::extract_projection_and_theta(
                    to,
                    &from.garner[i],
                    &scaling_factor.numerator,
                    &scaling_factor.denominator,
                    true,
                );
            for j in 0..to.moduli.len() {
                let qj = &to.moduli[j];
                omega[j][i] = qj.reduce(omega_i[j]);
                omega_shoup[j][i] = qj.shoup(omega[j][i]);
            }
            theta_omega_lo.push(theta_omega_i_lo);
            theta_omega_hi.push(theta_omega_i_hi);
            theta_omega_sign.push(theta_omega_i_sign);
        }

        // Determine the shift so that the sum of the scaled theta_garner fit on an U192
        // (shift + 1) + log(q * n) <= 192
        let theta_garner_shift = min(
            from.moduli_u64
                .iter()
                .map(|qi| {
                    192 - 1
                        - ilog2(
                            ((*qi as u128) * (from.moduli_u64.len() as u128)).next_power_of_two(),
                        )
                })
                .into_iter()
                .min()
                .unwrap() as usize,
            127,
        );
        // Finally, define theta_garner_i = from.garner_i / product, also scaled by
        // 2^127.
        let mut theta_garner_lo = Vec::with_capacity(from.garner.len());
        let mut theta_garner_hi = Vec::with_capacity(from.garner.len());
        for garner_i in &from.garner {
            let mut theta: BigUint =
                ((garner_i << theta_garner_shift) + (&from.product >> 1)) / &from.product;
            let theta_hi: BigUint = &theta >> 64;
            theta -= &theta_hi << 64;
            theta_garner_lo.push(theta.to_u64().unwrap());
            theta_garner_hi.push(theta_hi.to_u64().unwrap());
        }

        Self {
            from: from.clone(),
            to: to.clone(),
            scaling_factor,
            gamma: gamma.into_boxed_slice(),
            gamma_shoup: gamma_shoup.into_boxed_slice(),
            theta_gamma_lo,
            theta_gamma_hi,
            theta_gamma_sign,
            omega: omega.into_boxed_slice(),
            omega_shoup: omega_shoup.into_boxed_slice(),
            theta_omega_lo: theta_omega_lo.into_boxed_slice(),
            theta_omega_hi: theta_omega_hi.into_boxed_slice(),
            theta_omega_sign: theta_omega_sign.into_boxed_slice(),
            theta_garner_lo: theta_garner_lo.into_boxed_slice(),
            theta_garner_hi: theta_garner_hi.into_boxed_slice(),
            theta_garner_shift,
        }
    }

    // Let's define gamma = round(numerator * input / denominator)
    // and theta_gamma such that theta_gamma = numerator * input / denominator -
    // gamma. This function projects gamma in the RNS context, and scales
    // theta_gamma by 2**127 and rounds. It outputs the projection of gamma in the
    // RNS context, and theta_lo, theta_hi, theta_sign such that theta_gamma =
    // (-1)**theta_sign * (theta_lo + 2^64 * theta_hi).
    // Ref - https://github.com/Janmajayamall/fhe.rs/blob/6361fa3ce322b16551cfe4856a49e3933d85c872/crates/fhe-math/src/rns/mod.rs#L115
    fn extract_projection_and_theta(
        ctx: &RnsContext,
        input: &BigUint,
        numerator: &BigUint,
        denominator: &BigUint,
        round_up: bool,
    ) -> (Vec<u64>, u64, u64, bool) {
        let gamma = (numerator * input + (denominator >> 1)) / denominator;
        let projected = ctx.project(&gamma);

        let mut theta = (numerator * input) % denominator;
        let mut theta_sign = false;
        if denominator > &BigUint::one() {
            // If denominator is odd, flip theta if theta > (denominator >> 1)
            if denominator & BigUint::one() == BigUint::one() {
                if theta > (denominator >> 1) {
                    theta_sign = true;
                    theta = denominator - theta;
                }
            } else {
                // denominator is even, flip if theta >= (denominator >> 1)
                if theta >= (denominator >> 1) {
                    theta_sign = true;
                    theta = denominator - theta;
                }
            }
        }
        // theta = ((theta << 127) + (denominator >> 1)) / denominator;
        // We can now split theta into two u64 words.
        if round_up {
            if theta_sign {
                theta = (theta << 127) / denominator;
            } else {
                theta = ((theta << 127) + denominator - BigUint::one()) / denominator;
            }
        } else if theta_sign {
            theta = ((theta << 127) + denominator - BigUint::one()) / denominator;
        } else {
            theta = (theta << 127) / denominator;
        }
        let theta_hi_biguint: BigUint = &theta >> 64;
        theta -= &theta_hi_biguint << 64;
        let theta_lo = theta.to_u64().unwrap();
        let theta_hi = theta_hi_biguint.to_u64().unwrap();

        (projected, theta_lo, theta_hi, theta_sign)
    }

    /// Compute the RNS representation of the rests scaled by numerator *
    /// denominator, and either rounded or floored, and store the result in
    /// `out`.
    ///
    /// Aborts if the number of rests is different than the number of moduli in
    /// debug mode, or if the size of out is not in [1, ..., rests.len()].
    pub fn scale(
        &self,
        rests: ArrayView1<u64>,
        mut out: ArrayViewMut1<u64>,
        starting_index: usize,
    ) {
        debug_assert_eq!(rests.len(), self.from.moduli_u64.len());
        debug_assert!(!out.is_empty());
        debug_assert!(starting_index + out.len() <= self.to.moduli_u64.len());

        // First, let's compute the inner product of the rests with theta_omega.
        let mut sum_theta_garner = U192::ZERO;
        for (thetag_lo, thetag_hi, ri) in izip!(
            self.theta_garner_lo.iter(),
            self.theta_garner_hi.iter(),
            rests
        ) {
            let lo = (*ri as u128) * (*thetag_lo as u128);
            let hi = (*ri as u128) * (*thetag_hi as u128) + (lo >> 64);
            sum_theta_garner = sum_theta_garner.wrapping_add(&U192::from_words([
                lo as u64,
                hi as u64,
                (hi >> 64) as u64,
            ]));
        }
        // Let's compute v = round(sum_theta_garner / 2^theta_garner_shift)
        sum_theta_garner >>= self.theta_garner_shift - 1;
        let v = <[u64; 3]>::from(sum_theta_garner);
        let v = div_ceil((v[0] as u128) | ((v[1] as u128) << 64), 2);

        // If the scaling factor is not 1, compute the inner product with the
        // theta_omega
        let mut w_sign = false;
        let mut w = 0u128;
        if !self.scaling_factor.is_one {
            let mut sum_theta_omega = U256::zero();
            for (thetao_lo, thetao_hi, thetao_sign, ri) in izip!(
                self.theta_omega_lo.iter(),
                self.theta_omega_hi.iter(),
                self.theta_omega_sign.iter(),
                rests
            ) {
                let lo = (*ri as u128) * (*thetao_lo as u128);
                let hi = (*ri as u128) * (*thetao_hi as u128) + (lo >> 64);
                if *thetao_sign {
                    sum_theta_omega.wrapping_sub_assign(U256::from([
                        lo as u64,
                        hi as u64,
                        (hi >> 64) as u64,
                        0,
                    ]));
                } else {
                    sum_theta_omega.wrapping_add_assign(U256::from([
                        lo as u64,
                        hi as u64,
                        (hi >> 64) as u64,
                        0,
                    ]));
                }
            }

            // Let's subtract v * theta_gamma to sum_theta_omega.
            let vt_lo_lo = ((v as u64) as u128) * (self.theta_gamma_lo as u128);
            let vt_lo_hi = ((v as u64) as u128) * (self.theta_gamma_hi as u128);
            let vt_hi_lo = ((v >> 64) as u128) * (self.theta_gamma_lo as u128);
            let vt_hi_hi = ((v >> 64) as u128) * (self.theta_gamma_hi as u128);
            let vt_mi =
                (vt_lo_lo >> 64) + ((vt_lo_hi as u64) as u128) + ((vt_hi_lo as u64) as u128);
            let vt_hi = (vt_lo_hi >> 64) + (vt_mi >> 64) + ((vt_hi_hi as u64) as u128);
            if self.theta_gamma_sign {
                sum_theta_omega.wrapping_add_assign(U256::from([
                    vt_lo_lo as u64,
                    vt_mi as u64,
                    vt_hi as u64,
                    0,
                ]))
            } else {
                sum_theta_omega.wrapping_sub_assign(U256::from([
                    vt_lo_lo as u64,
                    vt_mi as u64,
                    vt_hi as u64,
                    0,
                ]))
            }

            // Let's compute w = round(sum_theta_omega / 2^127).
            w_sign = sum_theta_omega.msb() > 0;

            if w_sign {
                w = u128::from(&((!sum_theta_omega) >> 126)) + 1;
                w /= 2;
            } else {
                w = u128::from(&(sum_theta_omega >> 126));
                w = div_ceil(w, 2)
            }
        }

        unsafe {
            for i in 0..out.len() {
                debug_assert!(starting_index + i < self.to.moduli.len());
                debug_assert!(starting_index + i < self.omega.len());
                debug_assert!(starting_index + i < self.omega_shoup.len());
                debug_assert!(starting_index + i < self.gamma.len());
                debug_assert!(starting_index + i < self.gamma_shoup.len());
                let out_i = out.get_mut(i).unwrap();
                let qi = self.to.moduli.get_unchecked(starting_index + i);
                let omega_i = self.omega.get_unchecked(starting_index + i);
                let omega_shoup_i = self.omega_shoup.get_unchecked(starting_index + i);
                let gamma_i = self.gamma.get_unchecked(starting_index + i);
                let gamma_shoup_i = self.gamma_shoup.get_unchecked(starting_index + i);

                let mut yi = (qi.modulus() * 2
                    - qi.lazy_mul_shoup(qi.reduce_u128(v), *gamma_i, *gamma_shoup_i))
                    as u128;

                if !self.scaling_factor.is_one {
                    let wi = qi.lazy_reduce_u128(w);
                    yi += if w_sign { qi.modulus() * 2 - wi } else { wi } as u128;
                }

                debug_assert!(rests.len() <= omega_i.len());
                debug_assert!(rests.len() <= omega_shoup_i.len());
                for j in 0..rests.len() {
                    yi += qi.lazy_mul_shoup(
                        *rests.get(j).unwrap(),
                        *omega_i.get_unchecked(j),
                        *omega_shoup_i.get_unchecked(j),
                    ) as u128;
                }

                *out_i = qi.reduce_u128(yi)
            }
        }
    }
}
