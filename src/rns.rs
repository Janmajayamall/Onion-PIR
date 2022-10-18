use crate::poly::Modulus;
use crate::utils::ilog2;
use itertools::{izip, Itertools};
use num_bigint::BigUint;
use num_traits::{AsPrimitive, One, ToPrimitive, Zero};
use std::{cmp::min, sync::Arc};

#[derive(Debug, PartialEq, Clone)]
pub struct RnsContext {
    pub moduli_64: Vec<u64>,
    pub moduli: Vec<Modulus>,
    pub product: BigUint,
}

impl RnsContext {
    pub fn new(moduli_64: Vec<u64>) -> Self {
        //TODO: Check that moduli are coprime
        let moduli = moduli_64.iter().map(|m| Modulus::new(*m)).collect();

        let product = moduli_64.iter().fold(BigUint::one(), |acc, q| acc * *q);

        Self {
            moduli_64,
            moduli,
            product,
        }
    }

    pub fn project(&self, a: &BigUint) -> Vec<u64> {
        self.moduli_64
            .iter()
            .map(|q| (a % q).to_u64().unwrap())
            .collect()
    }
}

pub struct ScalingFactor {
    numerator: BigUint,
    denominator: BigUint,
}
impl ScalingFactor {
    pub fn new(numerator: BigUint, denominator: BigUint) -> Self {
        assert!(denominator != BigUint::zero());
        ScalingFactor {
            numerator,
            denominator,
        }
    }
}

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
}
