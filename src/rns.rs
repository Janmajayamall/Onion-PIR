use crate::poly::Modulus;

pub struct RnsContext {
    pub moduli_64: Vec<u64>,
    pub moduli: Vec<Modulus>,
}

impl RnsContext {
    pub fn new(moduli_64: Vec<u64>) -> Self {
        //TODO: Check that moduli are coprime
        let moduli = moduli_64.iter().map(|m| Modulus::new(*m)).collect();
        Self { moduli_64, moduli }
    }
}
