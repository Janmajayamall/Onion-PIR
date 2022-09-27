use super::bfv::BfvCipherText;

/// Structure:
/// Params: \beta & l
/// [
///     [
///         RLWE( 1/beta^1 * -s * m)
///         RLWE( 1/beta^2 * -s * m)
///         .
///         .
///         .
///         RLWE( 1/beta^l * -s * m)
///     ]
///     [
///         RLWE( 1/beta^1 * m)
///         RLWE( 1/beta^2 * m)
///         .
///         .
///         .
///         RLWE( 1/beta^l * m)
///     ]
/// ]
/// E (2l X 2)
struct Rgsw {
    cts: Vec<Vec<BfvCipherText>>,
    beta: u64,
    l: usize,
}

impl Rgsw {
    fn external_product(rgsw: Rgsw, rlwe: BfvCipherText) {
        // Steps
        // (1) Decompose B and A of rlwe into B` and A`
        // (2) Perform <B`, RGWS[0]> and <A`, RGSW[1]>
    }
}
