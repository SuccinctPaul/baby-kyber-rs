use crate::poly::Polynomial;
use crate::ring::PolynomialRingTrait;
use rand::{Rng, RngCore};

// a randomizer polynomial vector
// These polynomial vectors are freshly generated for every encryption.
pub fn small_poly_vector<P: PolynomialRingTrait>(
    rng: &mut impl RngCore,
    dimension: usize,
    degree: usize,
) -> Vec<P> {
    let mut s = vec![];
    // let value = rng.next_u64(); ;
    let value = 6;
    for _ in 0..dimension {
        // TODO: optimize here to use random one.
        s.push(P::from(P::PolyType::from_coefficients(vec![
            P::PolyCoeff::from(value);
            degree
        ])));
    }
    s
}

// a randomizer polynomial vector
// These polynomial vectors are freshly generated for every encryption.
pub fn random_poly_vector<P: Polynomial>(
    rng: &mut impl RngCore,
    dimension: usize,
    degree: usize,
) -> Vec<P> {
    let mut s = vec![];
    for _ in 0..dimension {
        s.push(P::rand(rng, degree));
    }
    s
}
