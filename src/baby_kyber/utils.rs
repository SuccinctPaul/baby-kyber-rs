use crate::poly::Polynomial;
use crate::ring::{PolynomialRingTrait, Ring};
use rand::{Rng, RngCore};
use std::ops::Neg;

// a randomizer small polynomial vector
pub fn random_small_poly_vector<P: PolynomialRingTrait>(
    rng: &mut impl RngCore,
    dimension: usize,
) -> Vec<P> {
    let mut s = vec![];
    let degree = 6;
    // let value = rng.next_u64() % 6; ;
    let value = 6;
    for _ in 0..dimension {
        // TODO: optimize here to use random one.
        let poly = {
            P::from_coefficients(
                (0..degree)
                    .map(|i| {
                        let value = value + i as u64;
                        P::PolyCoeff::from(value)
                    })
                    .collect(),
            )
            //     P::rand_with_bound_degree(rng)
        };
        println!("randomer poly: {:?}", poly.to_string());
        s.push(poly);
    }
    s
}

// a randomizer polynomial vector with `bound_degree`
pub fn random_poly_vector<P: PolynomialRingTrait>(
    rng: &mut impl RngCore,
    dimension: usize,
) -> Vec<P> {
    (0..dimension)
        .map(|_| P::rand_with_bound_degree(rng))
        .collect()
}
