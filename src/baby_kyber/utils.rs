use crate::poly::Polynomial;
use crate::ring::PolynomialRingTrait;
use rand::{Rng, RngCore};

// a randomizer polynomial vector
// These polynomial vectors are freshly generated for every encryption.
pub fn small_poly_vector<P: PolynomialRingTrait>(
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

// a randomizer polynomial vector
// These polynomial vectors are freshly generated for every encryption.
pub fn random_poly_vector<P: PolynomialRingTrait>(
    rng: &mut impl RngCore,
    dimension: usize,
) -> Vec<P> {
    let mut s = vec![];
    let degree = 6;
    // let value = rng.next_u64() % 6; ;
    let value = 5;
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
