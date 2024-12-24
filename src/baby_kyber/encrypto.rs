use crate::baby_kyber::constant::ERROR_POLY_DEGREE;
use crate::baby_kyber::keygen::PublickKey;
use crate::baby_kyber::{BabyKyber, Ciphtertexts};
use crate::matrix::Matrix;
use crate::matrix::poly_matrix::PolyMatrix;
use crate::poly::{Polynomial, random_poly_vector};
use crate::ring::Ring;
use crate::utils::{bits_normalize, bytes_to_le_bits};
use rand::RngCore;

impl<P: Polynomial> BabyKyber<P> {
    // Kyber ciphtertexts consist of those two values: (u,v)
    pub fn encrypto(
        &self,
        rng: &mut impl RngCore,
        msg: &[u8],
        public_key: &PublickKey<P>,
    ) -> Ciphtertexts<P> {
        // 1. generate msg poly
        let poly_m = {
            // 1.1 encode msg into bits
            let msg_bits = bytes_to_le_bits(&msg);
            let msg_bits = bits_normalize(&msg_bits);
            println!("msg_bits {:?}", msg_bits);
            // 1.2. generate msg poly according the msg_bits
            let msg_poly_coeffs = msg_bits
                .into_iter()
                .map(|bit| {
                    if bit {
                        P::Coefficient::one()
                    } else {
                        P::Coefficient::zero()
                    }
                })
                .collect::<Vec<_>>();
            let msg_poly = P::from_coefficients(msg_poly_coeffs);
            // 1.3 scale polynomial by the (q/2)
            let scalar = Self::compute_scalar();
            let scalared_coeffs = msg_poly
                .coefficients()
                .into_iter()
                .map(|coeff| coeff * scalar)
                .collect::<Vec<_>>();
            P::from_coefficients(scalared_coeffs)
        };

        // 2. generate random poly and error poly
        //      These polynomial vectors are freshly generated for every encryption.
        let (poly_e1, poly_e2, poly_r) = {
            let poly_e1 = random_poly_vector(rng, self.dimension, self.degree);
            let poly_e2 = P::rand(rng, self.degree);
            let poly_r = random_poly_vector(rng, self.dimension, self.degree);

            (poly_e1, poly_e2, poly_r)
        };

        // 3. compute u,v
        let (u, v) = {
            let A_transpose = public_key.A.transpose();
            let t = public_key.t.clone();

            // 3.1 comptue u=A^T * r + e1
            let u = {
                let A_transpose_dot_r = A_transpose.mul_vector(&poly_r);
                PolyMatrix::vec_add(&A_transpose_dot_r, &poly_e1)
            };
            assert_eq!(u.len(), self.dimension, "u.len != dimension");

            // 3.2 compute  v=t^T * r + e2 + m
            let v = {
                let t_dot_r = PolyMatrix::vec_mul(&t, &poly_r);
                // assert_eq!(t_dot_r.len(), 1, "t_dot_r.len() != 1");
                // TODO: t_dot_r needs to dig more.
                let t_mul_r_plus_e2 = t_dot_r[0].clone() + poly_e2;

                t_mul_r_plus_e2 + poly_m
            };
            (u, v)
        };

        Ciphtertexts { u, v }
    }

    // // a randomizer polynomial vector 'r'
    // // r = [-x^3 + x^2, x^3 + x^2 - 1]
    // pub(crate) fn random_r(rng: &mut impl RngCore) -> Vec<R> {
    //     let x = R::rand(rng);
    //
    //     //  f(x)= -x^3 + x^2
    //     let poly_f = RingPolynomial::from_coefficients(vec![R::zero(), R::one(), R::one().neg()]);
    //
    //     // g(x) = x^3 + x^2 - 1
    //     let poly_g = RingPolynomial::from_coefficients(vec![R::one().neg(), R::zero(), R::one()]);
    //
    //     vec![poly_f.evaluate(&x), poly_g.evaluate(&x)]
    // }
    //
    // // an error polynomial `e2`
    // // e = [-x^3 - x^2]
    // pub(crate) fn random_e2(rng: &mut impl RngCore) -> Vec<R> {
    //     let x = R::rand(rng);
    //
    //     //  f(x)= -x^3 - x^2
    //     let poly =
    //         RingPolynomial::from_coefficients(vec![R::zero(), R::one().neg(), R::one().neg()]);
    //
    //     vec![poly.evaluate(&x)]
    // }

    // To encrypt a message, using the message’s binary representation and it into a polynomial.
    // Every bit of the message is used as a coefficient
}
