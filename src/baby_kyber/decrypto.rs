use crate::baby_kyber::keygen::PrivateKey;
use crate::baby_kyber::{BabyKyber, Ciphtertexts};
use crate::matrix::poly_matrix::PolyMatrix;
use crate::poly::Polynomial;
use crate::ring::Ring;
use crate::utils::{bits_normalize, bytes_to_le_bits};
use rand::RngCore;
use std::ops::{Div, Sub};

impl<P: Polynomial> BabyKyber<P> {
    // Kyber ciphtertexts consist of those two values: (u,v)
    pub fn decrypto(
        msg: &[u8],
        ciphter_text: Ciphtertexts<P>,
        private_key: &PrivateKey<P>,
    ) -> bool {
        // 1. First, we compute a noisy result m_n
        //      m_n = v - s^T * u
        let m_n = {
            let s_tranpose_dot_u = PolyMatrix::vec_mul(&private_key.s, &ciphter_text.u);
            ciphter_text.v - s_tranpose_dot_u[0].clone()
        };
        // 2. round and scalar m_n
        let scalaed_m_n = {
            // 2.1 round m_n and
            let rounded_coeffs = Self::round(m_n.coefficients());
            // 2.2 mul 1/scalar to recover the msg_bin
            let scalar = Self::compute_scalar();
            let scalared_coeffs = rounded_coeffs
                .into_iter()
                .map(|m| m.div(scalar.clone()))
                .collect::<Vec<_>>();

            P::from_coefficients(scalared_coeffs)
        };

        // 4. check msg_bin,
        let expect_msg_bin = {
            // 4.1 compute msg_bin
            // 1.1 encode msg into bits
            let msg_bits = bytes_to_le_bits(&msg);
            let msg_bits = bits_normalize(&msg_bits);
            // 1.2. generate msg poly according the msg_bits
            let coeffs = msg_bits
                .into_iter()
                .map(|bit| {
                    if bit {
                        P::Coefficient::one()
                    } else {
                        P::Coefficient::zero()
                    }
                })
                .collect::<Vec<_>>();
            P::from_coefficients(coeffs)
        };

        // 5. final check
        scalaed_m_n == expect_msg_bin
    }

    // Round the value to (module/2, 0)
    pub fn round(coeffs: Vec<P::Coefficient>) -> Vec<P::Coefficient> {
        let scalar = Self::compute_scalar();

        coeffs
            .into_iter()
            .map(|c| {
                if scalar.sub(c.clone()) > c.sub(P::Coefficient::zero()) {
                    scalar.clone()
                } else {
                    P::Coefficient::zero()
                }
            })
            .collect::<Vec<_>>()
    }
}
