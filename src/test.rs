use crate::baby_kyber::BabyKyber;
use crate::param::{BABY_KYBER_BOUND_DEGREE, BABY_KYBER_DIMENSION};
use crate::poly::uni_poly::UniPolynomial;
use crate::ring::ring_poly::RingPolynomial;
use crate::ring::zq::Zq;
#[test]
fn test_baby_kyber() {
    let msg = 11;

    type BABY_KYBER_POLY_RING = RingPolynomial<UniPolynomial<Zq>, BABY_KYBER_BOUND_DEGREE>;

    let rng = &mut rand::thread_rng();

    let baby_kyber = BabyKyber::<BABY_KYBER_POLY_RING, BABY_KYBER_DIMENSION>::init();

    let (private_key, public_key) = baby_kyber.keygen(rng);

    println!("\n\n==encrypting");
    let ciphertext = baby_kyber.encrypto(rng, msg, &public_key);
    println!("\n\n==start decrypto");
    let is_verified = BabyKyber::<BABY_KYBER_POLY_RING, BABY_KYBER_DIMENSION>::decrypto(
        msg,
        ciphertext,
        &private_key,
    );
    assert!(is_verified, "Failed decryption");
    println!("==successfully verified decrypted message");
}
