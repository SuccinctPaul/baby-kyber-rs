use crate::baby_kyber::BabyKyber;
use crate::poly::uni_poly::UniPolynomial;
use crate::ring::zq::Zq;

#[test]
fn test_baby_kyber() {
    let msg = 11;

    let rng = &mut rand::thread_rng();

    let baby_kyber = BabyKyber::<UniPolynomial<Zq>>::init(2, 2);

    let (private_key, public_key) = baby_kyber.keygen(rng);

    let ciphertext = baby_kyber.encrypto(rng, msg, &public_key);
    let is_verified = BabyKyber::<UniPolynomial<Zq>>::decrypto(msg, ciphertext, &private_key);
    assert!(is_verified, "Failed decryption");
    println!("==successfully verified decrypted message");
}
