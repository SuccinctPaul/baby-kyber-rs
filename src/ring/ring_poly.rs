use crate::poly::Polynomial;
use crate::ring::{PolynomialRingTrait, Ring};
use rand::RngCore;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

// Z_q[x]/(x^n+1), N mean the module poly's degree
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct RingPolynomial<P: Polynomial, const N: u64> {
    pub inner: P,
}
impl<P: Polynomial, const N: u64> RingPolynomial<P, N> {
    pub fn new(inner: P) -> Self {
        Self { inner }
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> PolynomialRingTrait
    for RingPolynomial<P, DEGREE_BOUND>
{
    type PolyType = P;

    // x^n + 1
    fn modulus(&self) -> Self::PolyType {
        let coeffs = if DEGREE_BOUND == 0 {
            vec![P::Coefficient::from(2)]
        } else {
            let mut coeffs = vec![P::Coefficient::zero(); DEGREE_BOUND as usize + 1];
            coeffs[0] = P::Coefficient::one();
            coeffs[DEGREE_BOUND as usize] = P::Coefficient::one();
            coeffs
        };
        P::from_coefficients(coeffs)
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> Add for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.inner + rhs.inner)
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> AddAssign for RingPolynomial<P, DEGREE_BOUND> {
    fn add_assign(&mut self, rhs: Self) {
        self.inner += rhs.inner
    }
}

impl<'a, P: Polynomial, const DEGREE_BOUND: u64> Add<&'a Self> for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.inner + &rhs.inner)
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> Mul for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let product = self.inner.clone() * rhs.inner;
        let module_one = product % self.modulus();
        Self::new(module_one)
    }
}

impl<'a, P: Polynomial, const DEGREE_BOUND: u64> Mul<&'a Self> for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        let product = self.inner.clone() * &rhs.inner;
        let module_one = product % self.modulus();
        Self::new(module_one)
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> MulAssign for RingPolynomial<P, DEGREE_BOUND> {
    fn mul_assign(&mut self, rhs: Self) {
        let product = self.inner.clone() * &rhs.inner;
        let module_one = product % self.modulus();
        self.inner = module_one;
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> Sub for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.inner - rhs.inner)
    }
}

impl<'a, P: Polynomial, const DEGREE_BOUND: u64> Sub<&'a Self> for RingPolynomial<P, DEGREE_BOUND> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.inner + &rhs.inner)
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> SubAssign for RingPolynomial<P, DEGREE_BOUND> {
    fn sub_assign(&mut self, rhs: Self) {
        self.inner -= rhs.inner
    }
}

impl<P: Polynomial, const DEGREE_BOUND: u64> Display for RingPolynomial<P, DEGREE_BOUND> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        self.inner.to_string();
        Ok(())
    }
}
