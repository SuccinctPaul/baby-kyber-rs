pub mod zq;

use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use std::random::Random;

pub trait Ring:
    Add<Output = Self>
    + AddAssign
    + Mul<Output = Self>
    + MulAssign
    + Neg<Output = Self>
    + Sub<Output = Self>
    + SubAssign
    + Div<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + Sized
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> MulAssign<&'a Self>
    + Sized
    + Clone
    + Copy
    + Debug
    + PartialEq
    + Eq
    + Ord
    + PartialOrd
    + Random
    + From<u64>
    + Send
    + Sync
    // + ToString
    + Display
{
    const MODULUS: u64;

    fn rand(rng: &mut impl rand::RngCore) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn square(&self) -> Self;
    /// Computes self^exponent using exponentiation by squaring
    fn pow(&self, power: u64) -> Self;
}
