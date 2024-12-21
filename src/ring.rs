pub mod fq;

use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub};
use std::random::Random;

pub trait Ring:
    Add<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Sized
    + Clone
    + Copy
    + Debug
    + PartialEq
    + Eq
    + Ord
    + PartialOrd
    + Random
{
    const MODULUS: u64;

    fn rand(rng: &mut impl rand::RngCore) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
}
