use crate::ring::Ring;

pub struct VectorArithmatic<R: Ring> {
    phantom: std::marker::PhantomData<R>,
}

impl<R: Ring> VectorArithmatic<R> {
    // Dot product between two vectors.
    pub fn vec_dot_mul(a: &Vec<R>, b: &Vec<R>) -> R {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() * bi.clone())
            .fold(R::zero(), |acc, x| acc + x)
    }

    // Mul between two vectors.
    // It's a special case of matrix mul.
    pub fn vec_mul(a: &Vec<R>, b: &Vec<R>) -> Vec<R> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() * bi.clone())
            .collect::<Vec<_>>()
    }

    // Add between two vectors.
    pub fn vec_add(a: &Vec<R>, b: &Vec<R>) -> Vec<R> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() + bi.clone())
            .collect::<Vec<_>>()
    }

    // Add between two vectors.
    pub fn vec_sub(a: &Vec<R>, b: &Vec<R>) -> Vec<R> {
        assert_eq!(a.len(), b.len(), "Vectors must have the same length");

        a.iter()
            .zip(b.iter())
            .map(|(ai, bi)| ai.clone() - bi.clone())
            .collect::<Vec<_>>()
    }
}
