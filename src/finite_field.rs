use num_bigint::BigUint;
use thiserror::Error;
use num_traits::Zero;

pub struct FiniteField {
    pub p: BigUint,
}

#[derive(Error, Debug, PartialEq)]
pub enum FiniteFieldError {
    #[error("Cannot divide by zero in a finite field")]
    DivisionByZero,
    #[error("Zero does not have modular inverse")]
    ZeroHasNoInverse,
    #[error("Modular inverse does not exist (modulus might not be prime)")]
    NoInverse,
}

impl FiniteField {
    pub fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a + b) % &self.p
    }

    pub fn sub(&self, a: &BigUint, b: &BigUint) -> BigUint {
        // Optimization: a -b = ((a + p - (b mod p)) mod p
        (a + &self.p - (b % &self.p)) % &self.p
    }

    pub fn add_inv(&self, a: &BigUint) -> BigUint {
        // a + x = 0 mod p --> x = p - a
        // There are cases that needs to be aware of:
        // 1. a = 0, p - 0 = p, but should be 0.
        // 2. if a > p, substraction will crash.

        // 1. (a % p) ensures the number is within the field range, handles if a > p.
        // 2. (p - (a % p)) finds the difference
        // 3. The final % p ensures that if the result is p, it becomes 0
        (&self.p - (a % &self.p)) % &self.p
    }

    pub fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        (a * b) % &self.p
    }

    pub fn div(&self, a: &BigUint, b: &BigUint) -> Result<BigUint, FiniteFieldError> {
        if b.is_zero() {
            return Err(FiniteFieldError::DivisionByZero);
        }
        Ok((a * self.mul_inv(b)?) % &self.p)
    }

    pub fn pow(&self, n: &BigUint, a: &BigUint) -> BigUint {
        a.modpow(n, &self.p)
    }

    pub fn mul_inv(&self, a: &BigUint) -> Result<BigUint, FiniteFieldError> {
        if a.is_zero() {
            return Err(FiniteFieldError::ZeroHasNoInverse);
        }
        a.modinv(&self.p).ok_or(FiniteFieldError::NoInverse)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::{BigUint};

    fn bu(n: u64) -> BigUint {
        BigUint::from(n)
    }

    fn setup_field() -> FiniteField {
        FiniteField { p: bu(17) }
    }

    #[test]
    fn test_addition() {
        let field = setup_field();
        assert_eq!(field.add(&bu(10), &bu(10)), bu(3));
        assert_eq!(field.add(&bu(0), &bu(5)), bu(5));
        assert_eq!(field.add(&bu(12), &bu(5)), bu(0));
    }

    #[test]
    fn test_subtraction() {
        let field = setup_field();
        // Standard: 10 - 5 = 5
        assert_eq!(field.sub(&bu(10), &bu(5)), bu(5));
        // Underflow: 5 - 10 = -5. -5 mod 17 = 12
        assert_eq!(field.sub(&bu(5), &bu(10)), bu(12));
        // Large input: 20 - 5. (20 is 3 mod 17). 3 - 5 = -2. -2 mod 17 = 15
        assert_eq!(field.sub(&bu(20), &bu(5)), bu(15));
    }

    #[test]
    fn test_additive_inverse() {
        let field = setup_field();
        // Inv(5) = 12
        assert_eq!(field.add_inv(&bu(5)), bu(12));
        // Inv(0) = 0
        assert_eq!(field.add_inv(&bu(0)), bu(0));
        // Inv(18) -> Inv(1) = 16
        assert_eq!(field.add_inv(&bu(18)), bu(16));
    }

    #[test]
    fn test_multiplication() {
        let field = setup_field();
        assert_eq!(field.mul(&bu(3), &bu(10)), bu(13));
        assert_eq!(field.mul(&bu(0), &bu(5)), bu(0));
        assert_eq!(field.mul(&bu(11), &bu(14)), bu(1));
    }

    #[test]
    fn test_multiplicative_inverse() {
        let field = setup_field();
        // 3 * x = 1 mod 17 -> x = 6 (3*6=18, 18 mod 17 = 1)
        assert_eq!(field.mul_inv(&bu(3)).unwrap(), bu(6));

        // Test Error: Zero
        assert_eq!(
            field.mul_inv(&bu(0)).unwrap_err(),
            FiniteFieldError::ZeroHasNoInverse
        );

        // Test Error: No Inverse (using composite modulus)
        let composite_field = FiniteField { p: bu(10) };
        // gcd(2, 10) = 2, so 2 has no inverse mod 10
        assert_eq!(
            composite_field.mul_inv(&bu(2)).unwrap_err(),
            FiniteFieldError::NoInverse
        );
    }

    #[test]
    fn test_division() {
        let field = setup_field();
        // 10 / 2 = 5
        assert_eq!(field.div(&bu(10), &bu(2)).unwrap(), bu(5));

        // 1 / 3 mod 17 -> 1 * 6 = 6
        assert_eq!(field.div(&bu(1), &bu(3)).unwrap(), bu(6));

        // Test Error: Divide by zero
        assert_eq!(
            field.div(&bu(10), &bu(0)).unwrap_err(),
            FiniteFieldError::DivisionByZero
        );
    }

    #[test]
    fn test_power() {
        let field = setup_field();
        // 2^4 = 16 mod 17 = 16
        assert_eq!(field.pow(&bu(4), &bu(2)), bu(16));
        // 3^0 = 1
        assert_eq!(field.pow(&bu(0), &bu(3)), bu(1));
        // 0^5 = 0
        assert_eq!(field.pow(&bu(5), &bu(0)), bu(0));
    }
}