use num_bigint::BigUint;
use thiserror::Error;
use num_traits::Zero;

#[derive(Debug)]
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

    fn setup_field() -> FiniteField {
        FiniteField { p: 17u32.into() }
    }

    #[test]
    fn test_addition() {
        let field = setup_field();
        assert_eq!(field.add(&10u32.into(), &10u32.into()), 3u32.into());
        assert_eq!(field.add(&0u32.into(), &5u32.into()), 5u32.into());
        assert_eq!(field.add(&12u32.into(), &5u32.into()), 0u32.into());
    }

    #[test]
    fn test_subtraction() {
        let field = setup_field();
        // Standard: 10 - 5 = 5
        assert_eq!(field.sub(&10u32.into(), &5u32.into()), 5u32.into());
        // Underflow: 5 - 10 = -5. -5 mod 17 = 12
        assert_eq!(field.sub(&5u32.into(), &10u32.into()), 12u32.into());
        // Large input: 20 - 5. (20 is 3 mod 17). 3 - 5 = -2. -2 mod 17 = 15
        assert_eq!(field.sub(&20u32.into(), &5u32.into()), 15u32.into());
    }

    #[test]
    fn test_additive_inverse() {
        let field = setup_field();
        assert_eq!(field.add_inv(&5u32.into()), 12u32.into());
        assert_eq!(field.add_inv(&0u32.into()), 0u32.into());
        assert_eq!(field.add_inv(&18u32.into()), 16u32.into());
    }

    #[test]
    fn test_multiplication() {
        let field = setup_field();
        assert_eq!(field.mul(&3u32.into(), &10u32.into()), 13u32.into());
        assert_eq!(field.mul(&0u32.into(), &5u32.into()), 0u32.into());
        assert_eq!(field.mul(&11u32.into(), &14u32.into()), 1u32.into());
    }

    #[test]
    fn test_multiplicative_inverse() {
        let field = setup_field();
        assert_eq!(field.mul_inv(&3u32.into()).unwrap(), 6u32.into());

        // Test Error: Zero
        assert_eq!(
            field.mul_inv(&0u32.into()).unwrap_err(),
            FiniteFieldError::ZeroHasNoInverse
        );

        // Test Error: No Inverse (using composite modulus)
        let composite_field = FiniteField { p: 10u32.into() };
        // gcd(2, 10) = 2, so 2 has no inverse mod 10
        assert_eq!(
            composite_field.mul_inv(&2u32.into()).unwrap_err(),
            FiniteFieldError::NoInverse
        );
    }

    #[test]
    fn test_division() {
        let field = setup_field();
        assert_eq!(field.div(&10u32.into(), &2u32.into()).unwrap(), 5u32.into());
        assert_eq!(field.div(&1u32.into(), &3u32.into()).unwrap(), 6u32.into());

        // Test Error: Divide by zero
        assert_eq!(
            field.div(&10u32.into(), &0u32.into()).unwrap_err(),
            FiniteFieldError::DivisionByZero
        );
    }

    #[test]
    fn test_power() {
        let field = setup_field();
        assert_eq!(field.pow(&4u32.into(), &2u32.into()), 16u32.into());
        assert_eq!(field.pow(&0u32.into(), &3u32.into()), 1u32.into());
        assert_eq!(field.pow(&5u32.into(), &0u32.into()), 0u32.into());
    }
}