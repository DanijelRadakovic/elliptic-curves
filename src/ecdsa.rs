use crate::elliptic_curve::{EllipticCurve, EllipticCurveError, Point};
use crate::finite_field::FiniteField;
use num_bigint::{BigRng010, BigUint};
use num_traits::{One, Zero};
use rand::rng;
use thiserror::Error;

#[derive(Debug)]
pub struct ECDSA {
    // Elliptic Curve used for DSA
    ec: EllipticCurve,
    // Generator of the EC
    g: Point,
    // Order of the subgroup
    q: BigUint,
}

#[derive(Error, Debug, PartialEq)]
pub enum ECDSAErrors {
    #[error("Provided generator is not on the curve")]
    InvalidGenerator,
}

impl ECDSA {
    pub fn new(ec: EllipticCurve, g: Point, q: BigUint) -> Result<Self, ECDSAErrors> {
        if !ec.is_on_curve(&g) {
            return Err(ECDSAErrors::InvalidGenerator)?;
        }
        Ok(ECDSA { ec, g, q })
    }

    pub fn generate_key_pair(&self) -> Result<(BigUint, Point), EllipticCurveError> {
        let mut rng = rng();
        let d = rng.random_biguint_range(&BigUint::one(), &self.q);
        let pub_key = self.ec.scalar_mul(&d, &self.g)?;
        Ok((d, pub_key))
    }

    pub fn sign(
        &self,
        digest: &BigUint,
        private_key: &BigUint,
    ) -> Result<(BigUint, BigUint), EllipticCurveError> {
        let mut rng = rng();
        let f = FiniteField { p: self.q.clone() };

        loop {
            let k = rng.random_biguint_range(&BigUint::one(), &self.q);
            let r_point = self.ec.scalar_mul(&k, &self.g)?;

            if let Point::Coordinate { x, y: _ } = r_point {
                let r = x % self.q.clone();
                if r == BigUint::zero() {
                    continue;
                }

                let k_inv = f.mul_inv(&k)?;
                let s = f.mul(&f.add(&f.mul(private_key, &r), digest), &k_inv);
                if s.is_zero() {
                    continue;
                }
                return Ok((r, s));
            }
        }
    }

    pub fn verify(
        &self,
        digest: &BigUint,
        sig: &(BigUint, BigUint),
        public_key: &Point,
    ) -> Result<bool, EllipticCurveError> {
        let (r, s) = sig;
        let f = FiniteField { p: self.q.clone() };
        let w = f.mul_inv(s)?;
        let u1 = f.mul(digest, &w);
        let u2 = f.mul(&w, r);
        let point = self.ec.add(
            &self.ec.scalar_mul(&u1, &self.g)?,
            &self.ec.scalar_mul(&u2, public_key)?,
        )?;
        match point {
            Point::Coordinate { x, .. } => Ok((x % &self.q) == *r),
            Point::Identity => Ok(false),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Num;

    fn setup_ecdsa() -> ECDSA {
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into()).unwrap();
        ECDSA::new(
            ec,
            Point::Coordinate {
                x: 5u32.into(),
                y: 1u32.into(),
            },
            19u32.into(),
        )
        .unwrap()
    }

    #[test]
    fn test_invalid_generator() {
        let ec = EllipticCurve::new(2u32.into(), 2u32.into(), 17u32.into()).unwrap();
        let res = ECDSA::new(
            ec,
            Point::Coordinate {
                x: 3u32.into(),
                y: 3u32.into(),
            },
            19u32.into(),
        );
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), ECDSAErrors::InvalidGenerator);
    }

    #[test]
    fn test_generate_key_pair() {
        let ecdsa = setup_ecdsa();
        let (private_key, public_key) = ecdsa.generate_key_pair().unwrap();
        assert!(private_key.lt(&19u32.into()) && private_key.gt(&0u32.into()));
        assert!(ecdsa.ec.is_on_curve(&public_key));
    }

    #[test]
    fn test_sign_and_verify() {
        let ecdsa = setup_ecdsa();
        let (private_key, public_key) = ecdsa.generate_key_pair().unwrap();
        let msg = b"For the Emperor!";
        let digest = BigUint::from_str_radix(&sha256::digest(msg), 16).unwrap();
        let signature = ecdsa.sign(&digest, &private_key).unwrap();
        let result = ecdsa.verify(&digest, &signature, &public_key).unwrap();
        assert!(!signature.0.is_zero());
        assert!(!signature.1.is_zero());
        assert_eq!(result, true);
    }

    #[test]
    fn test_invalid_signature() {
        let ecdsa = setup_ecdsa();
        let (private_key, public_key) = ecdsa.generate_key_pair().unwrap();
        let msg = b"For the Emperor!";
        let digest = BigUint::from_str_radix(&sha256::digest(msg), 16).unwrap();
        let signature = ecdsa.sign(&digest, &private_key).unwrap();

        let invalid_msg = b"Heresy!";
        let invalid_digest = BigUint::from_str_radix(&sha256::digest(invalid_msg), 16).unwrap();
        let result = ecdsa.verify(&invalid_digest, &signature, &public_key).unwrap();
        assert_eq!(result, false);
    }
}
