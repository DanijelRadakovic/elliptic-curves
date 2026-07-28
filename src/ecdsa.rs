use crate::elliptic_curve::{EllipticCurve, EllipticCurveError, Point};
use crate::finite_field::{FiniteField, FiniteFieldError};
use num_bigint::{BigRng010, BigUint};
use num_traits::{One, Zero};
use rand::rng;
use thiserror::Error;

#[derive(Debug)]
pub struct ECDSA {
    /// Elliptic Curve used for DSA
    ec: EllipticCurve,
    /// Generator of the EC
    g: Point,
    /// Order of the subgroup
    q: BigUint,
}

#[derive(Debug)]
pub struct Signature {
    pub r: BigUint,
    pub s: BigUint,
}

#[derive(Error, Debug, PartialEq)]
pub enum ECDSAError {
    #[error("Provided generator is not on the curve")]
    InvalidGenerator,
    #[error("Invalid private key")]
    InvalidPrivateKey,
    #[error("Invalid public key")]
    InvalidPublicKey,
    #[error(transparent)]
    EllipticCurveError(#[from] EllipticCurveError),
    #[error(transparent)]
    FiniteFieldError(#[from] FiniteFieldError),
}

impl ECDSA {
    pub fn new(ec: EllipticCurve, g: Point, q: BigUint) -> Result<Self, ECDSAError> {
        if !ec.is_on_curve(&g) {
            return Err(ECDSAError::InvalidGenerator)?;
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
    ) -> Result<Signature, ECDSAError> {
        if private_key.is_zero() || private_key >= &self.q {
            return Err(ECDSAError::InvalidPrivateKey);
        }

        // Per SEC 1: If the hash is longer than the order bits, take the leftmost bits.
        let q_bits = self.q.bits();
        let truncated_digest = if digest.bits() > q_bits {
            digest >> (digest.bits() - q_bits)
        } else {
            digest.clone()
        };
        let truncated_digest = truncated_digest % &self.q;

        let mut rng = rng();
        let f = FiniteField { p: self.q.clone() };

        loop {
            let k = rng.random_biguint_range(&BigUint::one(), &self.q);
            let r_point = self.ec.scalar_mul(&k, &self.g)?;

            if let Point::Coordinate { x, y: _ } = r_point {
                let r = x % &self.q;
                if r == BigUint::zero() {
                    continue;
                }

                let k_inv = f.mul_inv(&k)?;
                let s = f.mul(&f.add(&f.mul(private_key, &r), &truncated_digest), &k_inv);
                if s.is_zero() {
                    continue;
                }
                return Ok(Signature { r, s });
            }
        }
    }

    pub fn verify(
        &self,
        digest: &BigUint,
        sig: &Signature,
        public_key: &Point,
    ) -> Result<bool, ECDSAError> {
        self.validate_public_key(public_key)?;

        let Signature { r, s } = sig;
        if r.is_zero() || r >= &self.q || s.is_zero() || s >= &self.q {
            return Ok(false);
        }

        // Per SEC 1: If the hash is longer than the order bits, take the leftmost bits.
        let q_bits = self.q.bits();
        let truncated_digest = if digest.bits() > q_bits {
            digest >> (digest.bits() - q_bits)
        } else {
            digest.clone()
        };
        let truncated_digest = truncated_digest % &self.q;

        let f = FiniteField { p: self.q.clone() };
        let w = f.mul_inv(s)?;
        let u1 = f.mul(&truncated_digest, &w);
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

    fn validate_public_key(&self, public_key: &Point) -> Result<(), ECDSAError> {
        match public_key {
            Point::Coordinate { x, y } => {
                if x >= &self.ec.p() || y >= &self.ec.p() {
                    return Err(ECDSAError::InvalidPublicKey);
                }
                if !self.ec.is_on_curve(public_key) {
                    return Err(EllipticCurveError::PointNotOnCurve {
                        x: x.clone(),
                        y: y.clone(),
                        a: self.ec.a().clone(),
                        b: self.ec.b().clone(),
                        p: self.ec.p().clone(),
                    }
                        .into());
                }
                Ok(())
            }
            Point::Identity => Err(ECDSAError::InvalidPublicKey),
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
        assert_eq!(res.unwrap_err(), ECDSAError::InvalidGenerator);
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
        assert!(!signature.r.is_zero());
        assert!(!signature.s.is_zero());
        assert_eq!(ecdsa.verify(&digest, &signature, &public_key).unwrap(), true);
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
        assert_eq!(ecdsa.verify(&invalid_digest, &signature, &public_key).unwrap(), false);
    }

    #[test]
    fn test_verify_invalid_public_key_bounds() {
        let ecdsa = setup_ecdsa();
        let digest = BigUint::from(10u32);
        let sig = Signature { r: BigUint::from(1u32), s: BigUint::from(1u32) };

        let invalid_pub = Point::Coordinate {
            x: BigUint::from(18u32),
            y: BigUint::from(5u32),
        };

        let result = ecdsa.verify(&digest, &sig, &invalid_pub);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), ECDSAError::InvalidPublicKey);
    }

    #[test]
    fn test_verify_invalid_public_key_not_on_curve() {
        let ecdsa = setup_ecdsa();
        let digest = BigUint::from(10u32);
        let sig = Signature { r: BigUint::from(1u32), s: BigUint::from(1u32) };

        let invalid_pub = Point::Coordinate {
            x: BigUint::from(1u32),
            y: BigUint::from(1u32),
        };

        let result = ecdsa.verify(&digest, &sig, &invalid_pub);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), ECDSAError::EllipticCurveError(EllipticCurveError::PointNotOnCurve {
            x: 1u32.into(),
            y: 1u32.into(),
            a: 2u32.into(),
            b: 2u32.into(),
            p: 17u32.into(),
        }));
    }

    #[test]
    fn test_verify_invalid_signature_values() {
        let ecdsa = setup_ecdsa();
        let digest = BigUint::from(10u32);
        let (_, public_key) = ecdsa.generate_key_pair().unwrap();

        // q = 19.
        // r is zero
        let sig_r_zero = Signature { r: BigUint::zero(), s: BigUint::from(1u32) };
        assert_eq!(ecdsa.verify(&digest, &sig_r_zero, &public_key).unwrap(), false);

        // s is zero
        let sig_s_zero = Signature { r: BigUint::from(1u32), s: BigUint::zero() };
        assert_eq!(ecdsa.verify(&digest, &sig_s_zero, &public_key).unwrap(), false);

        // r is equal to order q (19)
        let sig_r_large = Signature { r: BigUint::from(19u32), s: BigUint::from(1u32) };
        assert_eq!(ecdsa.verify(&digest, &sig_r_large, &public_key).unwrap(), false);

        // s is greater than order q (20)
        let sig_s_large = Signature { r: BigUint::from(1u32), s: BigUint::from(20u32) };
        assert_eq!(ecdsa.verify(&digest, &sig_s_large, &public_key).unwrap(), false);
    }

    #[test]
    fn test_sign_invalid_private_key() {
        let ecdsa = setup_ecdsa();
        let digest = BigUint::from(10u32);

        // Private key is zero
        let res_zero = ecdsa.sign(&digest, &BigUint::zero());
        assert!(res_zero.is_err());
        assert_eq!(res_zero.unwrap_err(), ECDSAError::InvalidPrivateKey);

        // Private key is equal to order q
        let res_large = ecdsa.sign(&digest, &BigUint::from(19u32));
        assert!(res_large.is_err());
        assert_eq!(res_large.unwrap_err(), ECDSAError::InvalidPrivateKey);
    }

    #[test]
    fn test_verify_identity_public_key() {
        let ecdsa = setup_ecdsa();
        let digest = BigUint::from(10u32);
        let sig = Signature { r: BigUint::from(1u32), s: BigUint::from(1u32) };

        // Identity point cannot be a public key
        let result = ecdsa.verify(&digest, &sig, &Point::Identity);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), ECDSAError::InvalidPublicKey);
    }
}
