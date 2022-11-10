use core::ops::Rem;

use crypto_bigint::{prelude::*, Wrapping, U1024};
use crypto_bigint::{NonZero, U512};
use num_bigint::BigUint;

use std::ops::{Add, Div};

pub trait NumberOps: Sized {
    fn one() -> &'static Self;
    fn modpow2(&self, modulo: &Self) -> Self;

    fn convert_usize(value: usize) -> Self;

    fn rem_short(&self, modulus: usize) -> Self;
    fn rem_long(&self, modulus: Self) -> Self;

    fn divmod(&self, d: usize) -> (Self, Self);

    fn div_assign(&mut self, d: usize);

    fn add_usize(self, other: usize) -> Self;

    fn to_varsize(self) -> BigUint;

    fn to_usize(self) -> usize;
}

pub type NumberType = U512;
pub type NonzeroType = NonZero<U512>;

static ONE: NumberType = NumberType::from_u32(1);

impl NumberOps for NumberType {
    #[inline(always)]
    fn one() -> &'static Self {
        &ONE
    }

    #[inline(always)]
    fn convert_usize(value: usize) -> Self {
        NumberType::from_u64(value as u64)
    }

    #[inline]
    fn rem_short(&self, modulus: usize) -> Self {
        self.rem(NonzeroType::new(Self::convert_usize(modulus)).unwrap())
    }

    #[inline]
    fn rem_long(&self, modulus: Self) -> Self {
        self.rem(NonzeroType::new(modulus).unwrap())
    }

    #[inline]
    fn divmod(&self, d: usize) -> (Self, Self) {
        self.div_rem(&Self::convert_usize(d)).unwrap()
    }

    #[inline]
    fn div_assign(&mut self, d: usize) {
        let data = *self;
        *self = Wrapping(data)
            .div(NonzeroType::new(Self::convert_usize(d)).unwrap())
            .0;
    }

    #[inline]
    fn to_varsize(self) -> BigUint {
        let buf = self.to_be_bytes();
        BigUint::from_bytes_be(&buf)
    }

    #[inline]
    fn add_usize(self, other: usize) -> Self {
        let data = Wrapping(self);
        let other = Wrapping(Self::convert_usize(other));
        data.add(&other).0
    }

    #[inline]
    fn modpow2(&self, modulo: &Self) -> Self {
        //we want to compute self^2 mod m. We will use multiplication with upcast

        let mut buf: [u64; 16] = [0; 16];
        let modulo = modulo.to_words();

        buf[0..8].copy_from_slice(&modulo);

        let modulo = U1024::from_words(buf);

        self.square().wrapping_rem(&modulo).split().1
    }

    fn to_usize(self) -> usize {
        self.to_words()[0] as usize
    }
}

#[cfg(test)]
mod tests {
    use super::{NumberOps, NumberType};

    #[test]
    fn number_uses_little_endian() {
        let original_number = 153;
        assert_eq!(
            NumberType::convert_usize(original_number).to_usize(),
            original_number
        )
    }
}
