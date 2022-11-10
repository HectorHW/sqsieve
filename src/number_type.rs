use crypto_bigint::prelude::*;
use crypto_bigint::{NonZero, U128, U256, U512, U64};
use num_bigint::BigUint;

pub trait NumberOps:
    Sized + crypto_bigint::Integer + crypto_bigint::Concat + crypto_bigint::Encoding
where
    Self: core::ops::ShrAssign<usize>,
{
    type NonzeroType = NonZero<Self>;

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

    fn wrapping_add(&self, other: &Self) -> Self;
    fn wrapping_sub(&self, other: &Self) -> Self;
    fn wrapping_div(&self, other: &Self) -> Self;
    fn wrapping_rem(&self, other: &Self) -> Self;
    fn wrapping_mul(&self, other: &Self) -> Self;
    fn bit_vartime(self, size: usize) -> u64;
    fn sqrt(&self) -> Self;

    fn size_in_bytes() -> usize;
}

macro_rules! impl_number_ops {
    ($t:ty, $size: expr) => {
        paste::item! {
            static [< $t _ONE >]: $t = <$t>::from_u32(1);
        }

        impl NumberOps for $t {
            #[inline(always)]
            fn one() -> &'static Self {
                &paste::item! { [< $t _ONE >] }
            }

            #[inline]
            fn convert_usize(value: usize) -> Self {
                <$t>::from_u64(value as u64)
            }

            #[inline]
            fn rem_short(&self, modulus: usize) -> Self {
                self.wrapping_rem(&Self::convert_usize(modulus))
            }

            #[inline]
            fn rem_long(&self, modulus: Self) -> Self {
                self.wrapping_rem(&modulus)
            }

            #[inline]
            fn divmod(&self, d: usize) -> (Self, Self) {
                self.div_rem(&Self::convert_usize(d)).unwrap()
            }

            #[inline]
            fn div_assign(&mut self, d: usize) {
                let data = *self;
                *self = data.wrapping_div(&Self::convert_usize(d))
            }

            #[inline]
            fn to_varsize(self) -> BigUint {
                let buf = self.to_be_bytes();
                BigUint::from_bytes_be(&buf)
            }

            #[inline]
            fn add_usize(self, other: usize) -> Self {
                self.wrapping_add(&Self::convert_usize(other))
            }

            #[inline]
            fn modpow2(&self, modulo: &Self) -> Self {
                //we want to compute self^2 mod m. We will use multiplication with upcast

                let mut buf: [u64; $size * 2] = [0; $size * 2];
                let modulo = modulo.to_words();

                buf[0..$size].copy_from_slice(&modulo);

                let modulo = <Self as Concat>::Output::from_words(buf);

                self.square().wrapping_rem(&modulo).split().1
            }

            fn to_usize(self) -> usize {
                self.to_words()[0] as usize
            }

            fn wrapping_add(&self, other: &Self) -> Self {
                <$t>::wrapping_add(self, other)
            }

            fn wrapping_sub(&self, other: &Self) -> Self {
                <$t>::wrapping_sub(self, other)
            }

            fn wrapping_div(&self, other: &Self) -> Self {
                <$t>::wrapping_div(self, other)
            }

            fn wrapping_mul(&self, other: &Self) -> Self {
                <$t>::wrapping_mul(self, other)
            }

            fn wrapping_rem(&self, other: &Self) -> Self {
                <$t>::wrapping_rem(self, other)
            }

            fn bit_vartime(self, value: usize) -> u64 {
                <$t>::bit_vartime(self, value)
            }

            fn sqrt(&self) -> Self {
                <$t>::sqrt(self)
            }

            fn size_in_bytes() -> usize {
                $size
            }
        }
    };
}

impl_number_ops!(U64, 1);
impl_number_ops!(U128, 2);
impl_number_ops!(U256, 4);
impl_number_ops!(U512, 8);

#[cfg(test)]
mod tests {
    use crypto_bigint::U512;

    use super::NumberOps;

    #[test]
    fn number_uses_little_endian() {
        let original_number = 153;
        assert_eq!(
            U512::convert_usize(original_number).to_usize(),
            original_number
        )
    }
}
