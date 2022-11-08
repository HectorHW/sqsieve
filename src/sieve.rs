use std::{error::Error, time::Instant};

use num_bigint::BigUint;
use num_traits::{Pow, ToPrimitive};

use crate::{
    number_type::{NumberOps, NumberType},
    numbers::factor_smooth,
};

pub fn compute_b_limit(n: &NumberType) -> usize {
    // as we are expected to work with rather small numbers (below 100 digits),
    // we can simply cast it to float to get somewhat decent result
    // (bc f64 will support up to 10^308)
    let float = n.to_varsize().to_f64().unwrap();
    let l = f64::exp(f64::sqrt(float.ln() * float.ln().ln()));
    let b = l.pow(1f64 / 2f64.sqrt());
    b.ceil() as usize
}

#[derive(Clone, Debug)]
pub struct SmoothNumber<NT> {
    pub number: NT,
    pub divisors: Vec<(usize, usize)>,
}

pub fn constant_b_limit<const N: usize>(_n: &NumberType) -> usize {
    N
}

/// number of smooth numbers we need to find while sieving
fn approximate_relations(_n: &NumberType, primes: &[usize]) -> usize {
    primes.len() + 7
}

pub type SmoothiesVec = Vec<SmoothNumber<NumberType>>;

pub fn sieve_for_smoothies<F: FnOnce(&NumberType, &[usize]) -> usize>(
    n: &NumberType,
    prime_table: &[usize],
    stop_criteria: F,
) -> Result<SmoothiesVec, Box<dyn Error>> {
    let lower_bound = n.sqrt().add_usize(1);

    let mut result = vec![];

    let mut next_number = lower_bound;

    let total_numbers = stop_criteria(n, prime_table);

    let mut numbers_to_find = total_numbers;

    let mut last_time = Instant::now();

    while numbers_to_find > 0 {
        let sq_mod = next_number.clone().modpow2(n);

        if let Some(mapping) = factor_smooth(&sq_mod, prime_table) {
            #[cfg(feature = "verbose")]
            println!(
                "found number {}^2 === {}",
                next_number.to_varsize(),
                sq_mod.to_varsize()
            );
            result.push(SmoothNumber {
                number: next_number.clone(),
                divisors: mapping,
            });
            numbers_to_find -= 1;

            let now = Instant::now();

            if (now - last_time).as_secs() >= 5 {
                last_time = now;
                println!(
                    "done {:.1}%",
                    (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
                );
            }
        }

        next_number = next_number.add_usize(1);
    }

    Ok(result)
}
