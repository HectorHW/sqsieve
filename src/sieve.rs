use std::{error::Error, iter::repeat_with, time::Instant};

use crypto_bigint::Zero;
use itertools::Itertools;
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

pub fn test_division_sieve<F: FnOnce(&NumberType, &[usize]) -> usize>(
    n: &NumberType,
    prime_table: &[usize],
    stop_criteria: F,
    continuation: Option<NumberType>,
) -> Result<(SmoothiesVec, NumberType), Box<dyn Error>> {
    let lower_bound = continuation.unwrap_or_else(|| n.sqrt().add_usize(1));

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

    Ok((result, next_number.wrapping_add(NumberType::one())))
}

const BLOCK_SIZE: usize = 5000;

pub fn block_division_sieve<F: FnOnce(&NumberType, &[usize]) -> usize>(
    n: &NumberType,
    prime_table: &[usize],
    stop_criteria: F,
) -> Result<SmoothiesVec, Box<dyn Error>> {
    let lower_bound = n.sqrt().add_usize(1);

    let mut result = vec![];

    let mut next_number = lower_bound;

    let total_numbers = stop_criteria(n, prime_table);

    let mut numbers_to_find = total_numbers as isize;

    let mut current_bound = lower_bound;

    let step = NumberType::convert_usize(BLOCK_SIZE);

    let mut last_time = Instant::now();

    while numbers_to_find > 0 {
        let mut produced_items = search_block(n, BLOCK_SIZE, current_bound, prime_table);

        #[cfg(feature = "verbose")]
        println!("block produced {} items", produced_items.len());

        current_bound = current_bound.wrapping_add(&step);

        numbers_to_find -= produced_items.len() as isize;

        result.append(&mut produced_items);

        let now = Instant::now();

        if (now - last_time).as_secs() >= 5 {
            last_time = now;
            println!(
                "done {:.1}%",
                (total_numbers as isize - numbers_to_find) as f64 / total_numbers as f64 * 100f64
            );
        }
    }
    Ok(result)
}

struct BlockEntry {
    original_number: NumberType,
    accumulator: NumberType,
    factorization: Vec<(usize, usize)>,
}

fn search_block(
    n: &NumberType,
    block_size: usize,
    mut start: NumberType,
    prime_table: &[usize],
) -> Vec<SmoothNumber<NumberType>> {
    assert!(block_size > prime_table.last().cloned().unwrap());
    // #[cfg(feature = "verbose")]
    println!(
        "working with block size {} starting at {}",
        block_size,
        start.to_varsize()
    );

    let mut block = repeat_with(|| {
        let number = start.clone();
        start = start.wrapping_add(NumberType::one());

        BlockEntry {
            original_number: number,
            accumulator: number.modpow2(n),
            factorization: Vec::with_capacity(16),
        }
    })
    .take(block_size)
    .collect_vec();

    'primeiter: for &prime in prime_table {
        //find start of sequence by trying different items
        let mut idx = 0;

        loop {
            unimplemented!("this search is wrong, use tonelli-shanks to find roots");
            let (d, r) = block[idx].accumulator.divmod(prime);
            if r.is_zero().unwrap_u8() == 1 {
                break;
            }
            idx += 1;
            if idx >= block.len() {
                continue 'primeiter;
            }
        }

        println!("prime is {prime}, idx is {idx}");

        while idx < block.len() {
            dbg!(block[idx].accumulator.rem_short(prime).to_varsize());
            dbg!(idx);
            let mut exponent = 0;
            loop {
                let (d, r) = block[idx].accumulator.divmod(prime);
                if r != NumberType::from_u32(0) {
                    break;
                }
                exponent += 1;
                block[idx].accumulator = d;
            }

            dbg!(exponent);

            block[idx].factorization.push((prime, exponent));
            idx += prime;
        }
    }

    block
        .into_iter()
        .filter_map(|item| {
            if &item.accumulator != NumberType::one() {
                return None;
            }

            #[cfg(feature = "verbose")]
            println!("found number {}", item.original_number.to_varsize());

            Some(SmoothNumber {
                number: item.original_number,
                divisors: item.factorization,
            })
        })
        .collect_vec()
}
