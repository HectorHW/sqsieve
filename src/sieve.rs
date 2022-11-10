use std::{error::Error, iter::repeat_with, ops::Rem, time::Instant};

use crypto_bigint::{subtle::ConditionallySelectable, Zero};
use itertools::Itertools;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{Pow, ToPrimitive};

use crate::{
    number_type::{NumberOps, NumberType},
    numbers::{factor_smooth, tonelli_shanks},
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

pub type SmoothiesVec = Vec<SmoothNumber<NumberType>>;

pub struct TestDivisionSieve {
    n: NumberType,
    factor_base: Vec<usize>,
    next_number: NumberType,
}

impl TestDivisionSieve {
    pub fn new(n: NumberType, factor_base: Vec<usize>) -> Self {
        let lower_bound = n.sqrt().add_usize(1);
        TestDivisionSieve {
            n,
            factor_base,
            next_number: lower_bound,
        }
    }

    pub fn run(&mut self, mut numbers_to_find: usize) -> SmoothiesVec {
        let mut result = vec![];

        let total_numbers = numbers_to_find;

        let mut last_time = Instant::now();

        while numbers_to_find > 0 {
            let sq_mod = self.next_number.clone().modpow2(&self.n);

            if let Some(mapping) = factor_smooth(&sq_mod, &self.factor_base) {
                #[cfg(feature = "verbose")]
                println!(
                    "found number {}^2 === {}",
                    self.next_number.to_varsize(),
                    sq_mod.to_varsize()
                );
                result.push(SmoothNumber {
                    number: self.next_number.clone(),
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

            self.next_number = self.next_number.add_usize(1);
        }

        result
    }
}

pub struct BlockSieve {
    n: NumberType,
    factor_base: Vec<usize>,
    roots: Vec<Option<(usize, usize)>>,
    block_size: usize,
    next_block: NumberType,
}

struct BlockEntry {
    original_number: NumberType,
    accumulator: NumberType,
    factorization: Vec<(usize, usize)>,
}

const BLOCK_SIZE: usize = 5000;

impl BlockSieve {
    pub fn new(n: NumberType, factor_base: Vec<usize>) -> Self {
        let roots = factor_base
            .iter()
            .map(|&factor| {
                if factor == 2 {
                    return None;
                }
                let n = n.to_varsize().rem(factor).to_u64().unwrap() as usize;

                let Some(s1) = tonelli_shanks(n, factor) else{
                    return None;
                };

                let s2 = factor - s1;
                Some((s1, s2))
            })
            .collect_vec();

        #[cfg(feature = "verbose")]
        {
            println!("factor base: {factor_base:?}");
            println!("roots: {roots:?}");
        }

        let block_size = usize::max(BLOCK_SIZE, factor_base.last().unwrap() * 2);

        BlockSieve {
            n,
            factor_base,
            roots,
            block_size,
            next_block: n.sqrt().add_usize(1),
        }
    }

    pub fn run(&mut self, total_numbers: usize) -> SmoothiesVec {
        let mut result = vec![];

        let total_numbers = total_numbers as isize;

        let mut numbers_to_find = total_numbers;

        let mut last_time = Instant::now();

        while numbers_to_find > 0 {
            let mut produced_items = self.search_block();

            #[cfg(feature = "verbose")]
            println!("block produced {} items", produced_items.len());

            self.next_block = self
                .next_block
                .wrapping_add(&NumberType::convert_usize(self.block_size));

            numbers_to_find -= produced_items.len() as isize;

            result.append(&mut produced_items);

            let now = Instant::now();

            if (now - last_time).as_secs() >= 5 {
                last_time = now;
                println!(
                    "done {:.1}%",
                    (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
                );
            }
        }
        result
    }

    fn search_block(&mut self) -> Vec<SmoothNumber<NumberType>> {
        #[cfg(feature = "verbose")]
        println!(
            "working with block size {} starting at {}",
            self.block_size,
            self.next_block.to_varsize()
        );

        let mut start = self.next_block;

        let mut block = repeat_with(|| {
            let number = start;
            start = start.wrapping_add(NumberType::one());

            BlockEntry {
                original_number: number,
                accumulator: number.modpow2(&self.n),
                factorization: Vec::with_capacity(16),
            }
        })
        .take(self.block_size)
        .collect_vec();

        //sieve for 2 manually, quick bit trick
        if self.factor_base[0] == 2 {
            let mut idx: usize = 0;

            if block[0].accumulator.bit_vartime(0) == 1 {
                idx += 1;
            }

            while idx < block.len() {
                let mut exp = 0;
                while block[idx].accumulator.bit_vartime(exp) == 0 {
                    exp += 1;
                }
                debug_assert!(exp > 0);

                block[idx].accumulator >>= exp;
                block[idx].factorization.push((2, exp));

                idx += 2;
            }
        }

        'primeiter: for (i, &prime) in self.factor_base.iter().enumerate() {
            //find start of sequence by trying different items

            let Some((s1, s2)) = self.roots[i] else{
                continue;
            };

            'rootiter: for root in [s1, s2] {
                //find closest value
                let mut idx = 0;
                loop {
                    let (d, r) = block[idx].original_number.divmod(prime);
                    if r == NumberType::convert_usize(root) {
                        break;
                    }
                    idx += 1;
                    if idx >= block.len() {
                        continue 'rootiter;
                    }
                }

                #[cfg(feature = "verbose")]
                println!("prime is {prime}, root is {root}, idx is {idx}");

                while idx < block.len() {
                    //dbg!(block[idx].accumulator.rem_short(prime).to_varsize());
                    //dbg!(idx);
                    let mut exponent = 0;
                    loop {
                        let (d, r) = block[idx].accumulator.divmod(prime);
                        if r != NumberType::from_u32(0) {
                            break;
                        }
                        exponent += 1;
                        block[idx].accumulator = d;
                    }

                    //dbg!(exponent);

                    block[idx].factorization.push((prime, exponent));
                    idx += prime;
                }
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
}

pub struct LogSieve {
    n: NumberType,
    factor_base: Vec<usize>,
    roots: Vec<Option<(usize, usize)>>,
    block_size: usize,
    next_block: NumberType,
    log_treshold: f64,
}

const LOGSIEVE_BLOCK_SIZE: usize = 15000;

impl LogSieve {
    //initalization is the same, what differs is the search algorithm
    pub fn new(n: NumberType, factor_base: Vec<usize>) -> Self {
        let roots = factor_base
            .iter()
            .map(|&factor| {
                let n = n.to_varsize().rem(factor).to_u64().unwrap() as usize;

                let Some(s1) = tonelli_shanks(n, factor) else{
                    return None;
                };

                let s2 = factor - s1;
                Some((s1, s2))
            })
            .collect_vec();

        #[cfg(feature = "verbose")]
        println!("roots: {roots:?}");

        let block_size = usize::max(LOGSIEVE_BLOCK_SIZE, factor_base.last().unwrap() * 2);

        LogSieve {
            n,
            factor_base,
            roots,
            block_size,
            next_block: n.sqrt().add_usize(1),
            log_treshold: n.to_varsize().to_f64().unwrap().log2() + (block_size as f64).log2(),
        }
    }

    pub fn run(&mut self, total_numbers: usize) -> SmoothiesVec {
        let mut result = vec![];

        let total_numbers = total_numbers as isize;

        let mut numbers_to_find = total_numbers;

        let mut last_time = Instant::now();

        let block_size = NumberType::convert_usize(self.block_size);

        while numbers_to_find > 0 {
            let mut produced_items = self.search_block();

            #[cfg(feature = "verbose")]
            println!("block produced {} items", produced_items.len());

            self.next_block = self.next_block.wrapping_add(&block_size);

            numbers_to_find -= produced_items.len() as isize;

            result.append(&mut produced_items);

            let now = Instant::now();

            if (now - last_time).as_secs() >= 5 {
                last_time = now;
                println!(
                    "done {:.1}%",
                    (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
                );
            }
        }
        result
    }

    fn search_block(&mut self) -> Vec<SmoothNumber<NumberType>> {
        #[cfg(feature = "verbose")]
        println!(
            "working with block size {} starting at {}",
            self.block_size,
            self.next_block.to_varsize()
        );

        let mut start = self.next_block;

        let (original_numbers, mut logs): (Vec<NumberType>, Vec<f64>) = repeat_with(|| {
            let number = start;
            start = start.wrapping_add(NumberType::one());

            let accumulator_log = number
                .modpow2(&self.n)
                .to_varsize()
                .to_f64()
                .unwrap()
                .log2();

            (number, accumulator_log)
        })
        .take(self.block_size)
        .unzip();

        'primeiter: for (i, &prime) in self.factor_base.iter().enumerate() {
            //find start of sequence by trying different items

            let Some((s1, s2)) = self.roots[i] else{
                continue;
            };

            'rootiter: for root in [s1, s2] {
                //find closest value
                let mut idx = 0;
                loop {
                    let (d, r) = original_numbers[idx].divmod(prime);
                    if r == NumberType::convert_usize(root) {
                        break;
                    }
                    idx += 1;
                    if idx >= original_numbers.len() {
                        continue 'rootiter;
                    }
                }

                let root_log_value = (prime as f64).log2();

                #[cfg(feature = "verbose")]
                println!("prime is {prime}, root is {root}, idx is {idx}");

                while idx < original_numbers.len() {
                    logs[idx] -= root_log_value;

                    idx += prime;
                }
            }
        }

        let mut buf = logs.clone();

        let dyn_treshold = *buf
            .select_nth_unstable_by(self.block_size / 1000, |a, b| a.partial_cmp(b).unwrap())
            .1;

        original_numbers
            .into_iter()
            .zip(logs.into_iter())
            .filter_map(|(n, ln)| if ln <= dyn_treshold { Some(n) } else { None })
            .filter_map(|n| {
                let Some(divisors) = factor_smooth(&n, &self.factor_base) else {
                return None;
            };

                Some(SmoothNumber {
                    number: n,
                    divisors,
                })
            })
            .collect_vec()
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Rem;

    use num_bigint::BigInt;

    #[test]
    fn negative_power_building() {
        let n = BigInt::from(-5);
        assert_eq!(
            n.modpow(&BigInt::from(1), &BigInt::from(8)),
            BigInt::from(3)
        );
    }
}
