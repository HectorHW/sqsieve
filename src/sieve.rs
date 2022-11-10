use std::{iter::repeat_with, ops::Rem, time::Instant};

use itertools::Itertools;

use num_traits::{Pow, ToPrimitive};

use crate::{
    number_type::NumberOps,
    numbers::{tonelli_shanks, trial_divide},
};

pub fn compute_b_limit<NT: NumberOps>(n: &NT) -> usize {
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

pub type SmoothiesVec<NT> = Vec<SmoothNumber<NT>>;

pub struct TestDivisionSieve<NT: NumberOps> {
    n: NT,
    factor_base: Vec<usize>,
    next_number: NT,
}

impl<NT: NumberOps> TestDivisionSieve<NT> {
    pub fn new(n: NT, factor_base: Vec<usize>) -> Self {
        let lower_bound = n.sqrt().add_usize(1);
        TestDivisionSieve {
            n,
            factor_base,
            next_number: lower_bound,
        }
    }

    pub fn run(&mut self, mut numbers_to_find: usize) -> SmoothiesVec<NT> {
        let mut result = vec![];

        let total_numbers = numbers_to_find;

        let mut last_time = Instant::now();

        while numbers_to_find > 0 {
            let sq_mod = self.next_number.clone().modpow2(&self.n);

            if let Some(mapping) = trial_divide(&sq_mod, &self.factor_base) {
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

pub struct BlockSieve<NT: NumberOps> {
    n: NT,
    factor_base: Vec<usize>,
    roots: Vec<Option<(usize, usize)>>,
    block_size: usize,
    next_block: NT,
}

struct BlockEntry<NT: NumberOps> {
    original_number: NT,
    accumulator: NT,
    factorization: Vec<(usize, usize)>,
}

const BLOCK_SIZE: usize = 5000;

impl<NT: NumberOps> BlockSieve<NT> {
    pub fn new(n: NT, factor_base: Vec<usize>) -> Self {
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

        let block_size = usize::max(BLOCK_SIZE, factor_base.last().unwrap() * 5);

        BlockSieve {
            n,
            factor_base,
            roots,
            block_size,
            next_block: n.sqrt().add_usize(1),
        }
    }

    pub fn run(&mut self, total_numbers: usize) -> SmoothiesVec<NT> {
        let mut result = vec![];

        let total_numbers = total_numbers as isize;

        let mut numbers_to_find = total_numbers;

        let mut last_time = Instant::now();

        let long_block_size = NT::convert_usize(self.block_size);

        while numbers_to_find > 0 {
            let mut produced_items = self.search_block();

            #[cfg(feature = "verbose")]
            println!(
                "block of size {} produced {} items",
                self.block_size,
                produced_items.len()
            );

            self.next_block = self.next_block.wrapping_add(&long_block_size);

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

    fn search_block(&mut self) -> Vec<SmoothNumber<NT>> {
        #[cfg(feature = "verbose")]
        println!(
            "working with block size {} starting at {}",
            self.block_size,
            self.next_block.to_varsize()
        );

        let mut start = self.next_block;

        let mut block = repeat_with(|| {
            let number = start;
            start = start.wrapping_add(NT::one());

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

        for (i, &prime) in self.factor_base.iter().enumerate() {
            //find start of sequence by trying different items

            let Some((s1, s2)) = self.roots[i] else{
                continue;
            };

            for root in [s1, s2] {
                //find closest value
                let mut idx;

                let long_root = NT::convert_usize(root);
                let long_prime = NT::convert_usize(prime);

                let mut closest_element = (block[0].original_number.wrapping_sub(&long_root))
                    .wrapping_div(&long_prime)
                    .wrapping_mul(&long_prime)
                    .wrapping_add(&long_root);

                if closest_element < block[0].original_number {
                    closest_element = closest_element.wrapping_add(&long_prime);
                }

                idx = closest_element
                    .wrapping_sub(&block[0].original_number)
                    .to_usize();

                debug_assert!({
                    let (_, r) = block[idx].original_number.divmod(prime);
                    r == NT::convert_usize(root)
                });

                #[cfg(feature = "verbose")]
                println!("prime is {prime}, root is {root}, idx is {idx}");

                while idx < block.len() {
                    //dbg!(block[idx].accumulator.rem_short(prime).to_varsize());
                    //dbg!(idx);
                    let mut exponent = 0;
                    loop {
                        let (d, r) = block[idx].accumulator.divmod(prime);
                        if r != NT::convert_usize(0) {
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
                if &item.accumulator != NT::one() {
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

pub struct LogSieve<NT: NumberOps> {
    n: NT,
    factor_base: Vec<usize>,
    roots: Vec<Option<(usize, usize)>>,
    block_size: usize,
    next_block: NT,
    log_treshold: f64,
}

const LOGSIEVE_BLOCK_SIZE: usize = 60_000;

impl<NT: NumberOps> LogSieve<NT> {
    //initalization is the same, what differs is the search algorithm
    pub fn new(n: NT, factor_base: Vec<usize>) -> Self {
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

        let block_size = usize::max(
            usize::min(LOGSIEVE_BLOCK_SIZE, factor_base.last().unwrap() * 5),
            factor_base.last().unwrap() * 2,
        );

        let n_f = n.to_varsize().to_f64().unwrap();

        let treshold = (block_size as f64).log2() + n_f.log2() * 0.5
            - Self::chose_t(n_f.log10()) * (*factor_base.last().unwrap() as f64).log2();

        LogSieve {
            n,
            factor_base,
            roots,
            block_size,
            next_block: n.sqrt().add_usize(1),
            log_treshold: treshold,
        }
    }

    fn chose_t(number_size: f64) -> f64 {
        if number_size <= 30.0 {
            return 1.5;
        }
        if number_size <= 45.0 {
            return 2.0;
        }
        if number_size <= 66.0 {
            return 2.6;
        }

        3.2
    }

    pub fn run(&mut self, total_numbers: usize) -> SmoothiesVec<NT> {
        println!("running log sieve with block size of {}", self.block_size);
        let mut result = vec![];

        let total_numbers = total_numbers as isize;

        let mut numbers_to_find = total_numbers;

        let mut last_time = Instant::now();

        let block_size = NT::convert_usize(self.block_size);

        while numbers_to_find > 0 {
            let mut produced_items = self.search_block(self.next_block, self.block_size);

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

        println!(
            "done {:.1}%",
            (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
        );

        result
    }

    fn search_block(&self, mut start: NT, size: usize) -> Vec<SmoothNumber<NT>> {
        #[cfg(feature = "verbose")]
        println!(
            "working with block size {} starting at {}",
            size,
            self.next_block.to_varsize()
        );

        let (original_numbers, mut accumulators): (Vec<NT>, Vec<NT>) = repeat_with(|| {
            let number = start;
            let accumulator = number.modpow2(&self.n);
            start = start.wrapping_add(NT::one());
            (number, accumulator)
        })
        .take(size)
        .unzip();

        let mut logs = vec![0f64; size];

        let acc2 = accumulators.clone();

        //sieve for 2 manually, quick bit trick
        if self.factor_base[0] == 2 {
            let mut idx: usize = 0;

            let n0 = accumulators[0];

            if n0.bit_vartime(0) == 1 {
                idx += 1;
            }

            while idx < logs.len() {
                let mut exp = 0;
                while accumulators[idx].bit_vartime(exp) == 0 {
                    exp += 1;
                }
                debug_assert!(exp > 0);

                logs[idx] += exp as f64;
                accumulators[idx] >>= exp;

                idx += 2;
            }
        }

        for (i, &prime) in self.factor_base.iter().enumerate() {
            //find start of sequence by trying different items

            let Some((s1, s2)) = self.roots[i] else{
                continue;
            };

            for root in [s1, s2] {
                //find closest value
                let mut idx = 0;

                let long_root = NT::convert_usize(root);
                let long_prime = NT::convert_usize(prime);

                let mut closest_element = (original_numbers[idx].wrapping_sub(&long_root))
                    .wrapping_div(&long_prime)
                    .wrapping_mul(&long_prime)
                    .wrapping_add(&long_root);

                if closest_element < original_numbers[0] {
                    closest_element = closest_element.wrapping_add(&long_prime);
                }

                idx = closest_element
                    .wrapping_sub(&original_numbers[0])
                    .to_usize();

                debug_assert!({
                    let (_, r) = original_numbers[idx].divmod(prime);
                    r == NT::convert_usize(root)
                });

                let root_log_value = (prime as f64).log2();

                #[cfg(feature = "verbose")]
                println!("prime is {prime}, root is {root}, idx is {idx}");

                while idx < original_numbers.len() {
                    logs[idx] += root_log_value;

                    idx += prime;
                }
            }
        }

        original_numbers
            .into_iter()
            .zip(logs.into_iter())
            .zip(acc2.into_iter())
            .filter_map(|((n, ln), acc)| {
                if ln >= self.log_treshold {
                    Some((n, acc))
                } else {
                    None
                }
            })
            .filter_map(|(number, acc)| {
                let Some(divisors) = trial_divide(&acc, &self.factor_base) else {
                return None;
            };

                Some(SmoothNumber { number, divisors })
            })
            .collect_vec()
    }

    pub fn run_parallel(&mut self, total_numbers: usize) -> Vec<SmoothNumber<NT>> {
        use rayon::prelude::*;
        println!("running parallel log sieve");
        let mut result = vec![];

        let total_numbers = total_numbers as isize;

        let mut numbers_to_find = total_numbers;

        let mut last_time = Instant::now();

        use std::env;

        let threads = env::var("THREADS")
            .map(|n| n.parse::<usize>().expect("failed to parse thread count"))
            .unwrap_or(1);
        let block_size = self.factor_base.last().cloned().unwrap() * 2;

        println!("block size is {block_size}, {threads} threads");

        while numbers_to_find > 0 {
            let jobs = repeat_with(|| {
                let block_start = self.next_block;

                self.next_block = self.next_block.add_usize(block_size);

                (block_start, block_size)
            })
            .take(threads)
            .collect_vec();

            let items = jobs
                .into_par_iter()
                .map(|(start, size)| self.search_block(start, size))
                .collect::<Vec<_>>();

            #[cfg(feature = "verbose")]
            println!("block produced {} items", produced_items.len());

            numbers_to_find -= items.iter().map(|v| v.len()).sum::<usize>() as isize;

            for mut item in items {
                result.append(&mut item);
            }

            let now = Instant::now();

            if (now - last_time).as_secs() >= 5 {
                last_time = now;
                println!(
                    "done {:.1}%",
                    (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
                );
            }
        }

        println!(
            "done {:.1}%",
            (total_numbers - numbers_to_find) as f64 / total_numbers as f64 * 100f64
        );

        result
    }
}

#[cfg(test)]
mod tests {

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
