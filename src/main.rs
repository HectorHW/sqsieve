#![feature(associated_type_defaults)]

use std::{env::args, error::Error, iter::repeat, str::FromStr};

use crypto_bigint::UInt;
use factor_building::find_factors_from_pivots;
use factorization::factorize;
use itertools::Itertools;
use num_bigint::BigUint;
use number_type::NumberOps;
use sieve::SmoothNumber;

use crate::{
    factor_building::{find_factor_exhaustive, find_factor_simple, find_factors_random},
    numbers::{build_factor_base, small_eratosphenes},
    sieve::{compute_b_limit, BlockSieve, LogSieve, TestDivisionSieve},
    solver::{produce_solution, CongruenceSystem},
};

mod number_type;
mod numbers;

mod factor_building;
mod factorization;
mod sieve;
mod solver;

fn main() {
    //let n = BigUint::from_str("1577271624417732056618338337651").unwrap();

    let args = args().collect::<Vec<String>>();

    if args.len() != 2 {
        println!("ERROR: please provide number as single argument");
        return;
    }

    match factorize(args[1].clone()) {
        Ok(r) => {
            println!("SUCCESS: {r:?}");
        }

        Err(e) => {
            println!("ERROR: {e}");
        }
    }
}
