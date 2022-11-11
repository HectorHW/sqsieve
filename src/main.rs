#![feature(associated_type_defaults)]

use std::env::args;

use factorization::factorize;
use itertools::Itertools;
use num_bigint::BigUint;

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
        Ok(items) => {
            let prod: BigUint = items.iter().product();
            let eq = items.iter().map(|n| n.to_string()).join(" * ");
            println!("{eq} = {prod}");

            println!("SUCCESS: {items:?}");
        }

        Err(e) => {
            println!("ERROR: {e}");
        }
    }
}
