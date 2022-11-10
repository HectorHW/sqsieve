#![feature(associated_type_defaults)]

use std::env::args;

use factorization::factorize;

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
