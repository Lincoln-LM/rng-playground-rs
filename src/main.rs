mod berkowitz;
mod gf2int;
mod gf2vec;
mod mat_builder;
mod pohlig_hellman;
mod rng;
mod xoroshiro;

use crate::xoroshiro::Xoroshiro128Plus;
use rand::{thread_rng, Rng};
use rng::RNG;
use std::time::Instant;

fn main() {
    let mut rngrng = thread_rng();
    for i in 1..127 {
        for _ in 0..1 {
            let adv = rngrng.gen_range(2u128.pow(i - 1), 2u128.pow(i));
            println!("{:} 2**{} - 2**{}", adv, i - 1, i);
            println!("{:}", adv);
            let mut rng = Xoroshiro128Plus::new(0x1234567887654321);
            let mut rng_test = Xoroshiro128Plus::new(0x1234567887654321);
            rng_test.jump(adv);
            let start = Instant::now();

            println!("{}", rng.distance(rng_test));

            let duration = start.elapsed();
            println!("Task took {:?}", duration);
        }
    }
}
