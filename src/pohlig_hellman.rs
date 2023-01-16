use crate::gf2vec::{compute_jump_poly, GF2Vec128};
use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Zero};
use rayon::{current_num_threads, prelude::*};
use std::collections::HashSet;

fn modular_inverse(a: &BigInt, m: &BigInt) -> BigInt {
    let mut a = a.to_owned();
    let mut m0 = m.to_owned();
    let mut y = BigInt::zero();
    let mut x = BigInt::one();
    if m0.is_one() {
        return BigInt::zero();
    }
    while a > BigInt::one() {
        let t = m0.to_owned() / a.to_owned();
        m0 = m0 - t.to_owned() * a.to_owned();
        y = y - t.to_owned() * x.to_owned();
        std::mem::swap(&mut a, &mut m0);
        std::mem::swap(&mut x, &mut y);
    }
    x
}

fn chinese_remainder_theorem(mods: Vec<u128>, rems: Vec<u128>) -> u128 {
    let mut total = BigInt::zero();
    let product: u128 = mods.iter().product();
    let product_biguint = product.to_bigint().unwrap();
    for i in 0..mods.len() {
        let modulus = mods[i].to_bigint().unwrap();
        let remainder = rems[i].to_bigint().unwrap();
        let p = &product_biguint / &modulus;
        let mut inv = modular_inverse(&p, &modulus);
        if inv < BigInt::zero() {
            inv += modulus.to_owned();
        }
        total += inv.to_owned() * remainder.to_owned() * p.to_owned();
    }
    let mut val = total;
    if val < BigInt::zero() {
        val -= BigInt::one()
    }
    (val % product_biguint).try_into().unwrap()
}

fn baby_step_giant_step(
    _gamma_poly: GF2Vec128,
    _h_poly: GF2Vec128,
    backwards_poly: GF2Vec128,
    char_poly: GF2Vec128,
    order: u128,
) -> u128 {
    let mut gamma_poly = GF2Vec128 {
        state_high: 0,
        state_low: 1,
    };
    let mut h_poly = _h_poly.clone();
    let step_size = (order as f64).sqrt().ceil() as u128;
    let backward_jump_poly = backwards_poly.modpow(step_size, char_poly);
    let result: Vec<u128> = (0..step_size).into_par_iter().collect();
    if step_size > current_num_threads().try_into().unwrap() {
        let lookup_vec: Vec<_> = result
            .par_chunks((step_size as usize) / current_num_threads())
            .map(|chunk| {
                let base = chunk[0];
                let mut base_poly = _gamma_poly.modpow(base, char_poly);
                chunk
                    .iter()
                    .map(|_| {
                        let state_low = base_poly.state_low;
                        base_poly = base_poly.mul(_gamma_poly).modulo(char_poly);
                        state_low
                    })
                    .collect::<Vec<u128>>()
            })
            .collect();
        let lookup_vec: Vec<_> = lookup_vec.into_iter().flat_map(|v| v.into_iter()).collect();
        let mut lookup_table = HashSet::<u128>::new();
        for item in &lookup_vec {
            lookup_table.insert(*item);
        }
        return result
            .par_chunks((step_size as usize) / current_num_threads())
            .find_map_any(|chunk| {
                let base = chunk[0];
                let mut base_poly = h_poly
                    .mul(backward_jump_poly.modpow(base, char_poly))
                    .modulo(char_poly);
                chunk.iter().find_map(|i| {
                    if lookup_table.contains(&base_poly.state_low) {
                        let j = lookup_vec
                            .iter()
                            .position(|&r| r == base_poly.state_low)
                            .unwrap() as u128;
                        return Some(i * step_size + j);
                    }
                    base_poly = base_poly.mul(backward_jump_poly).modulo(char_poly);
                    None
                })
            })
            .unwrap();
    } else {
        let mut lookup_table: Vec<u128> = vec![];
        for _ in 0..step_size {
            lookup_table.push(gamma_poly.state_low);
            gamma_poly = gamma_poly.mul(_gamma_poly).modulo(char_poly);
        }
        for i in 0..step_size {
            if lookup_table.contains(&h_poly.state_low) {
                let j = lookup_table
                    .iter()
                    .position(|&r| r == h_poly.state_low)
                    .unwrap() as u128;
                return i * step_size + j;
            }
            h_poly = h_poly.mul(backward_jump_poly).modulo(char_poly);
        }
    }

    panic!("Remainder Not Found");
}

pub fn pohlig_hellman(
    advance_poly: GF2Vec128,
    backwards_poly: GF2Vec128,
    jump_poly: GF2Vec128,
    char_poly: GF2Vec128,
    order: u128,
    primes: Vec<u128>,
) -> u128 {
    let mut remainders = vec![];
    let mut mods = vec![];
    for prime in &primes {
        let exp = order / prime;
        let g_i = advance_poly.modpow(exp, char_poly);
        let h_i = jump_poly.modpow(exp, char_poly);
        let b_i = backwards_poly.modpow(exp, char_poly);
        remainders.push(baby_step_giant_step(g_i, h_i, b_i, char_poly, *prime));
        mods.push(prime.to_owned());
        let jmp = chinese_remainder_theorem(mods.to_owned(), remainders.to_owned());
        let test_jump = compute_jump_poly(jmp, char_poly);
        if test_jump.state_low == jump_poly.state_low {
            return jmp;
        }
    }

    chinese_remainder_theorem(primes, remainders)
}
