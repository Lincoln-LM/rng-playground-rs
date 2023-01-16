use crate::{
    berkowitz::CharPoly,
    gf2int::GF2Int,
    gf2vec::{compute_jump_poly, GF2Vec128},
    mat_builder::{MatBuilder, MatInverse},
    pohlig_hellman::pohlig_hellman,
    rng::RNG,
};
use ndarray::{concatenate, Array2, Axis};
use std::ops::BitXorAssign;

#[derive(Copy, Clone, Debug)]
pub struct Xoroshiro128Plus {
    pub seed_0: u64,
    pub seed_1: u64,
}

impl Xoroshiro128Plus {
    pub fn new(seed_0: u64) -> Xoroshiro128Plus {
        Xoroshiro128Plus {
            seed_0,
            seed_1: 0x82A2B175229D6A5B,
        }
    }
}

impl BitXorAssign for Xoroshiro128Plus {
    fn bitxor_assign(&mut self, rhs: Self) {
        self.seed_0 ^= rhs.seed_0;
        self.seed_1 ^= rhs.seed_1;
    }
}

impl RNG for Xoroshiro128Plus {
    // https://xoshiro.di.unimi.it/xoroshiro128plus.c
    // uint64_t next(void) {
    //     const uint64_t s0 = s[0];
    //     uint64_t s1 = s[1];
    //     const uint64_t result = s0 + s1;
    //     s1 ^= s0;
    //     s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
    //     s[1] = rotl(s1, 37); // c
    //     return result;
    // }
    type AdvanceInt = u128;
    type MatrixInt = GF2Int;

    fn next_state(&mut self) {
        let s0 = self.seed_0;
        let mut s1 = self.seed_1;
        // let result = s0 + s1
        s1 ^= s0;
        self.seed_0 = s0.rotate_left(24) ^ s1 ^ (s1 << 16);
        self.seed_1 = s1.rotate_left(37);
    }

    fn matrix() -> Array2<Self::MatrixInt> {
        let mut s0_mat = MatBuilder::new(0, 64, 128);
        let mut s1_mat = MatBuilder::new(64, 64, 128);
        s1_mat ^= s0_mat.to_owned();
        s0_mat = s0_mat.rotate_left(24) ^ s1_mat.to_owned() ^ (s1_mat.to_owned() << 16);
        s1_mat.rotate_left_assign(37);

        concatenate![Axis(1), s0_mat.inner_matrix, s1_mat.inner_matrix]
    }

    fn advance(&mut self, adv: Self::AdvanceInt) {
        for _ in 0..adv {
            self.next_state();
        }
    }

    fn jump(&mut self, jmp: Self::AdvanceInt) {
        // can be precomputed
        let mut char_poly = Xoroshiro128Plus::matrix().compute_charpoly_coeffs();
        char_poly.reverse();
        let char_poly = GF2Vec128::new(char_poly);
        // can be precomputed via jump table
        let jump_poly = compute_jump_poly(jmp, char_poly);
        let mut final_state = Xoroshiro128Plus {
            seed_0: 0,
            seed_1: 0,
        };
        for bit in 0..128 {
            if ((jump_poly.state_low >> bit) & 1) != 0 {
                final_state ^= *self;
            }
            self.next_state();
        }
        self.seed_0 = final_state.seed_0;
        self.seed_1 = final_state.seed_1;
    }

    fn distance(&mut self, other: Self) -> Self::AdvanceInt {
        // can be precomputed
        let mut char_poly = Xoroshiro128Plus::matrix().compute_charpoly_coeffs();
        char_poly.reverse();
        let char_poly = GF2Vec128::new(char_poly);

        let mut start = self.clone();
        let mut end = other.clone();
        let mut jump_application_mat = Array2::<GF2Int>::zeros((128, 128));
        for i in 0..128 {
            let state = start.state();
            for j in 0..128 {
                jump_application_mat[[i, j]] = GF2Int::new(((state >> j) & 1) as u8);
            }
            start.next_state();
        }
        let jump_application_mat_inv = jump_application_mat.inverse();
        let mut jump_poly_mat = Array2::<GF2Int>::zeros((1, 128));
        let end_state = end.state();
        for i in 0..128 {
            jump_poly_mat[[0, i]] = GF2Int::new(((end_state >> i) & 1) as u8);
        }
        jump_poly_mat = jump_poly_mat.dot(&jump_application_mat_inv);
        let mut jump_poly = vec![];
        for i in 0..jump_poly_mat.shape()[1] {
            jump_poly.push(jump_poly_mat[[0, i]]);
        }
        let jump_poly = GF2Vec128::new(jump_poly);
        pohlig_hellman(
            GF2Vec128 {
                state_low: 0b10,
                state_high: 0,
            },
            // -1 ≡ mod 2**128 - 2 (mod 2**128 - 1)
            compute_jump_poly(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE, char_poly),
            jump_poly,
            char_poly,
            0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,
            // precomputed prime-factorization due to integer size limits
            vec![3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721],
        )
    }

    fn state(&mut self) -> Self::AdvanceInt {
        (self.seed_0 as u128) | ((self.seed_1 as u128) << 64u128)
    }
}
