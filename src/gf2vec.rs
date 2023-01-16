use crate::gf2int::GF2Int;
use std::cmp::min;

#[derive(Debug, Clone, Copy)]
pub struct GF2Vec128 {
    pub state_low: u128,
    pub state_high: u128,
}

impl GF2Vec128 {
    pub fn new(vec: Vec<GF2Int>) -> GF2Vec128 {
        let mut state_low: u128 = 0;
        let mut state_high: u128 = 0;
        for i in 0..min(128, vec.len()) {
            state_low |= (vec[i].val as u128) << i;
        }
        for i in 128..min(256, vec.len()) {
            state_high |= (vec[i].val as u128) << (i - 128);
        }
        GF2Vec128 {
            state_low,
            state_high,
        }
    }

    pub fn last_bit_pos(self) -> u32 {
        // common enough case to account for
        if self.state_high == 1 {
            return 128;
        } else if self.state_high != 0 {
            let mut state = self.state_high;
            let mut result = 128;
            while state != 0 {
                state >>= 1;
                result += 1;
            }
            return result;
        } else if self.state_low != 0 {
            let mut state = self.state_low;
            let mut result = 0;
            while state != 0 {
                state >>= 1;
                result += 1;
            }
            return result;
        } else {
            // should never happen
            panic!("Getting MSSB of 0 vector");
        }
    }

    pub fn is_zero(self) -> bool {
        return (self.state_high == 0) && (self.state_low == 0);
    }

    pub fn is_one(self) -> bool {
        return (self.state_high == 0) && (self.state_low == 1);
    }

    pub fn modulo(self, rhs: GF2Vec128) -> GF2Vec128 {
        let mut polynomial = self;
        let last_bit_pos = rhs.last_bit_pos();
        if polynomial.shr(last_bit_pos).is_zero() {
            return polynomial;
        }
        let poly_mssb = polynomial.last_bit_pos();
        let shift_num = poly_mssb - last_bit_pos;
        let mut modulus = rhs.shl(shift_num);
        for shift_pos in 0..shift_num + 1 {
            if polynomial.is_zero() {
                return GF2Vec128 {
                    state_low: 0,
                    state_high: 0,
                };
            }
            if polynomial.shr(poly_mssb - shift_pos).is_one() {
                polynomial = polynomial.bitxor(modulus);
            }
            modulus = modulus.shr(1);
        }
        polynomial
    }

    pub fn mul(self, rhs: GF2Vec128) -> GF2Vec128 {
        let mut multiplicand = self;
        let mut multiplier = rhs;
        let mut result = GF2Vec128 {
            state_low: 0,
            state_high: 0,
        };

        while !multiplier.is_zero() {
            if (multiplier.state_low & 1) != 0 {
                result = result.bitxor(multiplicand);
            }
            multiplicand = multiplicand.shl(1);
            multiplier = multiplier.shr(1);
        }
        result
    }

    pub fn modpow(self, rhs: u128, modulus: GF2Vec128) -> GF2Vec128 {
        let mut power = rhs;
        let mut base = self;
        let mut result = GF2Vec128 {
            state_low: 1,
            state_high: 0,
        };
        // exponentiation by squares
        while power > 0 {
            if (power & 1) == 1 {
                result = result.mul(base).modulo(modulus);
            }
            power >>= 1;
            base = base.mul(base).modulo(modulus);
        }

        result
    }

    pub fn bitxor(self, rhs: GF2Vec128) -> GF2Vec128 {
        return GF2Vec128 {
            state_low: self.state_low ^ rhs.state_low,
            state_high: self.state_high ^ rhs.state_high,
        };
    }

    pub fn shr(self, rhs: u32) -> GF2Vec128 {
        if rhs == 0 {
            return self;
        }
        if rhs >= 128 {
            return GF2Vec128 {
                state_low: self.state_high >> (rhs - 128),
                state_high: 0,
            };
        }
        let state_low = (self.state_low >> rhs) | (self.state_high << (128 - rhs));
        let state_high = self.state_high >> rhs;
        GF2Vec128 {
            state_low,
            state_high,
        }
    }

    pub fn shl(self, rhs: u32) -> GF2Vec128 {
        if rhs == 0 {
            return self;
        }
        if rhs >= 128 {
            return GF2Vec128 {
                state_low: 0,
                state_high: self.state_low << (rhs - 128),
            };
        }
        let state_high = (self.state_high << rhs) | (self.state_low >> (128 - rhs));
        let state_low = self.state_low << rhs;
        GF2Vec128 {
            state_low,
            state_high,
        }
    }
}

pub fn base_z_modpow(power: u128, modulus: GF2Vec128) -> GF2Vec128 {
    let base = GF2Vec128 {
        state_low: 0b10, // z ** 1
        state_high: 0,
    };
    base.modpow(power, modulus)
}

pub fn compute_jump_poly(jmp: u128, char_poly: GF2Vec128) -> GF2Vec128 {
    base_z_modpow(jmp, char_poly)
}
