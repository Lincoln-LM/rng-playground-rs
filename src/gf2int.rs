use num_traits::identities::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Copy)]
pub struct GF2Int {
    pub val: u8,
}

impl GF2Int {
    pub fn new(val: u8) -> GF2Int {
        GF2Int { val: val & 1 }
    }
}

impl fmt::Display for GF2Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}

impl fmt::Debug for GF2Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}

impl Zero for GF2Int {
    fn zero() -> GF2Int {
        GF2Int { val: 0 }
    }

    fn is_zero(&self) -> bool {
        self.val == 0
    }
}

impl One for GF2Int {
    fn one() -> GF2Int {
        GF2Int { val: 1 }
    }

    fn is_one(&self) -> bool {
        self.val == 1
    }
}

impl Add for GF2Int {
    type Output = GF2Int;

    #[inline(always)]
    fn add(self, other: GF2Int) -> GF2Int {
        GF2Int {
            val: self.val ^ other.val,
        }
    }
}

impl Sub for GF2Int {
    type Output = GF2Int;

    fn sub(self, other: GF2Int) -> GF2Int {
        GF2Int {
            val: self.val ^ other.val,
        }
    }
}

impl Mul for GF2Int {
    type Output = GF2Int;

    #[inline(always)]
    fn mul(self, other: GF2Int) -> GF2Int {
        GF2Int {
            val: self.val & other.val,
        }
    }
}

// trait needs to be implemented for LinalgScalar but should never be used
impl Div for GF2Int {
    type Output = GF2Int;

    fn div(self, _other: GF2Int) -> GF2Int {
        panic!("Division is undefined under GF2!")
    }
}
