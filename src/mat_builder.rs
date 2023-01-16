use crate::gf2int::GF2Int;
use ndarray::{Array2, Axis, Slice};
use num_traits::One;
use std::ops::{BitXor, BitXorAssign, Shl};

pub fn mat_shl(n: usize) -> Array2<GF2Int> {
    let mut eye = Array2::<GF2Int>::zeros((64, 64));
    for i in n..64 {
        eye[[i - n, i]] = GF2Int::one();
    }
    eye
}

pub fn mat_rotl(n: usize) -> Array2<GF2Int> {
    let mut eye = Array2::<GF2Int>::zeros((64, 64));
    for i in 0..64 {
        eye[[i, (i + n) % 64]] = GF2Int::one();
    }
    eye
}

pub trait MatInverse {
    fn inverse(self) -> Self;
}

impl MatInverse for Array2<GF2Int> {
    fn inverse(self) -> Self {
        // TODO: refactor ``for k in 0..width``
        let mut mat = self.to_owned();
        let (height, width) = mat.dim();
        let mut res = Array2::<GF2Int>::eye(height);
        let mut pivot = 0;
        for i in 0..width {
            let mut is_found = false;
            for j in i..height {
                if mat[[j, i]].val == 1 {
                    if is_found {
                        for k in 0..width {
                            mat[[j, k]] = mat[[j, k]] + mat[[pivot, k]];
                            res[[j, k]] = res[[j, k]] + res[[pivot, k]];
                        }
                    } else {
                        is_found = true;
                        let temp_mat = mat
                            .slice_axis(
                                Axis(0),
                                Slice::new(
                                    j.try_into().unwrap(),
                                    Some((j + 1).try_into().unwrap()),
                                    1,
                                ),
                            )
                            .to_owned();
                        let temp_res = res
                            .slice_axis(
                                Axis(0),
                                Slice::new(
                                    j.try_into().unwrap(),
                                    Some((j + 1).try_into().unwrap()),
                                    1,
                                ),
                            )
                            .to_owned();
                        for k in 0..width {
                            mat[[j, k]] = mat[[pivot, k]];
                            mat[[pivot, k]] = temp_mat[[0, k]];
                            res[[j, k]] = res[[pivot, k]];
                            res[[pivot, k]] = temp_res[[0, k]];
                        }
                    }
                }
            }
            if is_found {
                pivot += 1;
            }
        }

        for i in (1..width).rev() {
            for j in (0..i).rev() {
                if mat[[j, i]].val == 1 {
                    for k in 0..width {
                        mat[[j, k]] = mat[[j, k]] + mat[[i, k]];
                        res[[j, k]] = res[[j, k]] + res[[i, k]];
                    }
                }
            }
        }
        res
    }
}

pub struct MatBuilder {
    pub inner_matrix: Array2<GF2Int>,
}

impl MatBuilder {
    pub fn new(
        portion_position: usize,
        portion_size: usize,
        total_state_size: usize,
    ) -> MatBuilder {
        let mut inner_matrix = Array2::<GF2Int>::zeros((total_state_size, portion_size));
        for i in 0..portion_size {
            inner_matrix[[portion_position + i, i]] = GF2Int::one();
        }
        MatBuilder { inner_matrix }
    }

    pub fn rotate_left_assign(&mut self, rhs: usize) {
        self.inner_matrix = self.inner_matrix.dot(&mat_rotl(rhs));
    }

    pub fn rotate_left(&mut self, rhs: usize) -> Self {
        MatBuilder {
            inner_matrix: self.inner_matrix.dot(&mat_rotl(rhs)),
        }
    }

    pub fn to_owned(&mut self) -> Self {
        MatBuilder {
            inner_matrix: self.inner_matrix.to_owned(),
        }
    }
}

impl BitXorAssign for MatBuilder {
    fn bitxor_assign(&mut self, rhs: Self) {
        // + == ^ under GF2, ArrayBase doesnt impl element-wise xor
        self.inner_matrix = self.inner_matrix.to_owned() + rhs.inner_matrix;
    }
}
impl BitXor for MatBuilder {
    type Output = MatBuilder;
    fn bitxor(self, rhs: Self) -> Self {
        // + == ^ under GF2, ArrayBase doesnt impl element-wise xor
        MatBuilder {
            inner_matrix: self.inner_matrix.to_owned() + rhs.inner_matrix,
        }
    }
}
impl Shl<usize> for MatBuilder {
    type Output = MatBuilder;
    fn shl(self, rhs: usize) -> Self {
        MatBuilder {
            inner_matrix: self.inner_matrix.dot(&mat_shl(rhs)),
        }
    }
}
