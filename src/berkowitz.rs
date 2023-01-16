use crate::gf2int::GF2Int;
use ndarray::{arr2, s, Array2};

pub trait CharPoly {
    type NumberType;
    fn compute_charpoly_coeffs(&mut self) -> Vec<Self::NumberType>;
}

impl CharPoly for Array2<GF2Int> {
    type NumberType = GF2Int;

    fn compute_charpoly_coeffs(&mut self) -> Vec<Self::NumberType> {
        let berk_mat = berkowitz_vector(self);

        let mut coeffs = vec![];
        for i in 0..berk_mat.shape()[0] {
            coeffs.push(berk_mat[[i, 0]]);
        }

        coeffs
    }
}

fn berkowitz_toeplitz_matrix(matrix: &mut Array2<GF2Int>) -> (Array2<GF2Int>, Array2<GF2Int>) {
    // matrix = [
    //     [a, r],
    //     [c, sub_a]
    // ]
    let (a, r) = (
        matrix[[0, 0]],
        matrix.slice_mut(s![0usize, 1usize..]).to_owned(),
    );
    let (c, sub_a) = (
        matrix.slice_mut(s![1usize.., 0usize]).to_owned(),
        matrix.slice_mut(s![1usize.., 1usize..]).to_owned(),
    );
    let mut diag_vecs = vec![c];
    for i in 0..matrix.shape()[0] - 2 {
        diag_vecs.push(sub_a.dot(&diag_vecs[i]));
    }
    let mut diags = vec![GF2Int::new(1), a];
    for d in diag_vecs {
        diags.push(r.dot(&d));
    }
    let mut toeplitz = Array2::<GF2Int>::zeros((matrix.shape()[1] + 1, matrix.shape()[0]));
    for i in 0..toeplitz.shape()[0] {
        for j in 0..toeplitz.shape()[1] {
            toeplitz[[i, j]] = if j > i { GF2Int::new(0) } else { diags[i - j] }
        }
    }
    (sub_a, toeplitz)
}

fn berkowitz_vector(matrix: &mut Array2<GF2Int>) -> Array2<GF2Int> {
    if matrix.shape() == &[0, 0] {
        return arr2(&[[GF2Int::new(1)]]);
    } else if matrix.shape() == &[1, 1] {
        return arr2(&[[GF2Int::new(1)], [matrix[[0, 0]]]]);
    }
    let (mut submat, toeplitz) = berkowitz_toeplitz_matrix(matrix);
    toeplitz.dot(&berkowitz_vector(&mut submat))
}
