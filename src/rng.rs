use ndarray::Array2;

pub trait RNG {
    type AdvanceInt;
    type MatrixInt;
    fn next_state(&mut self);
    fn matrix() -> Array2<Self::MatrixInt>;
    fn advance(&mut self, adv: Self::AdvanceInt);
    fn jump(&mut self, jmp: Self::AdvanceInt);
    fn distance(&mut self, other: Self) -> Self::AdvanceInt;
    fn state(&mut self) -> Self::AdvanceInt;
}
