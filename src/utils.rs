use std::f64::consts::LOG2_E;
pub use std::collections::HashMap;
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;

pub type Char = u8;
pub type Seq = Vec<Char>;
pub type SeqSlice<'a> = &'a[Char];
pub type Prob = f64;
pub type LogProb = Prob;
pub type PartitionFunc = Prob;
pub type LogPf = PartitionFunc;
pub type Pos = usize;
pub type PosPair = (Pos, Pos);
pub type Energy = Prob;
pub type ExpPartOfTerm4LogPf = PartitionFunc;
pub type EpsOfTerms4LogPf = Vec<ExpPartOfTerm4LogPf>;
pub type SliceOfEpsOfTerms4LogPf<'a> = &'a[ExpPartOfTerm4LogPf];
pub type ExpPartOfTerm4LogProb = Prob;
pub type EpsOfTerms4LogProb = Vec<ExpPartOfTerm4LogProb>;
pub type SliceOfEpsOfTerms4LogProb<'a> = &'a[ExpPartOfTerm4LogProb];
pub type Hasher = BuildHasherDefault<FnvHasher>;
pub type CharPair = (Char, Char);

const INVERSE_LOG2_E: LogPf = 1. / LOG2_E;

#[inline]
pub fn logsumexp(xs: SliceOfEpsOfTerms4LogPf, max: ExpPartOfTerm4LogPf) -> LogPf {
  if !max.is_finite() {
    fast_ln(xs.iter().fold(0., |acc, &x| acc + x.exp()))
  } else {
    fast_ln(xs.iter().fold(0., |acc, &x| acc + (x - max).exp())) + max
  }
}

#[inline]
fn fast_ln(x: PartitionFunc) -> LogPf {
  x.log2() * INVERSE_LOG2_E
}
