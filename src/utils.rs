pub use std::collections::HashMap;
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
pub use std::f64::{INFINITY, NEG_INFINITY};

pub type Char = u8;
pub type Seq = Vec<Char>;
pub type SeqSlice<'a> = &'a[Char];
pub type Prob = f64;
pub type LogProb = Prob;
pub type PartitionFunc = Prob;
pub type LogPf = PartitionFunc;
pub type ExpPartOfTerm4LogPf = PartitionFunc;
pub type EpsOfTerms4LogPf = Vec<ExpPartOfTerm4LogPf>;
pub type SliceOfEpsOfTerms4LogPf<'a> = &'a[ExpPartOfTerm4LogPf];
pub type ExpPartOfTerm4LogProb = Prob;
pub type EpsOfTerms4LogProb = Vec<ExpPartOfTerm4LogProb>;
pub type SliceOfEpsOfTerms4LogProb<'a> = &'a[ExpPartOfTerm4LogProb];
pub type Hasher = BuildHasherDefault<FnvHasher>;
pub type CharPair = (Char, Char);
pub type Probs = Vec<Prob>;
pub type ProbMat = Vec<Probs>;
pub type LogProbs = Vec<LogProb>;
pub type LogProbMat = Vec<LogProbs>;
pub type LogPartialPfs = Vec<LogPf>;
pub type LogPpfMat = Vec<LogPartialPfs>;
pub type SsPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type SaScore = LogProb;
pub type CaScoreMat = HashMap<CharPair, SaScore, Hasher>;
pub struct SaScoringParams {
  pub ca_sm: CaScoreMat,
  pub base_opening_gap_penalty: SaScore,
  pub base_extending_gap_penalty: SaScore,
}
pub type ProbSeqPair = (Probs, Probs);
pub type Pos = usize;
pub type PosPair = (Pos, Pos);

#[inline]
pub fn logsumexp(xs: SliceOfEpsOfTerms4LogPf, max: ExpPartOfTerm4LogPf) -> LogPf {
  if !max.is_finite() {
    xs.iter().fold(0., |acc, &x| acc + x.exp()).ln()
  } else {
    xs.iter().fold(0., |acc, &x| acc + (x - max).exp()).ln() + max
  }
}
