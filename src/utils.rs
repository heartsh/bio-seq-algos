use std::f64::consts::LOG2_E;
pub use std::collections::HashMap;
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
pub use std::f64::NEG_INFINITY;

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
pub type ProbMatrix = Vec<Probs>;
pub type LogProbs = Vec<LogProb>;
pub type LogProbMatrix = Vec<LogProbs>;
pub type LogPartialPfs = Vec<LogPf>;
pub type LogPpfMatrix = Vec<LogPartialPfs>;
pub type SsPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type SaScore = LogProb;
pub type CaScoreMatrix = HashMap<CharPair, SaScore, Hasher>;
pub struct SaScoringParams {
  pub ca_sm: CaScoreMatrix,
  pub base_opening_gap_penalty: SaScore,
  pub base_extending_gap_penalty: SaScore,
}
pub type ProbSeqPair = (Probs, Probs);

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
