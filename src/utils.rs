// use std::f64::consts::LOG2_E;
pub use std::collections::HashMap;
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
pub use std::f64::{INFINITY, NEG_INFINITY};

pub type Char = u8;
pub type Seq = Vec<Char>;
pub type SeqSlice<'a> = &'a[Char];
pub type Prob = f64;
pub type LogProb = Prob;
pub type PartFunc = Prob;
/* pub type LogPf = PartitionFunc;
pub type ExpPartOfTerm4LogPf = PartitionFunc;
pub type EpsOfTerms4LogPf = Vec<ExpPartOfTerm4LogPf>;
pub type SliceOfEpsOfTerms4LogPf<'a> = &'a[ExpPartOfTerm4LogPf];
pub type ExpPartOfTerm4LogProb = Prob;
pub type EpsOfTerms4LogProb = Vec<ExpPartOfTerm4LogProb>;
pub type SliceOfEpsOfTerms4LogProb<'a> = &'a[ExpPartOfTerm4LogProb]; */
pub type Hasher = BuildHasherDefault<FnvHasher>;
pub type CharPair = (Char, Char);
pub type Probs = Vec<Prob>;
pub type ProbMat = Vec<Probs>;
/* pub type LogProbs = Vec<LogProb>;
pub type LogProbMat = Vec<LogProbs>; */
/* pub type LogPartialPfs = Vec<LogPf>;
pub type LogPpfMat = Vec<LogPartialPfs>; */
pub type PartFuncs = Vec<PartFunc>;
pub type PartFuncMat = Vec<PartFuncs>;
pub type ScaleParam = PartFunc;
pub type ScaleParams = Vec<ScaleParam>;
pub type ScaleParamMat = Vec<ScaleParams>;
pub type SsPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type SaScore = LogProb;
pub type CaScoreMat = HashMap<CharPair, SaScore, Hasher>;
#[derive(Clone)]
pub struct SaScoringParams {
  pub ca_sm: CaScoreMat,
  pub opening_gap_penalty: SaScore,
  pub extending_gap_penalty: SaScore,
}
pub type ProbSeqPair = (Probs, Probs);
pub type Pos = usize;
pub type PosPair = (Pos, Pos);

// const INVERSE_LOG2_E: LogPf = 1. / LOG2_E;

/* #[inline]
pub fn logsumexp(xs: SliceOfEpsOfTerms4LogPf, max: ExpPartOfTerm4LogPf) -> LogPf {
  if !max.is_finite() {
    // fast_ln(xs.iter().fold(0., |acc, &x| acc + x.exp()))
    xs.iter().fold(0., |acc, &x| acc + x.exp()).ln()
  } else {
    xs.iter().fold(0., |acc, &x| acc + (x - max).exp()).ln() + max
  }
} */

/* #[inline]
pub fn fast_ln(x: PartitionFunc) -> LogPf {
  x.log2() * INVERSE_LOG2_E
} */
