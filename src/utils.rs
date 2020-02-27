// pub use std::collections::HashMap;
// use fnv::FnvHasher;
pub use rustc_hash::FxHashMap;
// use std::hash::BuildHasherDefault;
pub use std::f64::NEG_INFINITY;

pub type Char = u8;
pub type Seq = Vec<Char>;
pub type SeqSlice<'a> = &'a[Char];
pub type Prob = f64;
pub type LogProb = Prob;
pub type PartFunc = Prob;
// pub type Hasher = BuildHasherDefault<FnvHasher>;
pub type CharPair = (Char, Char);
pub type Probs = Vec<Prob>;
pub type ProbMat = Vec<Probs>;
pub type PartFuncs = Vec<PartFunc>;
pub type PartFuncMat = Vec<PartFuncs>;
pub type ScaleParam = PartFunc;
pub type ScaleParams = Vec<ScaleParam>;
pub type ScaleParamMat = Vec<ScaleParams>;
pub type SsPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type SaScore = LogProb;
pub type CaScoreMat = FxHashMap<CharPair, SaScore>;
#[derive(Clone)]
pub struct SaScoringParams {
  pub ca_sm: CaScoreMat,
  pub opening_gap_penalty: SaScore,
  pub extending_gap_penalty: SaScore,
}
pub type ProbSeqPair = (Probs, Probs);
pub type Pos = usize;
pub type PosPair = (Pos, Pos);
