extern crate bio_seq_algos;
#[macro_use]
extern crate lazy_static;

use bio_seq_algos::utils::*;
use bio_seq_algos::durbin_algo::*;

type SeqPair = (Seq, Seq);

lazy_static! {
  static ref TEST_SEQ_PAIR: SeqPair = {
    (
      String::from("AUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCCAUGCAAGGGGGCUUUAACACAUGGGAUCC").into_bytes(),
      String::from("CAAGGCGGCUUAAACACUUGGGAUGCAUGCAAGGGCGCUUUGACACAAGGUAUCCAAGCAAGCGGGCUUAAACUCAUGC").into_bytes(),
    )
  };
}

#[test]
fn test_cap_mat_and_ucp_seq_pair() {
  let mut ca_sm = CaScoreMat::default();
  let alphabet = b"AUGC";
  for &base_1 in alphabet {
    for &base_2 in alphabet {
      ca_sm.insert((base_1, base_2), if base_1 == base_2 {0.1} else {-0.1});
    }
  }
  let sa_sps = SaScoringParams::new(&ca_sm, -1., -0.1);
  let (cap_mat, ucp_seq_pair) = get_cap_mat_and_unaligned_char_psp(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_sps);
  println!("The char. alignment mat for the seq. pair \"{}\" and \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &cap_mat);
  for caps in &cap_mat {
    for &cap in caps {assert!((0. <= cap && cap <= 1.));}
  }
  println!("The unaligned char. probs. for the 1st seq. \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), &ucp_seq_pair.0);
  for ucp in ucp_seq_pair.0 {assert!((0. <= ucp && ucp <= 1.));}
  println!("The unaligned char. probs. for the 2nd seq. \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &ucp_seq_pair.1);
  for ucp in ucp_seq_pair.1 {assert!((0. <= ucp && ucp <= 1.));}
}
