extern crate bio_seq_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use bio_seq_algos::utils::*;
use bio_seq_algos::durbin_algo::*;
use time::precise_time_s;

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
      ca_sm.insert((base_1, base_2), if base_1 == base_2 {1.} else {-1.});
    }
  }
  let sa_sps = SaScoringParams::new(&ca_sm, -11., -1.);
  let begin = precise_time_s();
  let _ = get_cap_mat_and_unaligned_char_psp(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_sps);
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time  = {} [s].", elapsed_time);
}
