use utils::*;

pub struct LogSaPpfMats {
  pub log_sa_forward_ppf_mat_4_char_align: LogPpfMat,
  pub log_sa_forward_ppf_mat_4_gap_1: LogPpfMat,
  pub log_sa_forward_ppf_mat_4_gap_2: LogPpfMat,
  pub log_sa_backward_ppf_mat_4_char_align: LogPpfMat,
  pub log_sa_backward_ppf_mat_4_gap_1: LogPpfMat,
  pub log_sa_backward_ppf_mat_4_gap_2: LogPpfMat,
}

impl LogSaPpfMats {
  fn new(slp: &(usize, usize)) -> LogSaPpfMats {
    let ni_mat = vec![vec![NEG_INFINITY; slp.1]; slp.0];
    LogSaPpfMats {
      log_sa_forward_ppf_mat_4_char_align: ni_mat.clone(),
      log_sa_forward_ppf_mat_4_gap_1: ni_mat.clone(),
      log_sa_forward_ppf_mat_4_gap_2: ni_mat.clone(),
      log_sa_backward_ppf_mat_4_char_align: ni_mat.clone(),
      log_sa_backward_ppf_mat_4_gap_1: ni_mat.clone(),
      log_sa_backward_ppf_mat_4_gap_2: ni_mat,
    }
  }
}

impl SaScoringParams {
  pub fn new(ca_sm: &CaScoreMat, base_opening_gap_penalty: SaScore, base_extending_gap_penalty: SaScore) -> SaScoringParams {
    SaScoringParams {
      ca_sm: ca_sm.clone(),
      base_opening_gap_penalty: base_opening_gap_penalty,
      base_extending_gap_penalty: base_extending_gap_penalty,
    }
  }
}

#[inline]
pub fn durbin_algo(seq_pair: &SsPair, sa_sps: &SaScoringParams) -> ProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let log_sa_ppf_mats = get_log_sa_ppf_mats(seq_pair, &seq_len_pair, sa_sps);
  let log_cap_mat = get_log_char_align_prob_mat(&log_sa_ppf_mats, &seq_len_pair);
  get_cap_mat(&log_cap_mat)
}

#[inline]
pub fn get_cap_mat(log_cap_mat: &LogProbMat) -> ProbMat {
  log_cap_mat.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()
}

#[inline]
pub fn get_log_cap_mat(seq_pair: &SsPair, sa_sps: &SaScoringParams) -> LogProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let log_sa_ppf_mats = get_log_sa_ppf_mats(seq_pair, &seq_len_pair, sa_sps);
  get_log_char_align_prob_mat(&log_sa_ppf_mats, &seq_len_pair)
}

#[inline]
pub fn get_cap_mat_and_unaligned_char_psp(sp: &SsPair, sa_sps: &SaScoringParams) -> (ProbMat, ProbSeqPair) {
  let slp = (sp.0.len(), sp.1.len());
  let log_sa_ppf_mats = get_log_sa_ppf_mats(sp, &slp, sa_sps);
  let log_cap_mat = get_log_char_align_prob_mat(&log_sa_ppf_mats, &slp);
  let mut ucp_seq_pair = (vec![NEG_INFINITY; slp.0], vec![NEG_INFINITY; slp.1]);
  for i in 0 .. slp.0 {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. slp.1 {
      let ep_of_term_4_log_prob = log_cap_mat[i][j];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    ucp_seq_pair.0[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  for i in 0 .. slp.1 {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. slp.0 {
      let ep_of_term_4_log_prob = log_cap_mat[j][i];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    ucp_seq_pair.1[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  (get_cap_mat(&log_cap_mat), ucp_seq_pair)
}

#[inline]
pub fn get_log_sa_ppf_mats(sp: &SsPair, slp: &(usize, usize), sa_sps: &SaScoringParams) -> LogSaPpfMats {
  let mut log_sa_ppf_mats = LogSaPpfMats::new(slp);
  log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[0][0] = sa_sps.ca_sm[&(sp.0[0], sp.1[0])];
  for i in 1 .. slp.0 {
    log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i][0] = sa_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sa_sps) + sa_sps.ca_sm[&(sp.0[i], sp.1[0])];
  }
  for i in 1 .. slp.1 {
    log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[0][i] = sa_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sa_sps) + sa_sps.ca_sm[&(sp.0[0], sp.1[i])];
  }
  for i in 1 .. slp.0 {
    for j in 1 .. slp.1 {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i - 1][j - 1];
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_1[i - 1][j - 1];
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_2[i - 1][j - 1];
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sa_sps.ca_sm[&(sp.0[i], sp.1[j])];
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i - 1][j] + sa_sps.base_opening_gap_penalty;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_1[i - 1][j] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_1[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i][j - 1] + sa_sps.base_opening_gap_penalty;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_2[i][j - 1] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_2[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
    }
  }
  log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[slp.0 - 1][slp.1 - 1] = 0.;
  log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_1[slp.0 - 1][slp.1 - 1] = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[slp.0 - 1][slp.1 - 1];
  log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_2[slp.0 - 1][slp.1 - 1] = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[slp.0 - 1][slp.1 - 1];
  for i in 0 .. slp.0 - 1 {
    log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i][slp.1 - 1] = sa_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.0 - 1), sa_sps);
  }
  for i in 0 .. slp.1 - 1 {
    log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[slp.0 - 1][i] = sa_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.1 - 1), sa_sps);
  }
  for i in (0 .. slp.0 - 1).rev() {
    for j in (0 .. slp.1 - 1).rev() {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let ca_score = sa_sps.ca_sm[&(sp.0[i + 1], sp.1[j + 1])];
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_1[i + 1][j] + sa_sps.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_2[i][j + 1] + sa_sps.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_1[i + 1][j] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_1[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_2[i][j + 1] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_mats.log_sa_backward_ppf_mat_4_gap_2[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
    }
  }
  log_sa_ppf_mats
}

#[inline]
fn get_log_char_align_prob_mat(log_sa_ppf_mats: &LogSaPpfMats, slp: &(usize, usize)) -> LogProbMat {
  let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
  let mut max_ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[slp.0 - 1][slp.1 - 1];
  eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
  let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_1[slp.0 - 1][slp.1 - 1];
  if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
  let ep_of_term_4_log_pf = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_gap_2[slp.0 - 1][slp.1 - 1];
  if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
  let log_sa_ppf = logsumexp(&eps_of_terms_4_log_pf, max_ep_of_term_4_log_pf);
  let mut log_cap_mat = vec![vec![NEG_INFINITY; slp.1]; slp.0];
  for i in 0 .. slp.0 {
    for j in 0 .. slp.1 {
      log_cap_mat[i][j] = log_sa_ppf_mats.log_sa_forward_ppf_mat_4_char_align[i][j] + log_sa_ppf_mats.log_sa_backward_ppf_mat_4_char_align[i][j] - log_sa_ppf;
    }
  }
  log_cap_mat
}

#[inline]
fn get_begp(pp: &PosPair, sa_sps: &SaScoringParams) -> SaScore {
  (pp.1 + 1 - pp.0) as SaScore * sa_sps.base_extending_gap_penalty
}
