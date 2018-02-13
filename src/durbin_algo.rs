use utils::*;

pub struct LogSaPpfMatrices {
  pub log_sa_forward_ppf_matrix_4_char_alignment: LogPpfMatrix,
  pub log_sa_forward_ppf_matrix_4_gap_1: LogPpfMatrix,
  pub log_sa_forward_ppf_matrix_4_gap_2: LogPpfMatrix,
  pub log_sa_backward_ppf_matrix_4_char_alignment: LogPpfMatrix,
  pub log_sa_backward_ppf_matrix_4_gap_1: LogPpfMatrix,
  pub log_sa_backward_ppf_matrix_4_gap_2: LogPpfMatrix,
}

impl LogSaPpfMatrices {
  fn new(slp: &(usize, usize)) -> LogSaPpfMatrices {
    let ni_matrix = vec![vec![NEG_INFINITY; slp.1]; slp.0];
    LogSaPpfMatrices {
      log_sa_forward_ppf_matrix_4_char_alignment: ni_matrix.clone(),
      log_sa_forward_ppf_matrix_4_gap_1: ni_matrix.clone(),
      log_sa_forward_ppf_matrix_4_gap_2: ni_matrix.clone(),
      log_sa_backward_ppf_matrix_4_char_alignment: ni_matrix.clone(),
      log_sa_backward_ppf_matrix_4_gap_1: ni_matrix.clone(),
      log_sa_backward_ppf_matrix_4_gap_2: ni_matrix,
    }
  }
}

impl SaScoringParams {
  pub fn new(ca_sm: &CaScoreMatrix, base_opening_gap_penalty: SaScore, base_extending_gap_penalty: SaScore) -> SaScoringParams {
    SaScoringParams {
      ca_sm: ca_sm.clone(),
      base_opening_gap_penalty: base_opening_gap_penalty,
      base_extending_gap_penalty: base_extending_gap_penalty,
    }
  }
}

#[inline]
pub fn durbin_algo(seq_pair: &SsPair, sa_sps: &SaScoringParams) -> ProbMatrix {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let log_sa_ppf_matrices = get_log_sa_ppf_matrices(seq_pair, &seq_len_pair, sa_sps);
  let log_cap_matrix = get_log_char_alignment_prob_matrix(&log_sa_ppf_matrices, &seq_len_pair);
  get_cap_matrix(&log_cap_matrix)
}

#[inline]
pub fn get_cap_matrix(log_cap_matrix: &LogProbMatrix) -> ProbMatrix {
  log_cap_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()
}

#[inline]
pub fn get_log_cap_matrix(seq_pair: &SsPair, sa_sps: &SaScoringParams) -> LogProbMatrix {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let log_sa_ppf_matrices = get_log_sa_ppf_matrices(seq_pair, &seq_len_pair, sa_sps);
  get_log_char_alignment_prob_matrix(&log_sa_ppf_matrices, &seq_len_pair)
}

#[inline]
pub fn get_cap_matrix_and_unaligned_char_psp(sp: &SsPair, sa_sps: &SaScoringParams) -> (ProbMatrix, ProbSeqPair) {
  let slp = (sp.0.len(), sp.1.len());
  let log_sa_ppf_matrices = get_log_sa_ppf_matrices(sp, &slp, sa_sps);
  let log_cap_matrix = get_log_char_alignment_prob_matrix(&log_sa_ppf_matrices, &slp);
  let mut ucp_seq_pair = (vec![NEG_INFINITY; slp.0], vec![NEG_INFINITY; slp.1]);
  for i in 0 .. slp.0 {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. slp.1 {
      let ep_of_term_4_log_prob = log_cap_matrix[i][j];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    ucp_seq_pair.0[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  for i in 0 .. slp.1 {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
    for j in 0 .. slp.0 {
      let ep_of_term_4_log_prob = log_cap_matrix[j][i];
      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
    }
    ucp_seq_pair.1[i] = 1. - logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob).exp();
  }
  (get_cap_matrix(&log_cap_matrix), ucp_seq_pair)
}

#[inline]
pub fn get_log_sa_ppf_matrices(sp: &SsPair, slp: &(usize, usize), sa_sps: &SaScoringParams) -> LogSaPpfMatrices {
  let mut log_sa_ppf_matrices = LogSaPpfMatrices::new(slp);
  log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[0][0] = sa_sps.ca_sm[&(sp.0[0], sp.1[0])];
  for i in 1 .. slp.0 {
    log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i][0] = sa_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sa_sps) + sa_sps.ca_sm[&(sp.0[i], sp.1[0])];
  }
  for i in 1 .. slp.1 {
    log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[0][i] = sa_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sa_sps) + sa_sps.ca_sm[&(sp.0[0], sp.1[i])];
  }
  for i in 1 .. slp.0 {
    for j in 1 .. slp.1 {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i - 1][j - 1];
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_1[i - 1][j - 1];
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_2[i - 1][j - 1];
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sa_sps.ca_sm[&(sp.0[i], sp.1[j])];
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i - 1][j] + sa_sps.base_opening_gap_penalty;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_1[i - 1][j] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_1[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i][j - 1] + sa_sps.base_opening_gap_penalty;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_2[i][j - 1] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_2[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
    }
  }
  log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[slp.0 - 1][slp.1 - 1] = 0.;
  log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_1[slp.0 - 1][slp.1 - 1] = 0.;
  log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_2[slp.0 - 1][slp.1 - 1] = 0.;
  for i in 0 .. slp.0 - 1 {
    log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i][slp.1 - 1] = sa_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.0 - 1), sa_sps);
  }
  for i in 0 .. slp.1 - 1 {
    log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[slp.0 - 1][i] = sa_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.1 - 1), sa_sps);
  }
  for i in (0 .. slp.0 - 1).rev() {
    for j in (0 .. slp.1 - 1).rev() {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let ca_score = sa_sps.ca_sm[&(sp.0[i + 1], sp.1[j + 1])];
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_1[i + 1][j] + sa_sps.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_2[i][j + 1] + sa_sps.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_1[i + 1][j] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_1[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i + 1][j + 1] + ca_score;
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_2[i][j + 1] + sa_sps.base_extending_gap_penalty;
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_gap_2[i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
    }
  }
  log_sa_ppf_matrices
}

#[inline]
fn get_log_char_alignment_prob_matrix(log_sa_ppf_matrices: &LogSaPpfMatrices, slp: &(usize, usize)) -> LogProbMatrix {
  let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
  let mut max_ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[slp.0 - 1][slp.1 - 1];
  eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
  let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_1[slp.0 - 1][slp.1 - 1];
  if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
  let ep_of_term_4_log_pf = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_gap_2[slp.0 - 1][slp.1 - 1];
  if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
  let log_sa_ppf = logsumexp(&eps_of_terms_4_log_pf, max_ep_of_term_4_log_pf);
  let mut log_cap_matrix = vec![vec![NEG_INFINITY; slp.1]; slp.0];
  for i in 0 .. slp.0 {
    for j in 0 .. slp.1 {
      log_cap_matrix[i][j] = log_sa_ppf_matrices.log_sa_forward_ppf_matrix_4_char_alignment[i][j] + log_sa_ppf_matrices.log_sa_backward_ppf_matrix_4_char_alignment[i][j] - log_sa_ppf;
    }
  }
  log_cap_matrix
}

#[inline]
fn get_begp(pp: &PosPair, sa_sps: &SaScoringParams) -> SaScore {
  (pp.1 + 1 - pp.0) as SaScore * sa_sps.base_extending_gap_penalty
}
