use utils::*;

pub struct SaPartFuncMats {
  pub sa_forward_part_func_mat: PartFuncMat,
  pub sa_forward_part_func_mat_4_char_align: PartFuncMat,
  pub sa_forward_part_func_mat_4_gap_1: PartFuncMat,
  pub sa_forward_part_func_mat_4_gap_2: PartFuncMat,
  pub sa_backward_part_func_mat_4_char_align: PartFuncMat,
  pub sa_backward_part_func_mat_4_gap_1: PartFuncMat,
  pub sa_backward_part_func_mat_4_gap_2: PartFuncMat,
}

impl SaPartFuncMats {
  fn new(slp: &(usize, usize)) -> SaPartFuncMats {
    let zero_mat = vec![vec![0.; slp.1 + 2]; slp.0 + 2];
    SaPartFuncMats {
      sa_forward_part_func_mat: zero_mat.clone(),
      sa_forward_part_func_mat_4_char_align: zero_mat.clone(),
      sa_forward_part_func_mat_4_gap_1: zero_mat.clone(),
      sa_forward_part_func_mat_4_gap_2: zero_mat.clone(),
      sa_backward_part_func_mat_4_char_align: zero_mat.clone(),
      sa_backward_part_func_mat_4_gap_1: zero_mat.clone(),
      sa_backward_part_func_mat_4_gap_2: zero_mat,
    }
  }
}

impl SaScoringParams {
  pub fn new(ca_sm: &CaScoreMat, opening_gap_penalty: SaScore, extending_gap_penalty: SaScore) -> SaScoringParams {
    SaScoringParams {
      ca_sm: ca_sm.clone(),
      opening_gap_penalty: opening_gap_penalty,
      extending_gap_penalty: extending_gap_penalty,
    }
  }
}

pub fn durbin_algo(seq_pair: &SsPair, sa_sps: &SaScoringParams) -> ProbMat {
  let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
  let (sa_part_func_mats, scale_param_mat) = get_sa_part_func_mats_and_scale_param_mat(seq_pair, &seq_len_pair, sa_sps);
  get_char_align_prob_mat(&sa_part_func_mats, &seq_len_pair, &scale_param_mat)
}

pub fn get_cap_mat_and_unaligned_char_psp(sp: &SsPair, sa_sps: &SaScoringParams) -> (ProbMat, ProbSeqPair) {
  let slp = (sp.0.len(), sp.1.len());
  let cap_mat = durbin_algo(sp, sa_sps);
  let mut ucp_seq_pair = (vec![0.; slp.0], vec![0.; slp.1]);
  for i in 0 .. slp.0 {
    let mut ucp = 1.;
    for j in 0 .. slp.1 {
      ucp -= cap_mat[i][j];
    }
    assert!(0. <= ucp && ucp <= 1.);
    ucp_seq_pair.0[i] = ucp;
  }
  for j in 0 .. slp.1 {
    let mut ucp = 1.;
    for i in 0 .. slp.0 {
      ucp -= cap_mat[i][j];
    }
    assert!(0. <= ucp && ucp <= 1.);
    ucp_seq_pair.1[j] = ucp;
  }
  (cap_mat, ucp_seq_pair)
}

pub fn get_sa_part_func_mats_and_scale_param_mat(sp: &SsPair, slp: &(usize, usize), sa_sps: &SaScoringParams) -> (SaPartFuncMats, ScaleParamMat) {
  let mut exp_sa_sps = sa_sps.clone();
  for ca_score in exp_sa_sps.ca_sm.values_mut() {
    *ca_score = ca_score.exp();
  };
  exp_sa_sps.opening_gap_penalty = exp_sa_sps.opening_gap_penalty.exp();
  exp_sa_sps.extending_gap_penalty = exp_sa_sps.extending_gap_penalty.exp();
  let mut sa_part_func_mats = SaPartFuncMats::new(slp);
  let mut scale_param_mat = vec![vec![0.; slp.1 + 2]; slp.0 + 2];
  for i in 0 .. slp.0 + 2 {
    for j in 0 .. slp.1 + 2 {
      if (i == 0 && j == 0) || (i == slp.0 + 1 && j == slp.1 + 1) {
        sa_part_func_mats.sa_forward_part_func_mat[i][j] = 1.;
        sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i][j] = 1.;
        scale_param_mat[i][j] = 1.;
        continue;
      } else if i == slp.0 + 1 || j == slp.1 + 1 {
        scale_param_mat[i][j] = 1.;
        continue;
      }
      let mut scale_param = 0.;
      if i > 0 && j > 0 {
        let ca_score = exp_sa_sps.ca_sm[&(sp.0[i - 1], sp.1[j - 1])];
        let forward_part_func = sa_part_func_mats.sa_forward_part_func_mat[i - 1][j - 1] * ca_score;
        scale_param += forward_part_func;
        sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i][j] = forward_part_func;
      }
      if i > 0 {
        let forward_part_func = (sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i - 1][j] + sa_part_func_mats.sa_forward_part_func_mat_4_gap_2[i - 1][j]) * exp_sa_sps.opening_gap_penalty + sa_part_func_mats.sa_forward_part_func_mat_4_gap_1[i - 1][j] * exp_sa_sps.extending_gap_penalty;
        scale_param += forward_part_func;
        sa_part_func_mats.sa_forward_part_func_mat_4_gap_1[i][j] = forward_part_func;
      }
      if j > 0 {
        let forward_part_func = (sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i][j - 1] + sa_part_func_mats.sa_forward_part_func_mat_4_gap_1[i][j - 1]) * exp_sa_sps.opening_gap_penalty + sa_part_func_mats.sa_forward_part_func_mat_4_gap_2[i][j - 1] * exp_sa_sps.extending_gap_penalty;
        scale_param += forward_part_func;
        sa_part_func_mats.sa_forward_part_func_mat_4_gap_2[i][j] = forward_part_func;
      }
      sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i][j] /= scale_param;
      sa_part_func_mats.sa_forward_part_func_mat_4_gap_1[i][j] /= scale_param;
      sa_part_func_mats.sa_forward_part_func_mat_4_gap_2[i][j] /= scale_param;
      sa_part_func_mats.sa_forward_part_func_mat[i][j] = 1.;
      scale_param_mat[i][j] = scale_param;
    }
  }
  for i in (0 .. slp.0 + 2).rev() {
    for j in (0 .. slp.1 + 2).rev() {
      if i == slp.0 + 1 && j == slp.1 + 1 {
        sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] = 1.;
        continue;
      } else if i == slp.0 + 1 || j == slp.1 + 1 || i == 0 || j == 0 {
        continue;
      }
      let ca_score = if i + 1 == slp.0 + 1 && j + 1 == slp.1 + 1 {1.} else if i + 1 == slp.0 + 1 || j + 1 == slp.1 + 1 {0.} else {exp_sa_sps.ca_sm[&(sp.0[i], sp.1[j])]};
      let part_func = sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i + 1][j + 1] * ca_score;
      sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] += part_func;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i][j] += part_func;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j] += part_func;
      sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i + 1][j] * exp_sa_sps.opening_gap_penalty;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i + 1][j] * exp_sa_sps.extending_gap_penalty;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i + 1][j] * exp_sa_sps.opening_gap_penalty;
      sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j + 1] * exp_sa_sps.opening_gap_penalty;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j + 1] * exp_sa_sps.opening_gap_penalty;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j] += sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j + 1] * exp_sa_sps.extending_gap_penalty;
      let scale_param = scale_param_mat[i][j];
      sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] /= scale_param;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_1[i][j] /= scale_param;
      sa_part_func_mats.sa_backward_part_func_mat_4_gap_2[i][j] /= scale_param;
    }
  }
  (sa_part_func_mats, scale_param_mat)
}

fn get_char_align_prob_mat(sa_part_func_mats: &SaPartFuncMats, slp: &(usize, usize), scale_param_mat: &ScaleParamMat) -> ProbMat {
  let mut cap_mat = vec![vec![0.; slp.1]; slp.0];
  for i in 1 .. slp.0 + 1 {
    for j in 1 .. slp.1 + 1 {
      let cap = sa_part_func_mats.sa_forward_part_func_mat_4_char_align[i][j] * sa_part_func_mats.sa_backward_part_func_mat_4_char_align[i][j] * scale_param_mat[i][j];
      assert!(0. <= cap && cap <= 1.);
      cap_mat[i - 1][j - 1] = cap;
    }
  }
  cap_mat
}
