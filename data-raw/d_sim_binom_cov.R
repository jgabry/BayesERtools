# Generate sample data `d_sim_binom_cov` for ER analysis with covariates
# for binary endpoint

set.seed(1234)

n_subj_sim <- 500

# Generate baseline glucose and HbA1c from multivariate normal
m_gluc_hba1c <- MASS::mvrnorm(
  n = n_subj_sim, c(6, 40),
  matrix(c(1, 2.5, 2.5, 50), 2, 2)
)

d_sim_cov <-
  dplyr::tibble(
    ID = 1:n_subj_sim,
    Dose_mg = rep(c(200, 400), each = n_subj_sim / 2),
    AUCss = stats::rlnorm(n_subj_sim, 8, 0.45),
    BAGE = stats::rnorm(n_subj_sim, 70, 8),
    BWT = stats::rlnorm(n_subj_sim, 4.5, 0.20),
    BGLUC = m_gluc_hba1c[, 1],
    BHBA1C = m_gluc_hba1c[, 2],
    RACE = sample(c("White", "Black", "Asian"), n_subj_sim,
      replace = TRUE,
      prob = c(0.7, 0.1, 0.2)
    ),
    VISC = sample(c("No", "Yes"), n_subj_sim,
      replace = TRUE,
      prob = c(0.8, 0.2)
    )
  ) |>
  dplyr::mutate(
    AUCss = AUCss * Dose_mg / 400,
    Cmaxss = AUCss / 24 * 2 * exp(stats::rnorm(n_subj_sim, 0, 0.3)),
    Cminss = AUCss / 24 / 2 * exp(stats::rnorm(n_subj_sim, 0, 0.3))
  ) |>
  dplyr::relocate(c(Cmaxss, Cminss), .after = AUCss)

# Dummy variables
d_sim_cov_dummy <-
  fastDummies::dummy_cols(d_sim_cov) |>
  dplyr::select(-c(RACE, VISC))

# Check covariates distribution
if (FALSE) {
  (
    d_sim_cov |>
      tidyr::pivot_longer(
        cols = c(AUCss, Cmaxss, Cminss, BGLUC, BHBA1C, BAGE, BWT),
        names_to = "var", values_to = "value"
      ) |>
      ggplot2::ggplot(ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = factor(Dose_mg))) +
      ggplot2::facet_wrap(~var, scales = "free")
  )
}

# Simulation with model

linpred_fun <-
  function(beta, aucss, bage, bwt, bgluc, bhba1c, race_asian, race_black,
           visc_yes) {
    beta[1] + beta[2] * aucss / 1000 + beta[3] * bage / 10 +
      beta[4] * bwt / 10 + beta[5] * bwt + beta[6] * bhba1c / 5 +
      beta[7] * race_asian + beta[8] * race_black + beta[9] * visc_yes
  }
beta_hgly2 <- c(-11, 0.5, 0, 0.1, 0.6, 0.5, 0, 0, -0.5)
beta_dr2 <- c(-4.3, 0.2, 0.3, 0, 0, 0, -0.7, 0, 0)
beta_ae_covsel_test <- c(-17, 0.5, 0, 0.1, 1, 1, 0, 0, -0.5)


d_sim_linpred <-
  d_sim_cov_dummy |>
  dplyr::mutate(
    lp_hgly2 = linpred_fun(
      beta_hgly2, AUCss, BAGE, BWT, BGLUC, BHBA1C,
      RACE_Asian, RACE_Black, VISC_Yes
    ),
    lp_dr2 = linpred_fun(
      beta_dr2, AUCss, BAGE, BWT, BGLUC, BHBA1C,
      RACE_Asian, RACE_Black, VISC_Yes
    ),
    lp_ae_covsel_test = linpred_fun(
      beta_ae_covsel_test, AUCss, BAGE, BWT, BGLUC, BHBA1C,
      RACE_Asian, RACE_Black, VISC_Yes
    )
  ) |>
  dplyr::mutate(
    p_hgly2 = boot::inv.logit(lp_hgly2),
    p_dr2 = boot::inv.logit(lp_dr2),
    p_ae_covsel_test = boot::inv.logit(lp_ae_covsel_test)
  ) |>
  dplyr::relocate(dplyr::starts_with("p_"), .after = ID)

# Check ER
if (FALSE) {
  (
    d_sim_linpred |>
      tidyr::pivot_longer(dplyr::starts_with("p_"),
        names_to = "AETYPE", values_to = "p"
      ) |>
      ggplot2::ggplot(ggplot2::aes(AUCss, p)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm") +
      ggplot2::facet_wrap(~AETYPE) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
  )
}

set.seed(1234)

d_sim_binom <-
  d_sim_linpred |>
  dplyr::mutate(
    AEFLAG_hgly2 = stats::rbinom(n_subj_sim, 1, p_hgly2),
    AEFLAG_dr2 = stats::rbinom(n_subj_sim, 1, p_dr2),
    AEFLAG_ae_covsel_test = stats::rbinom(n_subj_sim, 1, p_ae_covsel_test)
  )

d_sim_binom_cov <-
  dplyr::left_join(
    d_sim_binom |> dplyr::select(
      ID, AEFLAG_hgly2, AEFLAG_dr2,
      AEFLAG_ae_covsel_test
    ),
    d_sim_cov,
    by = "ID"
  ) |>
  tidyr::pivot_longer(
    cols = c(AEFLAG_hgly2, AEFLAG_dr2, AEFLAG_ae_covsel_test),
    names_to = "AETYPE", values_to = "AEFLAG", names_prefix = "AEFLAG_"
  ) |>
  dplyr::relocate(AETYPE, AEFLAG, .after = ID)

id_to_sample <- seq(1, max(d_sim_binom_cov$ID), by = n_subj_sim / 10)

d_sim_binom_cov |>
  dplyr::filter(ID %in% id_to_sample) |>
  readr::write_csv("data-raw/d_sim_binom_cov_sample.csv")
usethis::use_data(d_sim_binom_cov, overwrite = TRUE)

d_sim_binom_cov_hgly2 <-
  d_sim_binom_cov |>
  dplyr::filter(
    AETYPE == "hgly2"
  ) |>
  dplyr::mutate(
    AUCss_1000 = AUCss / 1000, BAGE_10 = BAGE / 10,
    BWT_10 = BWT / 10, BHBA1C_5 = BHBA1C / 5
  ) |>
  dplyr::mutate(VISC = factor(VISC, levels = c("Yes", "No")))

usethis::use_data(d_sim_binom_cov_hgly2, overwrite = TRUE)
