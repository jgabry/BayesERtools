if (.if_run_ex_eval_mod()) {
  d_sim_binom_cov <- d_sim_binom_cov

  id_to_sample <- seq(1, max(d_sim_binom_cov$ID), by = 5)

  df_er_ae_covsel_test <-
    d_sim_binom_cov |>
    dplyr::filter(
      AETYPE == "ae_covsel_test",
      ID %in% id_to_sample
    ) |>
    dplyr::mutate(AUCss_1000 = AUCss / 1000, BHBA1C_5 = BHBA1C / 5)

  data_er_cont <- rstanemax::exposure.response.sample
  data_er_cont_cov <- rstanemax::exposure.response.sample.with.cov
  data_er_bin <- rstanemax::exposure.response.sample.binary

  # Create  ermod_bin objects for testing ---------------------------------------
  set.seed(1234)
  ermod_bin <- dev_ermod_bin(
    data = df_er_ae_covsel_test,
    var_resp = "AEFLAG",
    var_exposure = "AUCss_1000",
    var_cov = "BHBA1C_5",
    verbosity_level = 0,
    # Below option to make the test fast
    chains = 2, iter = 1000
  )

  ermod_bin_wo_cov <- dev_ermod_bin(
    data = df_er_ae_covsel_test,
    var_resp = "AEFLAG",
    var_exposure = "AUCss_1000",
    verbosity_level = 0,
    # Below option to make the test fast
    chains = 2, iter = 1000
  )

  set.seed(1234)
  ermod_emax_w_cov <-
    dev_ermod_emax(
      data = data_er_cont_cov,
      var_exposure = "conc",
      var_resp = "resp",
      l_var_cov = list(emax = "cov2", ec50 = "cov3", e0 = "cov1"),
      verbosity_level = 0,
      chains = 2,
      iter = 1000,
      seed = 1
    )

  ermod_bin_emax <-
    suppressWarnings(dev_ermod_bin_emax(
      data = df_er_ae_covsel_test,
      var_resp = "AEFLAG",
      var_exposure = "AUCss_1000",
      verbosity_level = 0,
      # Increase iter for better convergence as exact reproducibility
      # over different machines doesn't seem realistic
      chains = 2,
      iter = 1000,
      seed = 1
    ))

  set.seed(1234)
  data(d_sim_binom_cov_hgly2)

  ermod_bin <- dev_ermod_bin(
    data = d_sim_binom_cov_hgly2 |> head(100),
    var_resp = "AEFLAG",
    var_exposure = "AUCss_1000",
    var_cov = "BHBA1C_5",
    chains = 2,
    iter = 1000
  )

  set.seed(1234)
  data_er_cont_cov <- rstanemax::exposure.response.sample.with.cov

  ermod_emax <-
    suppressWarnings(dev_ermod_emax(
      data = data_er_cont_cov,
      var_exposure = "conc",
      var_resp = "resp",
      l_var_cov = list(emax = "cov2", ec50 = "cov3", e0 = "cov1"),
      verbosity_level = 0,
      chains = 2,
      iter = 200,
      seed = 1
    ))

  # Test ----------------------------------------------------------------------
  test_that("loo", {
    loo_ermod_bin <- loo(ermod_bin)
    loo_ermod_emax_w_cov <- suppressWarnings(loo(ermod_emax_w_cov))
    loo_ermod_bin_emax <- loo(ermod_bin_emax)

    kfold_ermod_bin_rstanarm <-
      suppressMessages(rstanarm::kfold(extract_mod(ermod_bin)))
    kfold_ermod_bin <- run_kfold_cv(ermod_bin, k = 10)
    class(kfold_ermod_bin) <- c("kfold", "loo")

    kfold_ermod_bin_rstanarm
    kfold_ermod_bin
    loo::loo_compare(kfold_ermod_bin_rstanarm, kfold_ermod_bin)

    expect_equal(
      loo_ermod_bin$estimates[, 1],
      c(elpd_loo = -38.528662, p_loo = 3.325979, looic = 77.057323)
    )
    expect_equal(
      loo_ermod_emax_w_cov$estimates[, 1],
      c(elpd_loo = -216.539552, p_loo = 7.765598, looic = 433.079103),
      tolerance = 0.1
    )
    expect_equal(
      loo_ermod_bin_emax$estimates[, 1],
      c(elpd_loo = -60.787274, p_loo = 1.427765, looic = 121.574548),
      tolerance = 0.1
    )

    expect_silent(loo::loo_compare(loo_ermod_bin, loo_ermod_bin_emax))
  })

  if (requireNamespace("rsample")) {
    cv_results <- run_kfold_cv(ermod_bin, k = 3, seed = 123)
    cv_results_emax <- suppressWarnings(
      run_kfold_cv(ermod_emax, k = 3, seed = 123))

    # Test for other models are covered in test-eval_ermod.R with kfold-cv eval
    test_that("binary model", {
      expect_equal(
        names(cv_results$d_truth),
        c(".row", "truth", "fold_id")
      )
      expect_equal(
        names(cv_results$d_sim),
        c(".row", ".draw", "pred", "fold_id")
      )
      expect_equal(dim(cv_results$d_sim), c(100 * 1000, 4))

      out <- cli::cli_fmt({
        print(cv_results)
      })
      expect_true(any(grepl("k-fold Cross-Validation for ermod object", out)))
    })
  }
}
