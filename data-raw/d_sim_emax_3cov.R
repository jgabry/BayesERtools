
# core functions ----------------------------------------------------------

emax_fn <- function(exposure, emax, ec50, e0, gamma = 1) {
  e0 + emax * (exposure ^ gamma) / (ec50 ^ gamma + exposure ^ gamma)
}

generate_exposure <- function(dose, n, meanlog = 4, sdlog = 0.5) {
  dose * stats::qlnorm(
    p = stats::runif(n, min = .01, max = .99),
    meanlog = meanlog,
    sdlog = sdlog
  )
}

generate_covariate <- function(n) {
  stats::rbeta(n, 2, 2) * 10
}


# simulation functions ----------------------------------------------------


make_continuous_data <- function(seed = 123) {

  set.seed(seed)

  par <- list(
    emax   = 10,
    ec50   = 4000,
    e0     = 5,
    gamma  = 1,
    sigma  = .6,
    coef_a = .3,
    coef_b = .2,
    coef_c = 0
  )

  make_data <- function(dose, n, par) {
    tibble::tibble(
      dose = dose,
      exposure = generate_exposure(max(dose, .01), n = n),
      cov_a = generate_covariate(n = n),
      cov_b = generate_covariate(n = n),
      cov_c = generate_covariate(n = n),
      response = emax_fn(
        exposure,
        emax = par$emax,
        ec50 = par$ec50,
        e0 = par$e0,
        gamma = par$gamma
      ) +
        par$coef_a * cov_a +
        par$coef_b * cov_b +
        par$coef_c * cov_c +
        stats::rnorm(n, 0, par$sigma)
    )
  }

  dat <- dplyr::bind_rows(
    make_data(dose = 100, n = 100, par = par),
    make_data(dose = 200, n = 100, par = par),
    make_data(dose = 300, n = 100, par = par)
  )

  return(dat)
}


make_binary_data <- function(seed = 123) {

  set.seed(seed)

  par <- list(
    emax   = 5,
    ec50   = 8000,
    e0     = -3,
    gamma  = 1,
    coef_a = .2,
    coef_b = 0,
    coef_c = 0
  )

  make_data <- function(dose, n, par) {
    tibble::tibble(
      dose = dose,
      exposure = generate_exposure(max(dose, .01), n = n),
      cov_a = generate_covariate(n = n),
      cov_b = generate_covariate(n = n),
      cov_c = generate_covariate(n = n),
      pred = emax_fn( # non-linear predictor
        exposure,
        emax = par$emax,
        ec50 = par$ec50,
        e0 = par$e0,
        gamma = par$gamma
      ) +
        par$coef_a * cov_a +
        par$coef_b * cov_b +
        par$coef_c * cov_c,
      prob = 1 / (1 + exp(-pred)), # response probability
      response = as.numeric(stats::runif(n) < prob) # binary response
    ) |>
      dplyr::select(-pred, -prob)
  }

  dat <- dplyr::bind_rows(
    make_data(dose = 100, n = 100, par = par),
    make_data(dose = 200, n = 100, par = par),
    make_data(dose = 300, n = 100, par = par)
  )

  return(dat)
}


# generate data sets ------------------------------------------------------

d_sim_emax_3cov <- make_continuous_data()
d_sim_emax_bin_3cov <- make_binary_data()

readr::write_csv(d_sim_emax_3cov, "data-raw/d_sim_emax_3cov.csv")
readr::write_csv(d_sim_emax_bin_3cov, "data-raw/d_sim_emax_bin_3cov.csv")

usethis::use_data(d_sim_emax_3cov, overwrite = TRUE)
usethis::use_data(d_sim_emax_bin_3cov, overwrite = TRUE)
