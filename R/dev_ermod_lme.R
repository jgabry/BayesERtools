# Linear Mixed Effect Model Development
# `rstanarm` package required
#' @export
#' @examples
#' \donttest{
#' data(d_sim_lme)
#'
#' ermod_lme <- dev_ermod_lme(
#'   data = d_sim_lme,
#'   var_resp = "response",
#'   var_exposure = "AUCss",
#'   var_cov = c("SEX", "BAGE"),
#'   var_random = "SUBJECT",
#'   random_effect_type = "both"
#' )
#'
#' ermod_lme
#' }
#'
dev_ermod_lme <- function(
    data,
    var_resp,
    var_exposure,
    var_cov = NULL,
    var_random,
    random_effect_type = c("both", "intercept", "slope"),
    prior = rstanarm::default_prior_coef(),
    prior_intercept = rstanarm::default_prior_intercept(),
    prior_aux = rstanarm::exponential(autoscale = TRUE),
    verbosity_level = 1,
    chains = 4,
    iter = 2000) {
  stopifnot(verbosity_level %in% c(0, 1, 2, 3))
  random_effect_type <- match.arg(random_effect_type)
  refresh <- dplyr::if_else(verbosity_level >= 3, iter %/% 4, 0)

  input_args <- capture_selected_args(
    c("chains", "iter"),
    environment()
  )

  check_data_columns(
    data = data,
    var_exposure = var_exposure,
    var_resp = var_resp,
    var_cov = var_cov,
    var_random = var_random
  )

  var_full <- c(var_exposure, var_cov)

  # Build the formula for the fixed effects
  formula_fixed <- paste(var_resp, "~", paste(var_full, collapse = " + "))

  # Add the random effect based on the specified type
  if (random_effect_type == "intercept") {
    formula_full <- paste(formula_fixed, "+ (1 |", var_random, ")")
  } else if (random_effect_type == "slope") {
    formula_full <- paste(formula_fixed, "+ (0 +", var_exposure, "|", var_random, ")")
  } else if (random_effect_type == "both") {
    formula_full <- paste(formula_fixed, "+ (1 +", var_exposure, "|", var_random, ")")
  }

  # Convert to a formula object
  formula_final <- stats::as.formula(formula_full)

  # Use rlang::call2 to create the call to stan_glmer
  call_stan_glmer <- rlang::call2(
    rstanarm::stan_glmer,
    formula = formula_final,
    family = stats::gaussian(),
    data = quote(data),
    prior = prior,
    prior_intercept = prior_intercept,
    prior_aux = prior_aux,
    chains = chains,
    iter = iter,
    refresh = refresh
  )
  mod <- eval(call_stan_glmer)

  new_ermod_lme(
    mod = mod,
    data = data,
    var_resp = var_resp,
    var_exposure = var_exposure,
    var_cov = var_cov,
    var_random = var_random,
    random_effect_type = random_effect_type,
    input_args = input_args
  )
}
