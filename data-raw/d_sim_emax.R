
# data generation - no covariate ------------------------------------------

set.seed(1234)

n <- 20 # number of subjects
E0 <- 5 # effect at 0 concentration
Emax <- 10 # maximal effect
EC50 <- 20 # concentration at half maximal effect
h <- 2 # Hill coefficient

set.seed(130)
c.is <- 50 * stats::runif(n) # exposure

set.seed(130)
eps <- stats::rnorm(n) # residual error

y.is <- E0 + ((Emax * c.is^h) / (c.is^h + EC50^h)) + eps

d_sim_emax_nocov <- tibble::tibble(Conc = c.is, Y = y.is)

# data generation - covariate  --------------------------------------------

Ngp <- 2
N <- 20 * Ngp
GPid <- as.factor(rep(c("A", "B"), each = 20))

# set parameters
E0 <- 5
EC50 <- 15
h <- 2
beta1 <- .7

set.seed(12345)

d_sim_emax_cov <-
  tibble::tibble(GP = GPid) |>
  dplyr::mutate(c.is = 50 * runif(N), eps = rnorm(N)) |>
  dplyr::mutate(Emax.i = ifelse(GP == "A", 10, 15)) |>
  dplyr::mutate(y.is = E0 + ((Emax.i * c.is^h) / (c.is^h + EC50^h)) + eps) |>
  dplyr::mutate(Conc = c.is, Y = y.is)


# write data --------------------------------------------------------------

readr::write_csv(d_sim_emax_nocov, "data-raw/d_sim_emax_nocov.csv")
readr::write_csv(d_sim_emax_cov, "data-raw/d_sim_emax_cov.csv")

usethis::use_data(d_sim_emax_nocov, overwrite = TRUE)
usethis::use_data(d_sim_emax_cov, overwrite = TRUE)

