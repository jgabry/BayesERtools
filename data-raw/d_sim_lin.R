# Generate data for normal linear regression

set.seed(1234)
auc <- 0:100
n <- length(auc)
d_sim_lin <-
  dplyr::tibble(
    ID = seq_along(auc),
    AUCss = auc,
    Cmaxss = auc * rnorm(n, 1, 0.3) * 5,
    BAGE = rnorm(n, 50, 10),
    SEX_num = rbinom(n, 1, 0.5),
  ) |>
  dplyr::mutate(
    Cmaxss = dplyr::if_else(Cmaxss < 0, 0, Cmaxss),
    SEX = dplyr::if_else(SEX_num == 0, "M", "F"),
    mu_resp = 0.5 * auc + 0.5 * BAGE + 5 * SEX_num,
    response = rnorm(n, mu_resp, 10)
  ) |>
  dplyr::select(-mu_resp, -SEX_num)

d_sim_lin |>
  dplyr::filter(ID %% 10 == 1) |>
  readr::write_csv("data-raw/d_sim_lin.csv")
usethis::use_data(d_sim_lin, overwrite = TRUE)
