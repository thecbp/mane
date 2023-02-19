library(devtools)
load_all()
options(mc.cores = parallel::detectCores())
library(tidyverse)

# How does effect size affect the allocation probability to best treatment
betas = rbind(
  c(130, 0, 0),
  c(130, -2, 0),
  c(130, -4, 0),
  c(130, -6, 0),
  c(130, -8, 0),
  c(130, -10, 0),
  c(130, -12, 0)
) %>% as.matrix

chains = 4
iter = 2500 * 2

for (i in 94) {

  start = lubridate::now()

  sim = IND(n_subj = nrow(betas),
            n_trts = ncol(betas),
            n_periods = 50,
            n_obvs = 1,
            objective = "min",
            betas = betas,
            lag = 1,
            y_sigma = 10,
            chains = chains,
            iter = iter,
            seed = i)

  end = lubridate::now()
  dur = end - start
  print(paste0("Run took: ", dur))

  file = paste0("ES-v2-sim-run-", i, ".rds")
  write_rds(sim, file)

}
q()
