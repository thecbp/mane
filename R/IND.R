#' Title
#'
#' @param n_subj
#' @param n_trts
#' @param n_periods
#' @param n_obvs
#' @param betas
#' @param y_sigma
#' @param chains
#' @param warmup
#' @param iter
#' @param adapt_delta
#' @param max_treedepth
#'
#' @return
#' @export
#'
#' @examples
IND = function(n_subj, n_trts, n_periods, n_obvs, betas, y_sigma,
               chains, warmup, iter, adapt_delta = 0.99, max_treedepth = 15) {

  # Individualized design for multi-arm adaptive N-of-1
  # Same arguments as FRN()

  # Initial FRN phase data for burn-in
  current_data = generate_FRN_data(n_subj, n_trts, n_cycles = 1, n_obvs, betas, y_sigma)

  # Storing data in tidy format, keeping for storing incoming data
  trial_data = tidy_data(current_data, n_obvs = n_obvs)

  # Generating posterior samples after initial cycle
  stan_model = mcmc_ind(current_data)

  # Index to start the loop, since n_trts have been used in FRN phase
  adaptive_start = n_trts + 1 # FRN requires one period from every treatment

  # Allocation probabilities matris
  # Indexes: [period, subject, treatment probabilities]
  allocprobs = array(0, dim = c(n_periods - n_trts, n_subj, n_trts))
  # For each treatment cycle...
  for (p in adaptive_start:n_periods) {

    # Data structures for a single individual's period
    nn = n_subj * n_obvs                        # total obs. in period
    y = id = period = array(0, nn)              # data structures
    X = matrix(0, nrow = nn, ncol = n_trts)     # data treatment matrix
    global_iter = 1                                    # global counter

    post_samp = rstan::extract(stan_model)
    S = dim(post_samp$Betas)[1]           # number of posterior samples

    # ... for each subject, calculate posterior allocation probabilities,
    # assign treatment accordingly, and then generate data
    for (i in 1:n_subj) {

      # Thompson Sampling to update allocation probability from posterior betas
      I_trt = matrix(0, S, n_trts)      # storing winning arms among all S samples
      T_mat = create_trt_matrix(n_trts)

      # Identifying the highest reward among all of the posterior samples
      for (s in 1:S) {

        # NOTE: Possible optimization here: matrix multiplication over all samples

        # Calculate expected reward in each arm for each subject
        u = post_samp$Betas[s, i ,] %*% T_mat
        best_arm = which(u == max(u))
        I_trt[s, best_arm] = 1

      }

      # Calculate allocation probabilities from above posterior rewards
      probs = apply(I_trt, 2, mean) # average the wins for each treatment arm
      c = (0.5) * (p / n_periods) # tuning parameter
      stable_probs = stabilize_probs(probs, c)
      allocprobs[(p - n_trts), i,] = stable_probs


      # Create treatment vector based on regime (1 where active, 0 else)
      next_trt = sample(1:n_trts, size = 1, prob = stable_probs)
      x = array(0, n_trts)
      x[1] = 1              # intercept
      x[next_trt] = 1       # setting next treatment

      # Generate observations for subject based on selected treatment for cycle
      for (z in 1:n_obvs) {

        X[global_iter,] = x                                  # treatment vector
        y[global_iter] = x %*% betas[i,] + rnorm(1, 0, y_sigma)  # outcome
        id[global_iter] = i                                  # subject ID
        period[global_iter] = p                              # period

        global_iter = global_iter + 1
      }

    } # end of subject loop

    # Append data to ongoing trial information
    period_data = list(J = n_subj, K = n_trts, N = nn, X = X, y = y, id = id, period = period)
    period_data_tidy = tidy_data(period_data, n_obvs = n_obvs)
    trial_data = bind_rows(trial_data, period_data_tidy)

    current_data = list(J = n_subj, K = n_trts,
                        N = nrow(trial_data),
                        X = trial_data %>% select(contains("X")) %>% as.matrix,
                        y = trial_data$y,
                        id = trial_data$id,
                        period = trial_data$period)

    # Generating posterior samples after data collected for everyone in period
    # Model here gives everyone their own prior, no aggregation
    stan_model = mcmc_ind(current_data, chains = chains, warmup = warmup, iter = iter,
                          adapt_delta = adapt_delta, max_treedepth = max_treedepth)

  } # end of period loop

  # Final aggregated data on trial to understand population-level effects
  agg_stan_model = mcmc_agg(current_data, chains = chains, warmup = warmup, iter = iter,
                            adapt_delta = adapt_delta, max_treedepth = max_treedepth)

  # Format allocation probability matrix for long format
  prob_df = tibble()
  for (i in 1:dim(allocprobs)[1]) {
    period_probs = allocprobs[i,,] %>% as_tibble()
    colnames(period_probs) = paste0("p", 1:n_trts)
    period_probs$period = i
    period_probs$id = 1:n_subj
    prob_df = bind_rows(prob_df, period_probs)
  }

  list(
    n_subj = n_subj,
    n_periods = n_periods,
    n_obvs = n_obvs,
    betas = betas,
    y_sigma = y_sigma,
    data = trial_data,
    stan_model = stan_model,
    agg_stan_model = agg_stan_model,
    allocprobs = prob_df
  )

}
