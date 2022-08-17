#' Simulate an individual bandit design for a multi-arm N-of-1 trial
#'
#' `IND()` simulates an individual bandit design for a multi-arm N-of-1 trial.
#' This design comes in two phases: a burn-in phase where all the treatments are
#' gievn in a random order, and an adaptive phase where the next treatment
#' randomized is chosen via Thompson Sampling. The allocation probabilities are
#' calculated via the posterior probability of each treatment being optimal.
#'
#' In the individual bandit design, the posterior treatment effects are
#' updated after all individual have finished a treatment period. It is assumed
#' that estimating the individual-level effects are of interest.
#'
#' @inheritParams FRN
#'
#' @return List containing the simulation parameters and resulting trial data that came from the parameters
#' @export
#'
#' @examples
#' \dontrun{
#' # Generating data for simulation
#' set.seed(1)
#'
#' # Covariance between betas, factorized into scalar and correlation components
#' tau = 2 * diag(3)
#' Omega = diag(3)
#' Omega[1,2] = Omega[2,1] = 0.1
#' Omega[1,3] = Omega[3,1] = 0.1
#' Omega[3,2] = Omega[2,3] = 0.1
#' Sigma = tau %*% Omega %*% tau
#'
#' # Generating individual treatment effects from unique distributions
#' betas = array(0, dim = c(2, 3))
#' betas[1,] = MASS::mvrnorm(n = 1, mu = c(3, -2, 4), Sigma = Sigma)
#' betas[2,] = MASS::mvrnorm(n = 1, mu = c(1, 0, 3), Sigma = Sigma)
#' betas = round(betas, 1)
#' stanfile = 'stan/model_ind.stan'
#'
#' # ind = IND(n_subj = 2, n_periods = 6, n_obvs = 5, n_trts = 3, y_sigma = 2,
#' #           stanfile = stanfile, betas = betas, chains = 1, warmup = 1000,
#' #           iter = 3000)
#' }
IND = function(n_subj, n_trts, n_periods, n_obvs, betas, y_sigma, stanfile,
               chains, warmup, iter, adapt_delta = 0.99, max_treedepth = 15) {

  # Individualized design for multi-arm adaptive N-of-1
  # Same arguments as FRN()

  # Initial FRN phase data for burn-in
  current_data = generate_FRN_data(n_subj, n_trts, n_periods = 1, n_obvs, betas, y_sigma)

  # Storing data in tidy format, keeping for storing incoming data
  trial_data = tidy_data(current_data)

  # Generating posterior samples after initial cycle
  stan_model = mcmc(stanfile, current_data, chains = chains, warmup = warmup,
                    iter = iter, adapt_delta = adapt_delta,
                    max_treedepth = max_treedepth)

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
        y[global_iter] = x %*% betas[i,] + stats::rnorm(1, 0, y_sigma)  # outcome
        id[global_iter] = i                                  # subject ID
        period[global_iter] = p                              # period

        global_iter = global_iter + 1
      }

    } # end of subject loop

    # Append data to ongoing trial information
    period_data = list(J = n_subj, K = n_trts, N = nn, X = X, y = y, id = id, period = period)
    period_data_tidy = tidy_data(period_data)
    trial_data = dplyr::bind_rows(trial_data, period_data_tidy)

    current_data = list(J = n_subj, K = n_trts,
                        N = nrow(trial_data),
                        X = trial_data %>% dplyr::select(dplyr::contains("X")) %>% as.matrix,
                        y = trial_data$y,
                        id = trial_data$id,
                        period = trial_data$period)

    # Generating posterior samples after data collected for everyone in period
    # Model here gives everyone their own prior, no aggregation
    stan_model = mcmc(stanfile, current_data, chains = chains, warmup = warmup,
                      iter = iter, adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)

  } # end of period loop

  # Final aggregated data on trial to understand population-level effects
  # agg_stan_model = mcmc(stanfile, current_data, chains = chains, warmup = warmup,
  #                       iter = iter, adapt_delta = adapt_delta,
  #                       max_treedepth = max_treedepth)

  # Format allocation probability matrix for long format
  prob_df = tibble::tibble()
  for (i in 1:dim(allocprobs)[1]) {
    period_probs = allocprobs[i,,] %>% tibble::as_tibble()
    colnames(period_probs) = paste0("p", 1:n_trts)
    period_probs$period = i
    period_probs$id = 1:n_subj
    prob_df = dplyr::bind_rows(prob_df, period_probs)
  }

  list(
    n_subj = n_subj,
    n_periods = n_periods,
    n_obvs = n_obvs,
    betas = betas,
    y_sigma = y_sigma,
    data = trial_data,
    stan_model = stan_model,
    # agg_stan_model = agg_stan_model,
    allocprobs = prob_df
  )

}
