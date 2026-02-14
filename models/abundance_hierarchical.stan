data {
  int<lower=1> N_areas;               // Number of study areas (e.g., 6)
  int<lower=1> K;                     // Number of capture occasions
  int<lower=1> M_total;               // Total individuals across all areas (augmented)
  array[M_total, K] int<lower=0, upper=1> y; 
  array[M_total] int<lower=1, upper=N_areas> area_idx; // Area ID for each individual
  array[K - 1] real<lower=0> delta;   // Time intervals between sessions
}

parameters {
  // Global Hyper-parameters (shared across areas)
  real<lower=0> mu;                   // Global mortality rate
  
  // Area-specific parameters
  vector<lower=0, upper=1>[N_areas] psi_area; 
  array[N_areas, K] real alpha_p;     // Detection intercepts (logit scale)
  array[N_areas, K] real<lower=0, upper=1> beta; // Entry probabilities
}

transformed parameters {
  vector<lower=0, upper=1>[K - 1] phi; 
  for (k in 1:(K - 1))
    phi[k] = exp(-mu * delta[k]);     // Derived survival probabilities - besed on mu and taking into account
                                      // the exact time interval for the session
}

model {
  // Weakly Informative Priors
  mu ~ normal(0, 1);
  psi_area ~ beta(1, 1);
  
  // the structure with nested for loops reflects the spatial heterogeneity
  for (s in 1:N_areas) { 
    // each area has its own history of detectability and recruitment (the model does not force Area 1
    // to have the same detection probability of Area 6 and so on)
    alpha_p[s] ~ normal(0, 2);
    for (k in 1:K) {
      // within each area, the entry prob beta can be different at each session
      // so the model can capture peaks in recruitment at specific time points
      beta[s, k] ~ beta(1, 1);
    }
  }

  for (i in 1:M_total) {
    int s = area_idx[i]; // Current area index
    
    real lp0 = log1m(psi_area[s]);  // prob that an individual in the augmented dataset isn't part
                                    // of the real population
    for (k in 1:K) {
      if (y[i, k] == 1) lp0 = negative_infinity(); // Must be real if observed:
      // if an individual is seen at least once, lp0 = -inf, forcing the model to consider it real
    }

    // the three 'paths' for the status of an individual
    vector[K] log_alpha_alive;
    vector[K] log_alpha_not_entered; // not born or not immigrated yet
    vector[K] log_alpha_dead;

    // Initialization of the forward algorithm for k=1
    // in a state-space model, we define 'where' an individual is before calculating transitions
    real p_1 = inv_logit(alpha_p[s, 1]); // probability of detection of ind i in area s, first session
    if (y[i, 1] == 1) {
        log_alpha_not_entered[1] = negative_infinity();   // if seen, it can't be 'not entered', prob=0
    } else {
        log_alpha_not_entered[1] = log(1 - beta[s, 1]);   // if not seen, the prob not entered is 
                                                          // the complement to the entry prob
    }
    log_alpha_alive[1] = log(beta[s, 1]) + bernoulli_lpmf(y[i, 1] | p_1); // prob alive: an ind is alive
    // if it entered the system with prob beta and its observation id coherent with the detection probability
    // we are summing because we are on the log scale, multiply if probs
    log_alpha_dead[1] = negative_infinity(); // dead prob = 0; can't be dead if you weren't there :)

    // Forward Algorithm loop - change of status through time
    // starts with the second occazion
    for (k in 2:K) {
        real p_k = inv_logit(alpha_p[s, k]); // session specific detection prob
        real surv_prev = phi[k-1]; // survival prob that has been calculated in the transformed par block
        // because it is based on the global mortality and the temporal interval from the prev session
        
        if (y[i, k] == 1) { // if seen, cannot be 'not entered' nor 'dead'
            log_alpha_not_entered[k] = negative_infinity();
            log_alpha_dead[k] = negative_infinity();
        } else { // if not seen
            log_alpha_not_entered[k] = log_alpha_not_entered[k-1] + log(1 - beta[s, k]);
            // the prob of being not entered is the prob of remaining not entered with respect to
            // the previous k-1 session, and the prob of not entrying (1-beta)
            log_alpha_dead[k] = log_sum_exp(log_alpha_dead[k-1], log_alpha_alive[k-1] + log1m(surv_prev));
            // the prob of being dead is the prob of already been dead at k-1, or to die in between;
            // log_sum_exp used to sum this tow possibilities because they are mutually exclusive
        }

        log_alpha_alive[k] = log_sum_exp(
            log_alpha_alive[k-1] + log(surv_prev),        // an ind is alive if it was alive before (k-1) and it survived
            log_alpha_not_entered[k-1] + log(beta[s, k])  // it wasn't entered in k-1 but it enters now
        ) + bernoulli_lpmf(y[i, k] | p_k); // the outcome of the state is weighted with the likelihood of the observation process
        // that is the likelihood is updated at every step taking into account the coherence between estimated
        // state and the capture history
    }

    // Fixed the initialization of final_log_probs
    array[3] real final_log_probs; // contains the final probs calculated by the Forward Alg. at session K (last one)
    final_log_probs[1] = log_alpha_alive[K];
    final_log_probs[2] = log_alpha_not_entered[K];
    final_log_probs[3] = log_alpha_dead[K];
    
    target += log_sum_exp(lp0, log(psi_area[s]) + log_sum_exp(final_log_probs));
    // update of the total log-likelihood: sums the probability of mutually exclusive scenarios
    // lp0 - the ind does not exist (it is in the augmented pop but never entered the real one)
    // lp1 - given by log(psi_area[s]) + log_sum_exp(final_log_probs)
  }
}

generated quantities {
  matrix[N_areas, K] N = rep_matrix(0, N_areas, K); // Multi-site abundance matrix

  for (i in 1:M_total) {
    int s = area_idx[i];
    vector[K] q_alive;
    real q_not_entered_curr = psi_area[s] * (1 - beta[s, 1]);
    
    q_alive[1] = psi_area[s] * beta[s, 1];
    N[s, 1] += q_alive[1];

    for (k in 2:K) {
      q_alive[k] = q_alive[k-1] * phi[k-1] + q_not_entered_curr * beta[s, k];
      q_not_entered_curr *= (1 - beta[s, k]);
      N[s, k] += q_alive[k]; // Expected abundance
    }
  }
}
