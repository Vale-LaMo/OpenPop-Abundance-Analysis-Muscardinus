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
    phi[k] = exp(-mu * delta[k]); // Derived survival probabilities
}

model {
  // Weakly Informative Priors
  mu ~ normal(0, 1);
  psi_area ~ beta(1, 1);
  
  for (s in 1:N_areas) {
    alpha_p[s] ~ normal(0, 2);
    for (k in 1:K) {
      beta[s, k] ~ beta(1, 1);
    }
  }

  for (i in 1:M_total) {
    int s = area_idx[i]; // Current area index
    
    real lp0 = log1m(psi_area[s]);
    for (k in 1:K) {
      if (y[i, k] == 1) lp0 = negative_infinity(); // Must be real if observed
    }

    vector[K] log_alpha_alive;
    vector[K] log_alpha_not_entered;
    vector[K] log_alpha_dead;

    // Initialization for k=1
    real p_1 = inv_logit(alpha_p[s, 1]);
    if (y[i, 1] == 1) {
        log_alpha_not_entered[1] = negative_infinity();
    } else {
        log_alpha_not_entered[1] = log(1 - beta[s, 1]);
    }
    log_alpha_alive[1] = log(beta[s, 1]) + bernoulli_lpmf(y[i, 1] | p_1);
    log_alpha_dead[1] = negative_infinity();

    // Forward Algorithm loop
    for (k in 2:K) {
        real p_k = inv_logit(alpha_p[s, k]);
        real surv_prev = phi[k-1];
        
        if (y[i, k] == 1) {
            log_alpha_not_entered[k] = negative_infinity();
            log_alpha_dead[k] = negative_infinity();
        } else {
            log_alpha_not_entered[k] = log_alpha_not_entered[k-1] + log(1 - beta[s, k]);
            log_alpha_dead[k] = log_sum_exp(log_alpha_dead[k-1], log_alpha_alive[k-1] + log1m(surv_prev));
        }

        log_alpha_alive[k] = log_sum_exp(
            log_alpha_alive[k-1] + log(surv_prev),
            log_alpha_not_entered[k-1] + log(beta[s, k])
        ) + bernoulli_lpmf(y[i, k] | p_k);
    }

    // Fixed the initialization of final_log_probs
    array[3] real final_log_probs;
    final_log_probs[1] = log_alpha_alive[K];
    final_log_probs[2] = log_alpha_not_entered[K];
    final_log_probs[3] = log_alpha_dead[K];
    
    target += log_sum_exp(lp0, log(psi_area[s]) + log_sum_exp(final_log_probs));
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
