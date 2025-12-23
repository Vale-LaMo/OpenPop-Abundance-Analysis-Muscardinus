data {
  int<lower=1> M;                 // Number of individuals (observed + augmented)
  int<lower=1> K;                 // Number of capture occasions
  array[M, K] int<lower=0, upper=1> y; // Capture histories
  array[K - 1] real<lower=0> delta; // Time intervals between sessions
}

parameters {
  real<lower=0> mu;               // Mortality rate
  array[K] real alpha_p;          // Detection intercepts per session (on logit scale)
  array[K] real<lower=0, upper=1> beta; // Entry probability per session
  real<lower=0, upper=1> psi;     // Inclusion probability
}

transformed parameters {
  vector<lower=0, upper=1>[K - 1] phi; // Survival probabilities between sessions
  for (k in 1:(K - 1))
    phi[k] = exp(-mu * delta[k]);
}

model {
  // Weakly Informative Priors
  psi ~ beta(1, 1);        // Uniform(0,1) for inclusion probability
  for (k in 1:K) {
    beta[k] ~ beta(1, 1);  // Uniform(0,1) for entry probability per session
  }
  mu ~ normal(0, 1);       // Normal centered at 0, SD 1. Lower-bound 0 makes it a half-normal.
  alpha_p ~ normal(0, 2); // Normal centered at 0, SD 2.5 for logit-scale intercept.

  for (i in 1:M) {
    // log-probability if individual 'i' is NOT a real individual
    // This is only possible if the individual was NEVER observed.
    real lp0 = log1m(psi);
    for (k in 1:K) {
      if (y[i, k] == 1) {
        lp0 = negative_infinity(); // If observed, it must be a real individual
      }
    }

    // log-probability if individual 'i' IS a real individual (P(real) = psi)
    // The total log-likelihood for a "real" individual will be accumulated here.
    // This part is the sum of log-probabilities of all possible paths (entry time, survival/death)
    // that lead to the observed capture history.

    // Using a forward algorithm to track log-probabilities of states
    // log_alpha_alive[k]: log P(observed history up to k AND real AND alive at k)
    // log_alpha_not_entered[k]: log P(observed history up to k AND real AND not entered by k)
    // log_alpha_dead[k]: log P(observed history up to k AND real AND entered and died by k)

    vector[K] log_alpha_alive;
    vector[K] log_alpha_not_entered;
    vector[K] log_alpha_dead;

    // Initial conditions for k=1
    real p_1 = inv_logit(alpha_p[1]);

    // State: Not yet entered by occasion 1
    // Log-probability of being 'not entered' and not observed at k=1.
    // If y[i,1] == 1, this path is impossible.
    if (y[i,1] == 1) {
        log_alpha_not_entered[1] = negative_infinity();
    } else {
        log_alpha_not_entered[1] = log(1 - beta[1]);
    }

    // State: Alive at occasion 1 (entered at 1)
    // Log-probability of being 'alive' and observed/not-observed at k=1.
    log_alpha_alive[1] = log(beta[1]) + bernoulli_lpmf(y[i,1] | p_1);

    // State: Dead by occasion 1 (impossible, as entry is at 1 or later)
    log_alpha_dead[1] = negative_infinity();


    // Loop for k = 2 to K (propagating state probabilities forward)
    for (k in 2:K) {
        real p_k = inv_logit(alpha_p[k]);
        real surv_prev = phi[k-1]; // Survival from k-1 to k
        real death_prev = 1 - surv_prev; // Death from k-1 to k

        // Calculate log_alpha_not_entered[k]
        // This state is reached if not entered by k-1 AND does not enter at k.
        // It's impossible if the individual was observed at k.
        if (y[i,k] == 1) {
            log_alpha_not_entered[k] = negative_infinity();
        } else {
            log_alpha_not_entered[k] = log_alpha_not_entered[k-1] + log(1 - beta[k]);
        }

        // Calculate log_alpha_alive[k]
        // This state is reached if:
        // 1. Was alive at k-1 AND survived.
        // 2. Was not entered yet at k-1 AND enters now (at k).
        // Then, apply observation likelihood for y[i,k].
        log_alpha_alive[k] = log_sum_exp(
            log_alpha_alive[k-1] + log(surv_prev),      // Path 1: Was alive, survived
            log_alpha_not_entered[k-1] + log(beta[k])   // Path 2: Not entered yet, but enters now
        );
        log_alpha_alive[k] += bernoulli_lpmf(y[i,k] | p_k); // Add observation likelihood given alive

        // Calculate log_alpha_dead[k]
        // This state is reached if:
        // 1. Was already dead at k-1.
        // 2. Was alive at k-1 AND died between k-1 and k.
        // It's impossible if the individual was observed at k.
        if (y[i,k] == 1) {
            log_alpha_dead[k] = negative_infinity();
        } else {
            log_alpha_dead[k] = log_sum_exp(
                log_alpha_dead[k-1],                      // Path 1: Was already dead
                log_alpha_alive[k-1] + log(death_prev)    // Path 2: Was alive, but died
            );
        }
    }

    // The total log-likelihood for a "real" individual is the log_sum_exp of being in any of the
    // possible states (alive, not_entered, dead) at the *final* occasion K,
    // having produced the observed capture history.
    array[3] real final_log_probs; // Declare a temporary array
    final_log_probs[1] = log_alpha_alive[K];
    final_log_probs[2] = log_alpha_not_entered[K];
    final_log_probs[3] = log_alpha_dead[K];

    real log_lik_history_given_real = log_sum_exp(final_log_probs); // Pass the array to log_sum_exp

    // Add to target: The individual is either NOT real (lp0) OR IS real (log(psi) + log_lik_history_given_real)
    target += log_sum_exp(lp0, log(psi) + log_lik_history_given_real);
  }
}

generated quantities {
  array[K] real N;      // Estimated expected abundance per session
  array[M, K] real q_alive; // Probability that individual 'i' is alive at occasion 'k'
  array[M, K] real q_not_entered; // Probability that individual 'i' has not yet entered by occasion 'k'
  // array[M, K] real q_dead; // Probability that individual 'i' has entered and died by occasion 'k'

  for (i in 1:M) {
    // Initialize state probabilities for occasion 1
    q_alive[i,1] = psi * beta[1]; // P(real AND enters at 1)
    q_not_entered[i,1] = psi * (1 - beta[1]); // P(real AND does not enter at 1)
    // q_dead[i,1] = 0; // P(real AND dead at 1) - effectively 0 at first occasion

    for (k in 2:K) {
      real surv_k_minus_1 = phi[k-1]; // Survival from k-1 to k

      // Update P(real AND alive at k)
      // = P(real AND alive at k-1 AND survived)
      // + P(real AND not entered at k-1 AND enters at k)
      q_alive[i,k] = q_alive[i,k-1] * surv_k_minus_1 + q_not_entered[i,k-1] * beta[k];

      // Update P(real AND not yet entered at k)
      // = P(real AND not yet entered at k-1 AND does NOT enter at k)
      q_not_entered[i,k] = q_not_entered[i,k-1] * (1 - beta[k]);

      // Update P(real AND dead at k) - if you need it (not directly for N)
      // q_dead[i,k] = q_dead[i,k-1] + q_alive[i,k-1] * (1 - surv_k_minus_1);
    }
  }

  for (k in 1:K) {
    real sum_q_alive = 0;
    for (i in 1:M) {
      sum_q_alive += q_alive[i,k]; // Sum of probabilities of being alive for all augmented individuals
    }
    N[k] = sum_q_alive; // Expected abundance is the sum of these probabilities
  }
}
