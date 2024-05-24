# The code in this script should be run as a precursor to running the model.
# Computes a set of probabilites that are then used to determine in what state people leave isolation

isolation_period <- 10
sample_size <- 1000000

test_starting_state <- state_symptoms_testing(
  s_patches = 0,
  e_patches = 0,
  i_a_patches = 0,
  i_p_patches = 0,
  i_s_patches = 0,
  r_patches = 0,
  e_diag_patches = 0,
  i_a_diag_patches = 0,
  i_p_diag_patches = 0,
  i_s_diag_patches = 0,
  r_diag_patches = 0,
  s_false_patches = 0,
  r_false_patches = 0,
  birth_rates = 0,
  death_rates = 0,
  transmission_rates = beta,
  transmission_rates_pilgrims = beta * 1.5,
  transmission_rates_pilgrims_and_atrisk = beta,
  asymptomatic_infectiousness = theta_a,
  presymptomatic_infectiousness = theta_p,
  prop_symptomatic = p,
  infection_rates = delta,
  symptom_rates = gamma_p,
  recovery_rates_asym = gamma_a,
  recovery_rates_sym = gamma_s,
  testing_rate = 0.5,
  test_sensitivity = 0.9,
  test_specificity = 0.9,
  isolation_period = isolation_period,
  movement_rate = matrix(1)
)

possible_starting_compartments <- c("exposed_diagnosed", "infected_asymptomatic_diagnosed",
                                    "infected_presymptomatic_diagnosed", "infected_symptomatic_diagnosed"
)

# This generates a list of length(possible_starting_compartments).
# Each element of the list is a named vector of probabilities.
# The vector values give the probability that a person starting isolation in a
# given starting compartment (denoted by the element name) will end up in each
# of the possible final_states.
end_of_isolation_probabilities <- map(possible_starting_compartments, function(x) {
  
  simulate_isolation_period(x, test_starting_state,
                            sample_size, isolation_period)
  
})

names(end_of_isolation_probabilities) <- possible_starting_compartments

# NEXT: Think about how we sample from these vectors to get the final numbers.