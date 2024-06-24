# Function to compute the probability that a person entering isolation in compartment A, will be in
# compartment B, C, or D at the end of their isolation period
# Uses a numerical simulation approach to determine the probabilities.
# The argument sample_size is the number of simulations we run.


get_end_of_isolation_probabilities <- function(test_starting_state, isolation_period, sample_size) {
  
  test_starting_state <- state_symptoms_testing(
    s_patches = test_starting_state$s_patches,
    e_patches = test_starting_state$e_patches,
    i_a_patches = test_starting_state$i_a_patches,
    i_p_patches = test_starting_state$i_p_patches,
    i_s_patches = test_starting_state$i_s_patches,
    r_patches = test_starting_state$r_patches,
    diag_patches = test_starting_state$diag_patches,
    s_false_patches = test_starting_state$s_false_patches,
    r_false_patches = test_starting_state$r_false_patches,
    birth_rates = test_starting_state$birth_rates,
    death_rates = test_starting_state$death_rates,
    transmission_rates = test_starting_state$transmission_rates,
    transmission_rates_pilgrims = test_starting_state$transmission_rates_pilgrims,
    transmission_rates_pilgrims_and_atrisk = test_starting_state$transmission_rates_pilgrims_and_atrisk,
    asymptomatic_infectiousness = test_starting_state$asymptomatic_infectiousness,
    presymptomatic_infectiousness = test_starting_state$presymptomatic_infectiousness,
    prop_symptomatic = test_starting_state$prop_symptomatic,
    infection_rates = test_starting_state$infection_rates,
    symptom_rates = test_starting_state$symptom_rates,
    recovery_rates_asym = test_starting_state$recovery_rates_asym, # day^-1
    recovery_rates_sym = test_starting_state$recovery_rates_sym,
    testing_rates = test_starting_state$testing_rate,
    test_sensitivity = test_starting_state$test_sensitivity,
    test_specificity = test_starting_state$test_specificity,
    isolation_period = isolation_period,
    movement_rate = test_starting_state$movement_rate
  )
  
  # Define the compartments where it is possible for a person to start their diagnosed/isolation period
  possible_starting_compartments <- c("exposed_diagnosed", "infected_asymptomatic_diagnosed",
                                      "infected_presymptomatic_diagnosed", "infected_symptomatic_diagnosed"
  )
  
  # The code below generates a list of length(possible_starting_compartments).
  # Each element of the list is a named vector of probabilities.
  # The vector values give the probability that a person starting isolation in a
  # given starting compartment (denoted by the element name) will end up in each
  # of the possible final_states.
  end_of_isolation_probabilities <- purrr::map(possible_starting_compartments, function(x) {
    
    simulate_isolation_period(x, test_starting_state,
                              sample_size, isolation_period)
    
  })
  
  names(end_of_isolation_probabilities) <- possible_starting_compartments
  
  end_of_isolation_probabilities
  
}