# The code in this script should be run as a precursor to running the model.
# Computes a set of probabilities that are then used to determine in what state people leave isolation

test_model_parameters <- list(birth_rates = 0,
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
                              isolation_period = 10,
                              movement_rate = matrix(1)
)

test_initial_states <- list(s_patches = 0,
                            e_patches = 0,
                            i_a_patches = 0,
                            i_p_patches = 0,
                            i_s_patches = 0,
                            r_patches = 0,
                            diag_patches = 0,
                            s_false_patches = 0,
                            r_false_patches = 0)

test_starting_state <- c(test_initial_states, test_model_parameters)


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
  end_of_isolation_probabilities <- map(possible_starting_compartments, function(x) {
    
    simulate_isolation_period(x, test_starting_state,
                              sample_size, isolation_period)
    
  })
  
  names(end_of_isolation_probabilities) <- possible_starting_compartments
  
  end_of_isolation_probabilities
  
}