set_finished_isolators <- function(old_state,
                                   finished_isolating_numbers,
                                   idx) {
  
  if (is.null(old_state)) {
    
    finished_isolating <- 0
    
  } else {
    
    finished_isolating <- finished_isolating_numbers[idx]
    
  }
  
  finished_isolating
  
}

# finished_isolating_numbers argument should be a list of the different compartments in the
# isolated model section.
set_finished_isolators_infected <- function(old_state,
                                   finished_isolating_numbers,
                                   probability_vectors,
                                   idx) {
  
  if (is.null(old_state)) {
    
    finished_isolating <- rep(0, 5) # vector of zeros for each of the possible end of isolation states
    
  } else {
    
    if (!is.list(finished_isolating_numbers)) {
      stop("The argument finished_isolating_numbers must be a list.")
    } 
    
    # Finished_isolating_numbers is a list with the number of people entering each true positive
    # isolation compartment at the start of the isolation period.
    # For each of those starting compartments, probability_vectors gives the probability that an
    # individual finishes their isolation in each of the possible final isolation compartments. 
    finished_isolating <- map2(finished_isolating_numbers, probability_vectors,
                               function(x, y) {
      
      entered_isolation <- x[idx]
      
      # entered_isolation gives the numbers for the different starting compartments
      # we extract the corresponding finish probability vectors from y
      out <- as.vector(stats::rmultinom(1, entered_isolation, y))
      names(out) <- names(y)
      out

      ##DO NEXT: NEEDS FURTHER CHECKING!! PASSING THROUGH THE PROBABILITY VECTORS
    })
    
    finished_isolating <- Reduce(`+`, finished_isolating)
  }
  
  finished_isolating
  
}

# Separate function for determining the finished isolators in the diagnosed compartments
# In comparison to the isolators in the false_positive, the diagnosed isolators move through
# infection stages over the course of their isolation. Therefore we need to extract them from
# potentially later compartments.
distribute_finished_isolators_infected <- function(finished_isolating_numbers,
                                                   probability_vector) {
  
  as.vector(stats::rmultinom(1, finished_isolating_numbers, probability_vector))
  
  
}


simulate_isolation_period <- function(starting_compartment, input_state,
                                      sample_size, isolation_period) {
  
  # input_state should be a simple test state generated using state_symptoms_testing()
  # there should be a single patch
  # values for all compartments should be zero, and we then assign sample_size to starting_compartment
  test_starting_state <- input_state
  test_starting_state$patches[[1]][[starting_compartment]] <- sample_size
  
  nsim <- 1
  
  # Re-using existing code for multiple simulations so some of this is a little unnecessary, but works with later processing functions
  out <- vector(mode = "list", length = nsim)
  ## Run multiple simulations.
  ## Each simulation outputs a list of length dt
  ## Each element in this list is a state object.
  for (sim in seq_len(nsim)) {
    this_sim <- vector(mode = "list", length = isolation_period)
    this_sim[[1]] <- test_starting_state[[1]][[1]]
    
      for (time in 2:isolation_period) {
        print(time)
        this_sim[[time]] <- diagnosed_transitions(
          state = this_sim[[time - 1]],
          dt = 1
        )
      }
    
    out[[sim]] <- this_sim
    
  }
  
  out_df <- map2_dfr(out, seq_along(out), ~ {
    simulation_number <- .y
    time_list <- .x
    map2_dfr(time_list, seq_along(time_list), ~ {
      time <- .y
      variable_list <- .x
      tibble(
        simulation_number = simulation_number,
        time = time,
        final_compartment = names(variable_list),
        value = unlist(variable_list, use.names = FALSE)
      )
    }) %>% 
      filter(final_compartment %in% c("exposed_diagnosed", "infected_asymptomatic_diagnosed",
                             "infected_presymptomatic_diagnosed", "infected_symptomatic_diagnosed",
                             "recovered_diagnosed"))
  })
  
  out_probabilities <- out_df %>% 
    filter(time == isolation_period) %>% 
    group_by(final_compartment) %>% 
    summarise(total = sum(value)) %>% 
    mutate(probability = total/sample_size,
           starting_compartment = starting_compartment) %>% 
    select(starting_compartment, final_compartment, total, probability)
  
  probabilities_as_vector <- out_probabilities %>%
    select(final_compartment, probability) %>% deframe()
  
  probabilities_as_vector
  
}

diagnosed_transitions <- function(state, dt) {
  # browser()
  newly_infected_diagnosed <-  to_next_compartment(
    state$exposed_diagnosed, state$infection_rate, dt
  )
  
  newly_infected_diag_presymptomatic <- stats::rbinom(1, newly_infected_diagnosed, state$prop_symptomatic)
  newly_infected_diag_asymptomatic <- newly_infected_diagnosed - newly_infected_diag_presymptomatic
  
  state$exposed_diagnosed <- state$exposed_diagnosed -
    newly_infected_diagnosed -
    deaths(state$exposed_diagnosed, state$death_rate, dt)
  
  newly_recovered_diag_asymptomatic <-  to_next_compartment(
    state$infected_asymptomatic_diagnosed, state$recovery_rate_asym, dt
  ) 
  
  state$infected_asymptomatic_diagnosed <- state$infected_asymptomatic_diagnosed -
    newly_recovered_diag_asymptomatic -
    deaths(state$infected_asymptomatic_diagnosed, state$death_rate, dt) +
    newly_infected_diag_asymptomatic
  
  newly_diag_symptomatic <- to_next_compartment(
    state$infected_presymptomatic_diagnosed, state$symptom_rate, dt
  )
  
  state$infected_presymptomatic_diagnosed <- state$infected_presymptomatic_diagnosed -
    newly_diag_symptomatic -
    deaths(state$infected_presymptomatic_diagnosed, state$death_rate, dt) +
    newly_infected_diag_presymptomatic
  
  newly_recovered_diag_symptomatic <-  to_next_compartment(
    state$infected_symptomatic_diagnosed, state$recovery_rate_sym, dt
  )
  
  state$infected_symptomatic_diagnosed <- state$infected_symptomatic_diagnosed -
    newly_recovered_diag_symptomatic -
    deaths(state$infected_symptomatic_diagnosed, state$death_rate, dt) +
    newly_diag_symptomatic
  
  state$recovered_diagnosed <- state$recovered_diagnosed -
    deaths(state$recovered_diagnosed, state$death_rate, dt) +
    newly_recovered_diag_asymptomatic +
    newly_recovered_diag_symptomatic
  
  state
  
}