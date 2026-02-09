# Function to work with symptom compartments model and screening
# use this function when pilgrims are travelling to KSA and we want screening to occur

update_ksa_state_symptom_screening_incomingphase <- function(
    state,
    old_state,
    dt,
    moving_compartments = c("susceptible",
                            "exposed",
                            "infected_asymptomatic",
                            "infected_presymptomatic",
                            "infected_symptomatic",
                            "recovered"),
    screening_compartments = c("exposed_diagnosed",
                               "infected_asymptomatic_diagnosed",
                               "infected_presymptomatic_diagnosed",
                               "infected_symptomatic_diagnosed"),
    false_positive_compartments = c("susceptible_false_positive",
                                    "recovered_false_positive"),
    movement_type = c("probability", "rate"), # set default movement type to be prob
    relative_movement = c(1, 1, 1, 1, 1, 1),
    background_testing_rate = 0.61,
    symptomatic_testing_rate = 0.87,
    ksa_index, atrisk_index, end_of_isolation_probabilities
) {
  
  movement_type <- match.arg(movement_type)
  
  if (!is.null(old_state)) {
    
    # Get the numbers of people who entered the false positive compartments at t-isolation_period  
    finished_isolating_s_all_patches <- unlist(sapply(old_state$patches, function(x) x$new_susceptible_false_positive))
    finished_isolating_r_all_patches <- unlist(sapply(old_state$patches, function(x) x$new_recovered_false_positive))
    
    # Get the numbers of people who entered the true positive compartments at t-isolation_period
    finished_isolating_list <- list(
      entered_e = unlist(sapply(old_state$patches, function(x) x$new_exposed_diagnosed)),
      entered_ia = unlist(sapply(old_state$patches, function(x) x$new_infected_asymptomatic_diagnosed)),
      entered_ip = unlist(sapply(old_state$patches, function(x) x$new_infected_presymptomatic_diagnosed)),
      entered_is = unlist(sapply(old_state$patches, function(x) x$new_infected_symptomatic_diagnosed))
    )
  }
  
  n_moving <- get_number_migrating_symptoms(state, dt, moving_compartments, movement_type, relative_movement)
  
  # Symptomatic movers arriving in KSA are tested
  # a proportion (equal to symptomatic_testing_rate) are tested
  # calculate this by performing binomial draw for each element in matrix
  tested_compartments <- c("exposed", "infected_asymptomatic",
                           "infected_presymptomatic", "infected_symptomatic")
  
  movers_in_tested_compartments <- n_moving[tested_compartments]
  
  # Define testing rate and sensitivity
  # test_rate <- state[["patches"]][[1]][["testing_rate"]]
  test_sensitivity <- state[["patches"]][[1]][["test_sensitivity"]]
  
  # for each set of movers, draw from binomial distribution to get number that would be tested
  tested_on_arrival <- lapply(names(movers_in_tested_compartments), function(name) {
    mat <- movers_in_tested_compartments[[name]]
    
    # Set the test rate based on the name of the list element
    # test_rate <- if (name == "infected_symptomatic") {
    test_rate <- if (name %in% tested_compartments) {
      symptomatic_testing_rate
    } else {
      background_testing_rate
    }
    
    # Apply rbinom function to each element in the matrix
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_rate))
  })
    
  names(tested_on_arrival)  <- names(movers_in_tested_compartments)
    
  # for each set of tested ppl, draw from binomial distribution to get number that would be diagnosed
  diagnosed_on_arrival <- lapply(tested_on_arrival, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_sensitivity))
  })
  
  # Record how many false negatives there were
  missed_diagnosis <- purrr::map2(tested_on_arrival, diagnosed_on_arrival, \(x, y) x-y)
  
  # We can also get some pilgrims in S or R who are falsely diagnosed on arrival
  # We assume a certain false positive rate that depends on the test specificity and
  # apply a binomial draw as above
  
  compartment_sources_of_false_positives <- c("susceptible", "recovered")
  
  movers_in_false_pos_compartments <- n_moving[compartment_sources_of_false_positives]
  
  # single test specificity in all patches at the moment so we just extract from one of the patches
  test_specificity <- state[["patches"]][[1]][["test_specificity"]]
  false_positive_rate <- 1 - test_specificity
  
  # for each set of movers, draw from binomial distribution to get number that would be tested
  s_or_r_tested_on_arrival <- lapply(movers_in_false_pos_compartments, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, background_testing_rate))
  })
  
  falsely_diagnosed_on_arrival <- lapply(s_or_r_tested_on_arrival, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, false_positive_rate))
  })
  
  # Record how many true negatives there were
  true_negatives <- purrr::map2(s_or_r_tested_on_arrival,
                                falsely_diagnosed_on_arrival, \(x, y) x-y)
  
  # test_that("S test results = number of tested S",
  #           expect_identical(s_or_r_tested_on_arrival[[1]],
  #                            falsely_diagnosed_on_arrival[[1]] +
  #                                true_negatives[[1]]))
  
  diagnoses_on_arrival <- c(diagnosed_on_arrival, falsely_diagnosed_on_arrival)
  
  n_patches <- length(state[["patches"]])
  
  # Step 1. Move individuals
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    for (compartment in moving_compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx) -
        from_other_patches(diagnoses_on_arrival[[compartment]], idx) # diagnosed cases will go to separate compartments
      
      imported_compartment <- paste0("imported_", compartment)
      patch[[imported_compartment]] <- from_other_patches(n_moving[[compartment]], idx)
      
    }
    
    for (compartment in screening_compartments) {
      undiagnosed_compartment <- sub("_diagnosed$", "", compartment)
      # patch[[compartment]] <- patch[[compartment]] +
      #   from_other_patches(diagnoses_on_arrival[[undiagnosed_compartment]], idx) # screened cases
      
      new_diagnosed <- paste0("new_", compartment)
      patch[[new_diagnosed]] <- from_other_patches(diagnoses_on_arrival[[undiagnosed_compartment]], idx)
      
      new_false_neg <- paste0("new_false_neg_", undiagnosed_compartment)
      patch[[new_false_neg]] <- from_other_patches(missed_diagnosis[[undiagnosed_compartment]], idx)
      
    }
    
    # Trialing using exposed_diagnosed as a compartment to represent all isolating pilgrims
    
    # Identify elements that start with "new_" and end with "_diagnosed"
    elements_to_sum <- grep("^new_.*_diagnosed$", names(patch), value = TRUE)
    
    # Sum the values of these elements
    sum_of_elements <- sum(unlist(patch[elements_to_sum]))
    
    patch[["all_diagnosed"]] <- patch[["all_diagnosed"]] +
      sum_of_elements
    
    for (compartment in false_positive_compartments) {
      unscreened_compartment <- sub("_false_positive$", "", compartment)
      patch[[compartment]] <- patch[[compartment]] +
        from_other_patches(diagnoses_on_arrival[[unscreened_compartment]], idx) # screened cases
      
      new_false_positive <- paste0("new_", compartment)
      patch[[new_false_positive]] <- from_other_patches(diagnoses_on_arrival[[unscreened_compartment]], idx)
      
      new_true_neg <- paste0("new_true_neg_", unscreened_compartment)
      patch[[new_true_neg]] <- from_other_patches(
        true_negatives[[unscreened_compartment]], idx)
      
    }
    
    state[["patches"]][[idx]] <- patch
    
  }
  
  # Step 2. Update disease states.
  
  # Compute the exposure rates among pilgrims and "at risk" non-pilgrims  
  pilgrim_exposure_rate <- compute_ksa_exposure_rate(state, ksa_index, atrisk_index)
  atrisk_exposure_rate <- compute_atrisk_exposure_rate(state, ksa_index, atrisk_index)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    finished_isolating_s <- set_finished_isolators(old_state, finished_isolating_s_all_patches, idx)
    finished_isolating_r <- set_finished_isolators(old_state, finished_isolating_r_all_patches, idx)
    
    finished_isolating_infected <- set_finished_isolators_infected(old_state,
                                                                   finished_isolating_list,
                                                                   end_of_isolation_probabilities,
                                                                   idx)
    
    if (idx %in% ksa_index) {
      # this first modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      # if (!is.null(old_state)) browser()
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, pilgrim_exposure_rate,
                                                             finished_isolating_s,
                                                             finished_isolating_r,
                                                             finished_isolating_infected,
                                                             old_state,
                                                             screening = TRUE)
    } else if (idx %in% atrisk_index) {
      
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, atrisk_exposure_rate,
                                                             finished_isolating_s,
                                                             finished_isolating_r,
                                                             finished_isolating_infected,
                                                             old_state,
                                                             screening = TRUE)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt,
                                                         screening = TRUE)
    }
  }
  state
  
}
