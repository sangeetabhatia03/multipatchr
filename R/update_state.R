##' @title Update state
##' @param state state is a collection of patches and a matrix of
##' rates of movement between patches.
##' @param dt time for which state should be updated. It is the user's
##' responsibility to make sure that this number is consistent
##' with the units on rates. For instance, if the various rates are
##' per week, dt is assumed to be dt weeks.
##' @param compartments in case they are different from SEIR
##' @param movement_type select whether the movement matrix input
##' contains rates (with the caveat on appropriate units from
##' above remaining valid) or probabilities.
##' @param relative_movement vector specifying the relative movement
##' in each infection compartment. By default all values are set to 1.
##' We can state that movement in a given compartment should be reduced
##' to 10 percent the normal level by replacing that vector element with 0.1.
##' @return state updated
##' @author Sangeeta Bhatia
##' @export
update_state <- function(state,
                         dt,
                         compartments = c("susceptible",
                                          "exposed",
                                          "infected",
                                          "recovered"),
                         movement_type = c("probability", "rate"), # set default movement type to be prob
                         relative_movement = c(1, 1, 1, 1)
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating(state, dt, compartments, movement_type, relative_movement)
  n_patches <- length(state[["patches"]])
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    for (compartment in compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx)
      
    }
    
    # state[["n_moving"]] <- n_moving
    
    state[["patches"]][[idx]] <- update_patch(patch, dt)
    
  }
  state
}

# modify function to work this KSA example
update_ksa_state <- function(state,
                             dt,
                             compartments = c("susceptible",
                                              "exposed",
                                              "infected",
                                              "recovered"),
                             movement_type = c("probability", "rate"), # set default movement type to be prob
                             relative_movement = c(1, 1, 1, 1),
                             ksa_index
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating(state, dt, compartments, movement_type, relative_movement)
  n_patches <- length(state[["patches"]])
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many infections in KSA at this time, across all sub-patches
  ksa_infections <- ksa_patches %>% 
    purrr::map_dbl(~ .x$infected) %>% 
    purrr::reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    purrr::map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                  pluck(.x, "infected"), pluck(.x, "recovered"))) %>%
    purrr::reduce(`+`)
  
  # What is ksa exposure rate, calculated across all sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * ksa_infections / ksa_total,
                              0)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    for (compartment in compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx)
      
    }
    
    if (idx %in% ksa_index) {
      # this modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      state[["patches"]][[idx]] <- update_ksa_patch(patch, dt, ksa_exposure_rate)
    } else {
      state[["patches"]][[idx]] <- update_patch(patch, dt)
    }
    
  }
  state
}


# modify function to work with symptom compartments model
update_ksa_state_symptoms <- function(state,
                                      dt,
                                      compartments = c("susceptible",
                                                       "exposed",
                                                       "infected_asymptomatic",
                                                       "infected_presymptomatic",
                                                       "infected_symptomatic",
                                                       "recovered"),
                                      movement_type = c("probability", "rate"), # set default movement type to be prob
                                      relative_movement = c(1, 1, 1, 1, 1, 1),
                                      ksa_index
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating_symptoms(state, dt, compartments, movement_type, relative_movement)
  
  n_patches <- length(state[["patches"]])
  
  ksa_exposure_rate <- compute_ksa_exposure_rate(state, ksa_index)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    # Step 1. Move individuals
    for (compartment in compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx)
      
    }
    
    # Step 2. Update disease states
    if (idx %in% ksa_index) {
      # this modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, ksa_exposure_rate)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt)
    }
    
  }
  state
}

# modify function to work with symptom compartments model and screening
# use this function when pilgrims are travelling to KSA and we want screening to occur

update_ksa_state_screening_incomingphase <- function(state,
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
  
  # Movers arriving in KSA are tested
  # a proportion (equal to testing_rate) test positive
  # calculate this by performing binomial draw for each element in matrix
  ## SB: Make this an argument so that users can specify which compartments are 
  ## tested and which are not. This will allow for more flexibility in the model.
  tested_compartments <- c("exposed", "infected_asymptomatic",
                           "infected_presymptomatic", "infected_symptomatic")
  
  movers_in_tested_compartments <- n_moving[tested_compartments]
  
  # single testing rate and sensitivity at the moment so we just extract from one of the patches
  ## SB: Do you see a scenario where we would want to have multiple testing rates and sensitivities?
  ## If not, then note in the function documentation that this will be extracted
  ## from the first patch only.
  ## Also add a check here, you can extract from all patches and check if there
  ## is only one unique value, if not, throw a warning.
  test_rate <- state[["patches"]][[1]][["testing_rate"]]
  test_sensitivity <- state[["patches"]][[1]][["test_sensitivity"]]
  
  # for each set of movers, draw from binomial distribution to get number that would be tested
  tested_on_arrival <- lapply(movers_in_tested_compartments, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_rate))
  })
  
  # for each set of tested ppl, draw from binomial distribution to get number that would be diagnosed
  diagnosed_on_arrival <- lapply(tested_on_arrival, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_sensitivity))
  })
  
  # Record how many false negatives there were
  ## SB: Wouldn't tested_on_arrival - diagnosed_on_arrival achieve the same thing?
  ## Sorry if I am talking through my hat here, I am not sure if
  ## tested_on_arrival and diagnosed_on_arrival are matrices or lists of matrices.
  missed_diagnosis <- purrr::map2(tested_on_arrival, diagnosed_on_arrival, \(x, y) x-y)
  
  # We can also get some pilgrims in S or R who are falsely diagnosed on arrival
  # We assume a certain false positive rate that depends on the test specificity and
  # apply a binomial draw as above
  ## SB: Same comment as before. To avoid having an unmanageable number of 
  ## arguments, think about having a function that defines your total model
  ## structure, and read off of that in this function.
  compartment_sources_of_false_positives <- c("susceptible", "recovered")
  
  movers_in_false_pos_compartments <- n_moving[compartment_sources_of_false_positives]
  
  # single test specificity in all patches at the moment so we just extract from one of the patches
  test_specificity <- state[["patches"]][[1]][["test_specificity"]]
  false_positive_rate <- 1 - test_specificity
  
  # for each set of movers, draw from binomial distribution to get number that would be tested
  s_or_r_tested_on_arrival <- lapply(movers_in_false_pos_compartments, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_rate))
  })
  
  falsely_diagnosed_on_arrival <- lapply(s_or_r_tested_on_arrival, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, false_positive_rate))
  })
  
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
    ## SB: You don't want to rely on grepping the name of elements
    ## Could you not have a list of compartnment names here?
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
    }
    
    state[["patches"]][[idx]] <- patch
    
  }
  
  # Step 2. Update disease states.
  ## SB: Could this be moved to a separate function?
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

# use this function when pilgrims are either in KSA or travelling home (no screening occurs)
update_ksa_state_screening_otherphases <- function(state,
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
                                                   movement_type = c("probability", "rate"), # set default movement type to be prob
                                                   relative_movement = c(1, 1, 1, 1, 1, 1),
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
  
  n_patches <- length(state[["patches"]])
  
  
  # Step 1. Move individuals
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    for (compartment in moving_compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx) 
      
      imported_compartment <- paste0("imported_", compartment)
      patch[[imported_compartment]] <- from_other_patches(n_moving[[compartment]], idx)
      
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
