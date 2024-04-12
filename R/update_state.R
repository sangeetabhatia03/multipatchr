# Function for computing the ksa_exposure rate
compute_ksa_exposure_rate <- function(state, ksa_index) {
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many symptomatic, pre-symptomatic, and asymptomatic
  # infections in KSA at this time, across all sub-patches
  ksa_infections_symptomatic <- ksa_patches %>%
    purrr::map_dbl(~ .x$infected_symptomatic) %>% 
    purrr::reduce(`+`)
  
  ksa_infections_presymptomatic <- ksa_patches %>% 
    purrr::map_dbl(~ .x$infected_presymptomatic) %>% 
    purrr::reduce(`+`)
  
  ksa_infections_asymptomatic <- ksa_patches %>% 
    purrr::map_dbl(~ .x$infected_asymptomatic) %>% 
    purrr::reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    purrr::map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                  pluck(.x, "infected_asymptomatic"), pluck(.x, "infected_presymptomatic"),
                  pluck(.x, "infected_symptomatic"), pluck(.x, "recovered"))) %>%
    purrr::reduce(`+`)
  
  # What is ksa exposure rate, calculated across all sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  asymptomatic_infectivity <- ksa_patches[[1]]$asymptomatic_infectiousness
  presymptomatic_infectivity <- ksa_patches[[1]]$presymptomatic_infectiousness
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * (
                                ksa_infections_symptomatic + 
                                  presymptomatic_infectivity * ksa_infections_presymptomatic +
                                  asymptomatic_infectivity * ksa_infections_asymptomatic
                              )  / ksa_total,
                              0)
  
}



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
                                                     ksa_index
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating_symptoms(state, dt, moving_compartments, movement_type, relative_movement)
  
  # Movers arriving in KSA are tested
  # a proportion (equal to testing_rate) test positive
  # calculate this by performing binomial draw for each element in matrix
  tested_compartments <- c("exposed", "infected_asymptomatic",
                           "infected_presymptomatic", "infected_symptomatic")
  
  movers_in_tested_compartments <- n_moving[tested_compartments]
  
  # single testing rate at the moment so we just extract from one of the patches
  test_rate <- state[["patches"]][[1]][["testing_rate"]]
  
  # for each set of movers, draw from binomial distribution to get number that would be diagnosed
  diagnosed_on_arrival <- lapply(movers_in_tested_compartments, function(mat) {
    apply(mat, c(1, 2), function(x) stats::rbinom(1, x, test_rate))
  })
  
  # Create a list that also includes the compartments that were not tested (S, R, etc)
  # use this to record which of the movers need to go to test compartments
  diagnosed_all <- lapply(n_moving, function(mat) {
    apply(mat, c(1, 2), function(x) 0L)
  })
  
  for (name in names(diagnosed_on_arrival)) {
    diagnosed_all[[name]] <- diagnosed_on_arrival[[name]]
  }
  
  n_patches <- length(state[["patches"]])
  
  ksa_exposure_rate <- compute_ksa_exposure_rate(state, ksa_index)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    # Step 1. Move individuals
    for (compartment in moving_compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx) -
        from_other_patches(diagnosed_all[[compartment]], idx) # diagnosed cases will go to separate compartments
    }
    
    for (compartment in screening_compartments) {
      undiagnosed_compartment <- sub("_diagnosed$", "", compartment)
      patch[[compartment]] <- patch[[compartment]] +
        from_other_patches(diagnosed_all[[undiagnosed_compartment]], idx) # screened cases
    }
    
    # Step 2. Update disease states.
    if (idx %in% ksa_index) {
      # this first modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, ksa_exposure_rate,
                                                             screening = TRUE)
      # state[["patches"]][[idx]] <- update_ksa_patch_screening(patch, dt)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt,
                                                         screening = TRUE)
    }
    
  }
  state
}

# use this function when pilgrims are either in KSA or travelling home (no screening occurs)
update_ksa_state_screening_otherphases <- function(state,
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
                                                   ksa_index
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating_symptoms(state, dt, moving_compartments, movement_type, relative_movement)
  
  # MASK FOR NOW - DELETE LATER
  # # Movers arriving in KSA are tested
  # # a proportion (equal to testing_rate) test positive
  # # calculate this by performing binomial draw for each element in matrix
  # tested_compartments <- c("exposed", "infected_asymptomatic",
  #                          "infected_presymptomatic", "infected_symptomatic")
  # 
  # movers_in_tested_compartments <- n_moving[tested_compartments]
  # 
  # # single testing rate at the moment so we just extract from one of the patches
  # test_rate <- state[["patches"]][[1]][["testing_rate"]]
  # 
  # # for each set of movers, draw from binomial distribution to get number that would be diagnosed
  # diagnosed_on_arrival <- lapply(movers_in_tested_compartments, function(mat) {
  #   apply(mat, c(1, 2), function(x) rbinom(1, x, test_rate))
  # })
  # 
  # # Create a list that also includes the compartments that were not tested (S, R, etc)
  # # use this to record which of the movers need to go to test compartments
  # diagnosed_all <- lapply(n_moving, function(mat) {
  #   apply(mat, c(1, 2), function(x) 0L)
  # })
  # 
  # for (name in names(diagnosed_on_arrival)) {
  #   diagnosed_all[[name]] <- diagnosed_on_arrival[[name]]
  # }
  # 
  
  n_patches <- length(state[["patches"]])
  
  ksa_exposure_rate <- compute_ksa_exposure_rate(state, ksa_index)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]
    
    for (compartment in moving_compartments) {
      patch[[compartment]] <- patch[[compartment]] -
        to_other_patches(n_moving[[compartment]],  idx) +
        from_other_patches(n_moving[[compartment]], idx) #-
      # from_other_patches(diagnosed_all[[compartment]], idx) # diagnosed cases will go to separate compartments
      
    }
    
    # for (compartment in screening_compartments) {
    #   undiagnosed_compartment <- sub("_diagnosed$", "", compartment)
    #   patch[[compartment]] <- patch[[compartment]] +
    #     from_other_patches(diagnosed_all[[undiagnosed_compartment]], idx) # screened cases
    # }
    
    if (idx %in% ksa_index) {
      # this first modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, ksa_exposure_rate,
                                                             screening = TRUE)
      # state[["patches"]][[idx]] <- update_ksa_patch_screening(patch, dt)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt,
                                                         screening = TRUE)
    }
    
  }
  state
}
