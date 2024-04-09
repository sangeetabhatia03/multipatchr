deaths <- function(n, death_rate, dt) {
  
  stats::rbinom(1, size = n, prob = death_rate * dt)
}

births <- function(n, birth_rate, dt) {
  
  stats::rbinom(1, size = n, prob = birth_rate * dt)
}

to_next_compartment <- function(n_current, rate, dt) {
  
  prob <- 1 - rate_to_probability(rate, dt)
  
  out <- stats::rbinom(1, size = n_current, prob = prob)
  
  if (is.na(out)) {
    
    # arg_values <- lapply(substitute(list(n_current, rate, dt)), deparse)
    # print(arg_values)
    # print(paste0(
      # "Number=", n_current, "; Prob=", prob))
    stop("Value moving to compartment cannot be NA", call. = FALSE)
  }
  
  out
  
}

update_patch <- function(patch, dt) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  # conditional statement means that we can handle patches that become empty
  exposure_rate <- ifelse(patch$susceptible +
                            patch$exposed +
                            patch$infected +
                            patch$recovered > 0,
                          patch$transmission_rate *
                            patch$infected / (patch$susceptible +
                                                patch$exposed +
                                                patch$infected +
                                                patch$recovered),
                          0)
  
  newly_exposed <- to_next_compartment(
    patch$susceptible, exposure_rate, dt
  )
  
  patch$susceptible <- patch$susceptible -
    newly_exposed -
    deaths(patch$susceptible, patch$death_rate, dt) +
    births(patch$susceptible, patch$birth_rate, dt)
  
  newly_infected <-  to_next_compartment(
    patch$exposed, patch$infection_rate, dt
  )
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered <-  to_next_compartment(
    patch$infected, patch$recovery_rate, dt
  )
  patch$infected <- patch$infected -
    newly_recovered -
    deaths(patch$infected, patch$death_rate, dt) +
    newly_infected
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered
  
  patch
}

# Function to allow mixing between all the sub-patches in ksa
update_ksa_patch <- function(patch, dt, patch_exposure_rate) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  exposure_rate <- patch_exposure_rate
  
  newly_exposed <- to_next_compartment(
    patch$susceptible, exposure_rate, dt
  )
  
  patch$susceptible <- patch$susceptible -
    newly_exposed -
    deaths(patch$susceptible, patch$death_rate, dt) +
    births(patch$susceptible, patch$birth_rate, dt)
  
  newly_infected <-  to_next_compartment(
    patch$exposed, patch$infection_rate, dt
  )
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered <-  to_next_compartment(
    patch$infected, patch$recovery_rate, dt
  )
  patch$infected <- patch$infected -
    newly_recovered -
    deaths(patch$infected, patch$death_rate, dt) +
    newly_infected
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered
  
  patch
}

# Modified functions for a model w/ symptom compartments
update_patch_symptoms <- function(patch, dt) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  # conditional statement means that we can handle patches that become empty
  exposure_rate <- ifelse(patch$susceptible +
                            patch$exposed +
                            patch$infected_asymptomatic +
                            patch$infected_presymptomatic +
                            patch$infected_symptomatic +
                            patch$recovered > 0,
                          patch$transmission_rate * (
                            (patch$infected_presymptomatic + patch$infected_symptomatic) +
                              (patch$asymptomatic_infectiousness * patch$infected_asymptomatic)) /
                            (patch$susceptible +
                               patch$exposed +
                               patch$infected_asymptomatic +
                               patch$infected_presymptomatic +
                               patch$infected_symptomatic +
                               patch$recovered
                            ),
                          0)
  
  newly_exposed <- to_next_compartment(
    patch$susceptible, exposure_rate, dt
  )
  
  patch$susceptible <- patch$susceptible -
    newly_exposed -
    deaths(patch$susceptible, patch$death_rate, dt) +
    births(patch$susceptible, patch$birth_rate, dt)
  
  newly_infected <-  to_next_compartment(
    patch$exposed, patch$infection_rate, dt
  )
  
  newly_infected_presymptomatic <- round(patch$prop_symptomatic * newly_infected)
  newly_infected_asymptomatic <- newly_infected - newly_infected_presymptomatic
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic, patch$recovery_rate_asym, dt
  )
  
  patch$infected_asymptomatic <- patch$infected_asymptomatic -
    newly_recovered_asymptomatic -
    deaths(patch$infected_asymptomatic, patch$death_rate, dt) +
    newly_infected_asymptomatic
  
  newly_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic <- patch$infected_presymptomatic -
    newly_symptomatic -
    deaths(patch$infected_presymptomatic, patch$death_rate, dt) +
    newly_infected_presymptomatic
  
  newly_recovered_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic <- patch$infected_symptomatic -
    newly_recovered_symptomatic -
    deaths(patch$infected_symptomatic, patch$death_rate, dt) +
    newly_symptomatic
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered_asymptomatic +
    newly_recovered_symptomatic
  
  patch
}

update_ksa_patch_symptoms <- function(patch, dt, patch_exposure_rate) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  # First apply transitions to diagnosed compartments
  patch <- update_ksa_patch_screening(patch, dt)
  
  # Now apply transitions to undiagnosed compartments
  
  exposure_rate <- patch_exposure_rate
  
  newly_exposed <- to_next_compartment(
    patch$susceptible, exposure_rate, dt
  )
  
  patch$susceptible <- patch$susceptible -
    newly_exposed -
    deaths(patch$susceptible, patch$death_rate, dt) +
    births(patch$susceptible, patch$birth_rate, dt)
  
  newly_infected <-  to_next_compartment(
    patch$exposed, patch$infection_rate, dt
  )
  
  newly_infected_presymptomatic <- round(patch$prop_symptomatic * newly_infected)
  newly_infected_asymptomatic <- newly_infected - newly_infected_presymptomatic
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic, patch$recovery_rate_asym, dt
  )
  
  patch$infected_asymptomatic <- patch$infected_asymptomatic -
    newly_recovered_asymptomatic -
    deaths(patch$infected_asymptomatic, patch$death_rate, dt) +
    newly_infected_asymptomatic
  
  newly_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic <- patch$infected_presymptomatic -
    newly_symptomatic -
    deaths(patch$infected_presymptomatic, patch$death_rate, dt) +
    newly_infected_presymptomatic
  
  newly_recovered_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic <- patch$infected_symptomatic -
    newly_recovered_symptomatic -
    deaths(patch$infected_symptomatic, patch$death_rate, dt) +
    newly_symptomatic
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered_asymptomatic +
    newly_recovered_symptomatic
  
  patch
}

update_ksa_patch_screening <- function(patch, dt) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  newly_infected_diagnosed <-  to_next_compartment(
    patch$exposed_diagnosed, patch$infection_rate, dt
  )
  
  newly_infected_diag_presymptomatic <- round(patch$prop_symptomatic * newly_infected_diagnosed)
  newly_infected_diag_asymptomatic <- newly_infected_diagnosed - newly_infected_diag_presymptomatic
  
  patch$exposed_diagnosed <- patch$exposed_diagnosed -
    newly_infected_diagnosed -
    deaths(patch$exposed_diagnosed, patch$death_rate, dt)
  
  newly_recovered_diag_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic_diagnosed, patch$recovery_rate_asym, dt
  ) 
  
  patch$infected_asymptomatic_diagnosed <- patch$infected_asymptomatic_diagnosed -
    newly_recovered_diag_asymptomatic -
    deaths(patch$infected_asymptomatic_diagnosed, patch$death_rate, dt) +
    newly_infected_diag_asymptomatic
  
  newly_diag_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic_diagnosed, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic_diagnosed <- patch$infected_presymptomatic_diagnosed -
    newly_diag_symptomatic -
    deaths(patch$infected_presymptomatic_diagnosed, patch$death_rate, dt) +
    newly_infected_diag_presymptomatic
  
  newly_recovered_diag_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic_diagnosed, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic_diagnosed <- patch$infected_symptomatic_diagnosed -
    newly_recovered_diag_symptomatic -
    deaths(patch$infected_symptomatic_diagnosed, patch$death_rate, dt) +
    newly_diag_symptomatic
  
  patch$recovered_diagnosed <- patch$recovered_diagnosed -
    deaths(patch$recovered_diagnosed, patch$death_rate, dt) +
    newly_recovered_diag_asymptomatic +
    newly_recovered_diag_symptomatic
  
  patch
}

rate_to_probability <- function(rate, dt) {
  
  exp(-rate * dt)
  
}

get_number_migrating <- function(state, dt, compartments, movement_type, relative_movement) {
  
  if (movement_type == "rate") {
    
    pmat <- 1 - rate_to_probability(state$movement_rate, dt)
    
  } else {
    
    pmat <- state$movement_rate
    
  }
  
  # *** Not currently using this variation of the code, but can be added back in
  
  ## Include any compartment effects on movement
  ## e.g. infected people might move less than others
  # compartment_moving <- purrr::imap(compartments, function(comp, index) {
  #
  #   out <- pmat * relative_movement[index]
  #   diag(out) <- 0
  #   diag(out) <- 1 - rowSums(out)
  #   diag(out)[diag(out) < 0] <- 0  # to correct floating-point errors that can occur
  #   out
  #
  # })
  #
  # names(compartment_moving) <- compartments
  
  ## For each compartment, get the number of people moving in and
  ## out of patches.
  # n_moving <- vector(
  #   mode = "list", length = length(compartments)
  # )
  # names(n_moving) <- compartments
  n_patches <- length(state[["patches"]])
  # browser()
  
  # Extract the "patches" sublist
  patches_list <- state[["patches"]]
  
  # Use sapply to sum the first four elements (i.e. compartments) from each sublist
  sum_compartments <- sapply(patches_list, function(sublist) sum(unlist(sublist[1:4])))
  
  # Use conditional statement here
  # This is dummy code: stops us getting an error for NA for proportions if sum_compartments == 0
  # Also need to set n_movers==0 below so that we do not move anybody else
  prop_s <- ifelse(sum_compartments == 0,
                   0.25,
                   sapply(state[["patches"]], '[[', "susceptible") / sum_compartments)
  
  prop_e <- ifelse(sum_compartments == 0,
                   0.25,
                   sapply(state[["patches"]], '[[', "exposed") / sum_compartments)
  prop_i <- ifelse(sum_compartments == 0,
                   0.25,
                   sapply(state[["patches"]], '[[', "infected") / sum_compartments)
  prop_r <- ifelse(sum_compartments == 0,
                   0.25,
                   sapply(state[["patches"]], '[[', "recovered") / sum_compartments)
  
  out <- array(0, dim = c(n_patches, n_patches, 4))
  
  for (i in seq_len(n_patches)) {
    
    prop_compartments <- c(prop_s[i], prop_e[i], prop_i[i], prop_r[i])
    
    # Get the proportion of people moving to each destination from origin i
    # This can allows us to distribute the movers accordingly when we have small
    # numbers of movers (not currently needed).
    prop_movers <- state$movement_rate[i,] / sum(state$movement_rate[i,])
    
    if(sum_compartments[i] < sum(state$movement_rate[i,])) {
      
      all_movers <- stats::rmultinom(n = 1,
                                     size = sum_compartments[i],
                                     prob = prop_movers)
      
    }
    
    for (j in seq_len(n_patches)) {
      
      # Some code to handle the occasions where we have fewer people in a patch than we expect to move
      # In this model, movement from an origin will only be to a single destination, so
      # distributing the last remaining movers between different j values is not essential. But I
      # have set the code up so that it should be able to handle other versions where we have
      # two destinations.
      
      if(sum_compartments[i] >= sum(state$movement_rate[i,])) {
        n_movers <- state$movement_rate[i,j]
      } else {
        n_movers <- all_movers[j]
      }
      
      if(n_movers != 0) {
        out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
        
        max_compartment_values <- c(state$patches[[i]]$susceptible,
                                    state$patches[[i]]$exposed,
                                    state$patches[[i]]$infected,
                                    state$patches[[i]]$recovered)
        
        # Re-run the multinomial draw if any of the movement numbers exceed people in that compartment
        # This is an imperfect solution, but the problem only occurs in the final movement stage, so does not have substantial effect on results
        
        n <- 1 # set counter for monitoring number of re-draws needed
        
        while(sum(out[i,j,] > max_compartment_values) > 0) {
          
          warning(paste0(
            "In draw ", n, " the randomly drawn movers exceeded the number
            of people in one of the compartments. Another draw was made."
          ))
          n <- n + 1 # counter to see how many times we re-draw
          out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
          
        }
        
       } else {
        out[i,j,] <- c(0,0,0,0)
      }
    }
  }
  
  # convert array to a list of 4 elements, 1 per compartment
  n_moving <- lapply(1:dim(out)[3], function(i) array(out[, , i], dim = dim(out)[1:2]))
  names(n_moving) <- compartments
  n_moving
  
}

# Modified function for computing movement numbers in SEIR model with symptom compartments

get_number_migrating_symptoms <- function(state, dt, compartments, movement_type, relative_movement) {
  
  if (movement_type == "rate") {
    
    pmat <- 1 - rate_to_probability(state$movement_rate, dt)
    
  } else {
    
    pmat <- state$movement_rate
    
  }
  
  n_patches <- length(state[["patches"]])
  
  n_compartments <- length(compartments)

    # Extract the "patches" sublist
  patches_list <- state[["patches"]]
  
  # Use sapply to sum the first six elements (i.e. compartments) from each sublist
  sum_compartments <- sapply(patches_list, function(sublist) sum(unlist(sublist[compartments])))
  
  # Use conditional statement here
  # This is dummy code: stops us getting an error for NA for proportions if sum_compartments == 0
  # Also need to set n_movers==0 below so that we do not move anybody else
  # To do in future: find a way to set these more dynamically
  
  numbers_in_compartments <- vector(mode = "list", length = n_compartments)
  names(numbers_in_compartments) <- compartments
  
  proportions_in_compartments <- vector(mode = "list", length = n_compartments)
  names(proportions_in_compartments) <- compartments
  
  for (compartment in compartments) {
    
    numbers_in_compartments[[compartment]] <- sapply(state[["patches"]], '[[', compartment)
    
    proportions_in_compartments[[compartment]] <- ifelse(sum_compartments == 0,
                                                         1/n_compartments,
                                                         sapply(state[["patches"]], '[[', compartment) / sum_compartments)
    
  }
  
  # prop_s <- ifelse(sum_compartments == 0,
  #                  1/n_compartments,
  #                  sapply(state[["patches"]], '[[', "susceptible") / sum_compartments)
  # prop_e <- ifelse(sum_compartments == 0,
  #                  1/n_compartments,
  #                  sapply(state[["patches"]], '[[', "exposed") / sum_compartments)
  # prop_i_a <- ifelse(sum_compartments == 0,
  #                    1/n_compartments,
  #                  sapply(state[["patches"]], '[[', "infected_asymptomatic") / sum_compartments)
  # prop_i_p <- ifelse(sum_compartments == 0,
  #                    1/n_compartments,
  #                    sapply(state[["patches"]], '[[', "infected_presymptomatic") / sum_compartments)
  # prop_i_s <- ifelse(sum_compartments == 0,
  #                    1/n_compartments,
  #                    sapply(state[["patches"]], '[[', "infected_symptomatic") / sum_compartments)
  # prop_r <- ifelse(sum_compartments == 0,
  #                  1/n_compartments,
  #                  sapply(state[["patches"]], '[[', "recovered") / sum_compartments)
  
  out <- array(0, dim = c(n_patches, n_patches, n_compartments))
  
  for (i in seq_len(n_patches)) {
    
    # Extract the numbers and compartment proportions for this patch
    patch_number_compartments <- sapply(numbers_in_compartments, `[`, i)
    patch_prop_compartments <- sapply(proportions_in_compartments, `[`, i)
    # prop_compartments <- c(prop_s[i], prop_e[i], prop_i_a[i],
    #                        prop_i_p[i], prop_i_s[i], prop_r[i])
    
    # Get the proportion of people moving to each destination from origin i
    # This can allows us to distribute the movers accordingly when we have small
    # numbers of movers (not currently needed).
    prop_movers <- state$movement_rate[i,] / sum(state$movement_rate[i,])
    
    if(sum_compartments[i] < sum(state$movement_rate[i,])) {
      
      all_movers <- stats::rmultinom(n = 1,
                                     size = sum_compartments[i],
                                     prob = prop_movers)
      
    }
    
    for (j in seq_len(n_patches)) {
      
      # Some code to handle the occasions where we have fewer people in a patch than we expect to move
      # In this model, movement from an origin will only be to a single destination, so
      # distributing the last remaining movers between different j values is not essential. But I
      # have set the code up so that it should be able to handle other versions where we have
      # two destinations.
      
      if(sum_compartments[i] >= sum(state$movement_rate[i,])) {
        # when there are enough people in a patch/sub-patch, we move the fixed number denoted in movement_rate
        n_movers <- state$movement_rate[i,j]
      } else {
        # when we do not have enough ppl in compartment we use the numbers defined above
        n_movers <- all_movers[j]
      }
      
      if(n_movers != 0) {
        
        # Numbers in individual compartments
        people_in_compartments <- patch_number_compartments
        
        #Total number of people
        total_people <- sum(patch_number_compartments)
        
        # Calculate probabilities for sampling from each category
        probabilities <- people_in_compartments / total_people
        
        # Initialize a vector to store the compartments of the sampled movers
        sampled_movers <- character(n_movers)
        
        # Perform sampling without replacement
        for (each_mover in 1:n_movers) {
          # Sample a category based on probabilities
          sampled_compartment <- sample(names(people_in_compartments), 1, prob = probabilities)
          
          # Check if there are people left in the sampled compartment
          while (people_in_compartments[sampled_compartment] == 0) {
            sampled_compartment <- sample(names(people_in_compartments), 1, prob = probabilities)
          }
          
          # Reduce the count of people in the sampled compartment
          people_in_compartments[sampled_compartment] <- people_in_compartments[sampled_compartment] - 1
          
          # Store the sampled category
          sampled_movers[each_mover] <- sampled_compartment
        }
        
        sampled_movers_factor <- factor(sampled_movers, levels = compartments)
        sampled_counts <- as.vector(table(sampled_movers_factor))
        
        out[i,j,] <- sampled_counts
        
        # browser()
        
        # out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
        # 
        # max_compartment_values <- c(state$patches[[i]]$susceptible,
        #                             state$patches[[i]]$exposed,
        #                             state$patches[[i]]$infected_a,
        #                             state$patches[[i]]$infected_p,
        #                             state$patches[[i]]$infected_s,
        #                             state$patches[[i]]$recovered)
        # 
        # # Re-run the multinomial draw if any of the movement numbers exceed people in that compartment
        # # This is an imperfect solution, but the problem only occurs in the final movement stage, so does not have substantial effect on results
        # 
        # n <- 1 # set counter for monitoring number of re-draws needed
        # 
        # while(sum(out[i,j,] > max_compartment_values) > 0) {
        #   
        #   warning(paste0(
        #     "In draw ", n, " the randomly drawn movers exceeded the number
        #     of people in one of the compartments. Another draw was made."
        #   ))
        #   n <- n + 1 # counter to see how many times we re-draw
        #   out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
        #   
        # }
        
      } else {
        # This is for times where n_movers is zero
        out[i,j,] <- rep(0, times = n_compartments)
      }
    }
  }
  
  # convert array to a list of n elements, 1 per compartment
  n_moving <- lapply(1:dim(out)[3], function(i) array(out[, , i], dim = dim(out)[1:2]))
  names(n_moving) <- compartments
  n_moving
  
}

get_number_migrating_symptoms_copy <- function(state, dt, compartments, movement_type, relative_movement) {
  
  if (movement_type == "rate") {
    
    pmat <- 1 - rate_to_probability(state$movement_rate, dt)
    
  } else {
    
    pmat <- state$movement_rate
    
  }
  
  # *** Not currently using this variation of the code, but can be added back in
  
  ## Include any compartment effects on movement
  ## e.g. infected people might move less than others
  # compartment_moving <- purrr::imap(compartments, function(comp, index) {
  #
  #   out <- pmat * relative_movement[index]
  #   diag(out) <- 0
  #   diag(out) <- 1 - rowSums(out)
  #   diag(out)[diag(out) < 0] <- 0  # to correct floating-point errors that can occur
  #   out
  #
  # })
  #
  # names(compartment_moving) <- compartments
  
  ## For each compartment, get the number of people moving in and
  ## out of patches.
  # n_moving <- vector(
  #   mode = "list", length = length(compartments)
  # )
  # names(n_moving) <- compartments
  n_patches <- length(state[["patches"]])
  
  n_compartments <- length(compartments)
  
  # Extract the "patches" sublist
  patches_list <- state[["patches"]]
  
  # Use sapply to sum the first six elements (i.e. compartments) from each sublist
  sum_compartments <- sapply(patches_list, function(sublist) sum(unlist(sublist[1:n_compartments])))
  
  # Use conditional statement here
  # This is dummy code: stops us getting an error for NA for proportions if sum_compartments == 0
  # Also need to set n_movers==0 below so that we do not move anybody else
  # To do in future: find a way to set these more dynamically
  prop_s <- ifelse(sum_compartments == 0,
                   1/n_compartments,
                   sapply(state[["patches"]], '[[', "susceptible") / sum_compartments)
  prop_e <- ifelse(sum_compartments == 0,
                   1/n_compartments,
                   sapply(state[["patches"]], '[[', "exposed") / sum_compartments)
  prop_i_a <- ifelse(sum_compartments == 0,
                     1/n_compartments,
                     sapply(state[["patches"]], '[[', "infected_asymptomatic") / sum_compartments)
  prop_i_p <- ifelse(sum_compartments == 0,
                     1/n_compartments,
                     sapply(state[["patches"]], '[[', "infected_presymptomatic") / sum_compartments)
  prop_i_s <- ifelse(sum_compartments == 0,
                     1/n_compartments,
                     sapply(state[["patches"]], '[[', "infected_symptomatic") / sum_compartments)
  prop_r <- ifelse(sum_compartments == 0,
                   1/n_compartments,
                   sapply(state[["patches"]], '[[', "recovered") / sum_compartments)
  
  out <- array(0, dim = c(n_patches, n_patches, n_compartments))
  
  for (i in seq_len(n_patches)) {
    
    prop_compartments <- c(prop_s[i], prop_e[i], prop_i_a[i],
                           prop_i_p[i], prop_i_s[i], prop_r[i])
    
    # Get the proportion of people moving to each destination from origin i
    # This can allows us to distribute the movers accordingly when we have small
    # numbers of movers (not currently needed).
    prop_movers <- state$movement_rate[i,] / sum(state$movement_rate[i,])
    
    if(sum_compartments[i] < sum(state$movement_rate[i,])) {
      
      all_movers <- stats::rmultinom(n = 1,
                                     size = sum_compartments[i],
                                     prob = prop_movers)
      
    }
    
    for (j in seq_len(n_patches)) {
      
      # Some code to handle the occasions where we have fewer people in a patch than we expect to move
      # In this model, movement from an origin will only be to a single destination, so
      # distributing the last remaining movers between different j values is not essential. But I
      # have set the code up so that it should be able to handle other versions where we have
      # two destinations.
      
      if(sum_compartments[i] >= sum(state$movement_rate[i,])) {
        n_movers <- state$movement_rate[i,j]
      } else {
        n_movers <- all_movers[j]
      }
      
      if(n_movers != 0) {
        out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
        
        max_compartment_values <- c(state$patches[[i]]$susceptible,
                                    state$patches[[i]]$exposed,
                                    state$patches[[i]]$infected_a,
                                    state$patches[[i]]$infected_p,
                                    state$patches[[i]]$infected_s,
                                    state$patches[[i]]$recovered)
        
        # Re-run the multinomial draw if any of the movement numbers exceed people in that compartment
        # This is an imperfect solution, but the problem only occurs in the final movement stage, so does not have substantial effect on results
        
        n <- 1 # set counter for monitoring number of re-draws needed
        
        while(sum(out[i,j,] > max_compartment_values) > 0) {
          
          warning(paste0(
            "In draw ", n, " the randomly drawn movers exceeded the number
            of people in one of the compartments. Another draw was made."
          ))
          n <- n + 1 # counter to see how many times we re-draw
          out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
          
        }
        
      } else {
        out[i,j,] <- rep(0, times = n_compartments)
      }
    }
  }
  
  # convert array to a list of n elements, 1 per compartment
  n_moving <- lapply(1:dim(out)[3], function(i) array(out[, , i], dim = dim(out)[1:2]))
  names(n_moving) <- compartments
  n_moving
  
}

## n_moving is a matrix such that n_moving[i, j] is the
## number of people moving from i to j. Number of people
## moving out of patch i then sum(n_moving[i, j], i != j)

to_other_patches <- function(n_moving, patch_idx) {
  
  sum(n_moving[patch_idx, ]) - n_moving[patch_idx, patch_idx]
  
}

## n_moving is a matrix such that n_moving[i, j] is the
## number of people moving from i to j. Number of people
## moving into patch i sum(n_moving[, i], i != j)
from_other_patches <- function(n_moving, patch_idx) {
  
  sum(n_moving[, patch_idx]) - n_moving[patch_idx, patch_idx]
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
                             relative_movement = c(1, 1, 1, 1)
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating(state, dt, compartments, movement_type, relative_movement)
  n_patches <- length(state[["patches"]])
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many infections in KSA at this time, across all sub-patches
  ksa_infections <- ksa_patches %>% 
    map_dbl(~ .x$infected) %>% 
    reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                  pluck(.x, "infected"), pluck(.x, "recovered"))) %>%
    reduce(`+`)
  
  # What is ksa exposure rate, calculated across all sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * ksa_infections / ksa_total,
                              0)
  # browser()
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
                             relative_movement = c(1, 1, 1, 1, 1, 1)
) {
  
  movement_type <- match.arg(movement_type)
  
  n_moving <- get_number_migrating_symptoms(state, dt, compartments, movement_type, relative_movement)
  
  n_patches <- length(state[["patches"]])
  # browser()
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many symptomatic (incl. pre-symptomatic) and asymptomatic
  # infections in KSA at this time, across all sub-patches
  ksa_infections_symptomatic <- ksa_patches %>%
    map_dbl(~ sum(pluck(.x, "infected_presymptomatic"),
                  pluck(.x, "infected_symptomatic"))) %>%
    reduce(`+`)
  
  ksa_infections_asymptomatic <- ksa_patches %>% 
    map_dbl(~ .x$infected_asymptomatic) %>% 
    reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                  pluck(.x, "infected_asymptomatic"), pluck(.x, "infected_presymptomatic"),
                  pluck(.x, "infected_symptomatic"), pluck(.x, "recovered"))) %>%
    reduce(`+`)
  
  # What is ksa exposure rate, calculated across all sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  asymptomatic_infectivity <- ksa_patches[[1]]$asymptomatic_infectiousness
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * (
                                ksa_infections_symptomatic + 
                                  asymptomatic_infectivity * ksa_infections_asymptomatic
                              )  / ksa_total,
                              0)
  # browser()
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
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, ksa_exposure_rate)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt)
    }
    
  }
  state
}

# modify function to work with symptom compartments model and screening
update_ksa_state_symptoms_screening <- function(state,
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
                                      relative_movement = c(1, 1, 1, 1, 1, 1)
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
    apply(mat, c(1, 2), function(x) rbinom(1, x, test_rate))
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
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many symptomatic (incl. pre-symptomatic) and asymptomatic
  # infections in KSA at this time, across all sub-patches
  ksa_infections_symptomatic <- ksa_patches %>%
    map_dbl(~ sum(pluck(.x, "infected_presymptomatic"),
                  pluck(.x, "infected_symptomatic"))) %>%
    reduce(`+`)
  
  ksa_infections_asymptomatic <- ksa_patches %>% 
    map_dbl(~ .x$infected_asymptomatic) %>% 
    reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                  pluck(.x, "infected_asymptomatic"), pluck(.x, "infected_presymptomatic"),
                  pluck(.x, "infected_symptomatic"), pluck(.x, "recovered"))) %>%
    reduce(`+`)
  
  # What is ksa exposure rate, calculated across all ksa sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  asymptomatic_infectivity <- ksa_patches[[1]]$asymptomatic_infectiousness
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * (
                                ksa_infections_symptomatic + 
                                  asymptomatic_infectivity * ksa_infections_asymptomatic
                              )  / ksa_total,
                              0)
  
  for (idx in seq_len(n_patches)) {
    
    patch <- state[["patches"]][[idx]]

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
    
    if (idx %in% ksa_index) {
      # this first modified function uses the pre-specified exposure rate for KSA sub-patches
      # this was computed above
      state[["patches"]][[idx]] <- update_ksa_patch_symptoms(patch, dt, ksa_exposure_rate)
      # state[["patches"]][[idx]] <- update_ksa_patch_screening(patch, dt)
    } else {
      state[["patches"]][[idx]] <- update_patch_symptoms(patch, dt)
    }
    
  }
  state
}