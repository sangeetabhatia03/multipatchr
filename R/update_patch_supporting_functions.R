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
    
    stop("Value moving to compartment cannot be NA", call. = FALSE)
  }
  
  out
  
}

rate_to_probability <- function(rate, dt) {
  
  exp(-rate * dt)
  
}

# TO DO: update this
get_number_migrating <- function(state, dt, compartments, movement_type, relative_movement) {
  
  if (movement_type == "rate") {
    
    pmat <- 1 - rate_to_probability(state$movement_rate, dt)
    
  } else {
    
    pmat <- state$movement_rate
    
  }
  
  n_patches <- length(state[["patches"]])
  
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
  
  # Use sapply to sum across all compartments in each patch
  # Returns a vector of the total people in each patch
  sum_compartments <- sapply(patches_list, function(sublist) sum(unlist(sublist[compartments])))
  
  # Use conditional statements here to set proportions in compartments
  # This is a bit of a fudge: stops us getting an error for NA for proportions if sum_compartments == 0
  # Also need to set n_movers==0 below so that we do not move anybody else
  
  numbers_in_compartments <- vector(mode = "list", length = n_compartments)
  names(numbers_in_compartments) <- compartments
  
  proportions_in_compartments <- vector(mode = "list", length = n_compartments)
  names(proportions_in_compartments) <- compartments
  
  for (compartment in compartments) {
    
    # Gives a list of length n_compartments; each element of the list is a vector
    # The elements of the vector give the number of ppl in that compartment for each patch
    numbers_in_compartments[[compartment]] <- sapply(state[["patches"]], '[[', compartment)
    
    # Gives a list of length n_compartments; each element of the list is a vector
    # The elements of the vector give the proportion of ppl in that compartment for each patch
    # When there are no people in the patch, we can assign a uniform fraction
    proportions_in_compartments[[compartment]] <- ifelse(sum_compartments == 0,
                                                         1/n_compartments,
                                                         sapply(state[["patches"]], '[[', compartment) / sum_compartments)
    
  }
  
  out <- array(0, dim = c(n_patches, n_patches, n_compartments))
  
  for (i in seq_len(n_patches)) {
    
    # Extract the numbers and compartment proportions for this patch
    patch_number_compartments <- sapply(numbers_in_compartments, `[`, i)
    patch_prop_compartments <- sapply(proportions_in_compartments, `[`, i)
    
    # Get the proportion of people moving to each destination from origin i
    # This can allows us to assign the detsinations of movers accordingly when we have small
    # numbers of movers (not currently needed as we only have single destinations -
    # depending on model phase these are KSA or home)
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
        # when we do not have enough ppl in compartment we use the all_movers numbers defined above
        n_movers <- all_movers[j]
      }
      
      if(n_movers != 0) {
        
        # Numbers in individual disease compartments
        people_in_compartments <- patch_number_compartments
        
        # Total number of people in this patch across all compartments
        total_people <- sum(patch_number_compartments)
        
        # Calculate probabilities for sampling from each category
        probabilities <- people_in_compartments / total_people
        
        # Draw from hypergeometric distribution to get the compartments of our movers
        # https://search.r-project.org/CRAN/refmans/extraDistr/html/MultiHypergeometric.html
        out[i,j,] <- extraDistr::rmvhyper(1, people_in_compartments, n_movers)
        
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