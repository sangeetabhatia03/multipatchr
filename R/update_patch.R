deaths <- function(n, death_rate, dt) {
  
  stats::rbinom(1, size = n, prob = death_rate * dt)
}

births <- function(n, birth_rate, dt) {
  
  stats::rbinom(1, size = n, prob = birth_rate * dt)
}

to_next_compartment <- function(n_current, rate, dt) {
  
  prob <- 1 - rate_to_probability(rate, dt)
  
  stats::rbinom(1, size = n_current, prob = prob)
  
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
    # This allows us to distribute the movers accordingly when we have small
    # numbers of movers.
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
        n_movers <- all_movers[j] # TO DO: how can we edit this so that the final people are able to move?
        # we will need this in our model because we want all pilgrims to leave KSA after Hajj.]
      }
      
      if(n_movers != 0) {
        out[i,j,] <- stats::rmultinom(n = 1, size = n_movers, prob = prop_compartments)
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
