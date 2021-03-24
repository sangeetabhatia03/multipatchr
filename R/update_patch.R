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

    exposure_rate <- patch$transmission_rate *
        patch$infected / (patch$susceptible +
                          patch$exposed +
                          patch$infected +
                          patch$recovered)

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
  
  ## Include any compartment effects on movement
  ## e.g. infected people might move less than others
  compartment_moving <- purrr::imap(compartments, function(comp, index) {
    
    out <- movement * relative_movement[index]
    diag(out) <- 0
    diag(out) <- 1 - rowSums(out)
    out
    
  })
  
  names(compartment_moving) <- compartments
  
  ## For each compartment, get the number of people moving in and
  ## out of patches.
  n_moving <- vector(
    mode = "list", length = length(compartments)
  )
  names(n_moving) <- compartments
  n_patches <- length(state[["patches"]])
  for (compartment in compartments) {
    n_current <- sapply(state[["patches"]], '[[', compartment)
    out <- matrix(
      NA, ncol = n_patches, nrow = n_patches
    )
    
    for (idx in seq_len(n_patches)) {
      
      out[idx, ] <- stats::rmultinom(
        n = 1,
        size = n_current[idx],
        prob = compartment_moving[[compartment]][idx, ]
      )[,1]
      
    }
    n_moving[[compartment]] <- out
  }
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
##' to 10% the normal level by replacing that vector element with 0.1.  
##' @return state updated
##' @author Sangeeta Bhatia
##' @export
update_state <- function(state,
                         dt,
                         compartments = c("susceptible",
                                          "exposed",
                                          "infected",
                                          "recovered"),
                         movement_type = c("rate", "probability"),
                         relative_movement = c(1, 1, 1, 1)
                        ) {
  
    movement_type <- match.arg(movement_type)
    
    n_moving <- get_number_migrating(state, dt, compartments, movement_type)
    n_patches <- length(state[["patches"]])
    for (idx in seq_len(n_patches)) {
      
      patch <- state[["patches"]][[idx]]
      
      for (compartment in compartments) {
        patch[[compartment]] <- patch[[compartment]] -
          to_other_patches(n_moving[[compartment]],  idx) +
          from_other_patches(n_moving[[compartment]], idx)
        
      }
      
      state[["patches"]][[idx]] <- update_patch(patch, dt)
      
    }
    state
  }

