deaths <- function(n, death_rate, dt) {

    rbinom(1, size = n, prob = death_rate * dt)
}

births <- function(n, birth_rate, dt) {

    rbinom(1, size = n, prob = birth_rate * dt)
}

to_next_compartment <- function(n_current, rate, dt) {

    prob <- 1 - exp(-rate * dt)
    rbinom(1, size = n_current, prob = prob)

}


to_other_patches <- function(n_current, rates, dt) {

    prob <- 1 - exp(-rates * dt)
    rmultinom(n = 1, size = n_current, prob = prob)
}

from_other_patches <- function() {
}




update_patch <- function(patch, dt) {

    if (! inherits(patch, "patch")) {
        stop(
            "Error in updating patch. Argument not of class patch.",
            call. = FALSE
        )
    }

    exposure_rate <- beta *
        patch$infected / (patch$susceptible +
                          patch$exposed +
                          patch$infected +
                          patch$recovered)

    newly_exposed <- to_next_compartment(
          patch$susceptible, exposure_rate, dt
      )

    patch$susceptible <- patch$susceptible -
     to_other_patches(
       patch$susceptible, rates, dt
     ) -
     newly_exposed -
     deaths(patch$birth_rate, dt) +
     births(patch$birth_rate, dt) +
     from_other_patches(rates, dt)

    newly_infected <-  to_next_compartment(
        patch$exposed, infection_rate, dt
    )
    patch$exposed <- patch$exposed -
        newly_infected
     to_other_patches(
       patch$exposed, rates, dt
     ) -
     deaths(patch$birth_rate, dt) +
     newly_exposed +
     from_other_patches(rates, dt)

    newly_recovered <-  to_next_compartment(
          patch$infected, recovery_rate, dt
    )
    patch$infected <- patch$infected -
        newly_recovered -
     to_other_patches(
       patch$infected, rates, dt
     ) -
     deaths(patch$birth_rate, dt) +
     newly_infected +
     from_other_patches(rates, dt)


    patch$recovered <- patch$recovered -
     to_other_patches(
       patch$infected, rates, dt
     ) -
        deaths(patch$birth_rate, dt) +
        newly_recovered +
    from_other_patches(rates, dt)

    patch



}

## state is a collection of patches and a matrix of rates of movement
## between patches.
update_state <- function(state, dt) {

    pmat <- rate_to_probability(state$rate_matrix)
    comaprtments <- c("susceptible", "exposed", "infected", "recovered")
    ## For each compartment, get the number of people moving in and
    ## out of patches.
    n_moving <- vector(mode = "list", length = length(compartments))
    names(n_moving) <- compartments

    n_patches <- length(state[["patches"]])
    for (compartment in comaprtments) {
        n_current <- sapply(state[["patches"]], '[[', compartment)
        out <- matrix(
            NA, ncol = n_patches, nrow = n_patches
        )

        for (idx in seq_len(n_patches)) {

            out[idx, ] <- rmultinom(
                n = 1,
                size = n_current[idx],
                prob = pmat[idx, ]
            )[,1]

        }
        n_moving[[compartment]] <- out
    }

    for (i in seq_len(length(state))) {

        state[[i]] <- update_patch(patch, dt)

    }
    state
}
