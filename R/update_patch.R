update_patch <- function(patch, dt) {

    if (! inherits(patch, "patch")) {
        stop(
            "Error in updating patch. Argument not of class patch."
            call. = FALSE
        )
    }

    patch$suspectible <- patch$suspectible - to_exposed(patch$suspectible, exposure_rate, dt)
    patch$suspectible <- patch$suspectible + n_births(patch$birth_rate, dt)
    patch$suspectible <- patch$suspectible + to_other_patches(patch$suspectible, rates, dt)
    patch$suspectible <- patch$suspectible + from_other_patches(rates, dt)


}

## state is a collection of patches and a matrix of rates of movement
## between patches.
update_state <- function(state, dt) {

    for (i in seq_len(length(state))) {

        new_state[[i]] <- update_patch(patch, dt)

    }
    state
}
