##' Constructor for state class
##'
##' Build a state, which is a collection of
##' patches where each patch has SEIR
##' compartments and patch-specific birth and
##' death rates.
##'
##' @param s_patches Vector of number of
##' susceptibles in each patch,
##' @param e_patches Vector of number of
##' exposed in each patch
##' @param i_patches Vector of number of
##' infected in each patch
##' @param r_patches Vector of number of
##' recovered in each patch
##' @param birth_rates patch-specific birth rates
##' @param death_rates patch-specific death rates
##' @param rate_matrix M rate of movement between
##' patches is defined the rate matrix M where
##' m[i, j] is the rate of movement from i to j
##' conditional on moving out of i per unit time.
##' Must be n X n matrix where is the number of
##' patches i.e, max of lengths of all input
##' vectors s_patches, e_patches, i_patches,
##' r_patches, birth_rates, and death_rates.
##' @return List with class `state`. Each element
##' of `state` is a list of class `patch`.
##' In addition, `state` contains a `rate_matrix`
##' which is a non-negative matrix of rates of
##' movement bwteen patches.
##' @seealso [make_patch()]
##' @export
##' @author Sangeeta Bhatia
make_state <- function(s_patches,
                       e_patches,
                       i_patches,
                       r_patches,
                       birth_rates,
                       death_rates,
                       rate_matrix) {
    args <- list(
        s_patches = s_patches,
        e_patches = e_patches,
        i_patches = i_patches,
        r_patches = r_patches,
        birth_rates = birth_rates,
        death_rates = death_rates
    )

    nonnumeric <- unlist(
        lapply(args, is.numeric)
    )

    if (! any(nonnumeric)) {
        stop(
            "when trying to make a state.
             At least one argument must be numeric",
            call. = FALSE
        )
    }

    if (! is.matrix(rate_matrix) & any(rate_matrix < 0)) {
        stop(
            "when trying to make a state.
             rate_matrix should be a matrix of non-negative rates.",
            call. = FALSE
        )
    }

    n <- unlist(lapply(args, length))
    ## Check if all vectors have the same length
    ## If not, give warning a
    if (max(n) !=  min(n)) {
        warning(
            "Not all input vectors have the same length. Shorter vectors will be recycled."
        )
    }

    npatches <- max(n)
    args <- lapply(args, rep, length.out = npatches)

    state <- vector(
        mode = "list", length = 2
    )

    state[["patches"]] <- vector(
        mode = "list", length = npatches
    )

    for (idx in seq_len(npatches)) {

        patch_args <- lapply(args, '[[', idx)
        state[["patches"]][[idx]] <- make_patch(
            s_patch = patch_args$s_patches,
            e_patch = patch_args$e_patches,
            i_patch = patch_args$i_patches,
            r_patch = patch_args$r_patches,
            birth_rate = patch_args$birth_rates,
            death_rate = patch_args$death_rates
        )

    }

    if (
        nrow(rate_matrix) != npatches ||
        ncol(rate_matrix) != npatches
    ) {
        stop(
            "when trying to make a state.
             rate_matrix should be a ", npatches,
            " X ", npatches,
           " matrix",
            call. = FALSE
        )
    }

    state[["rate_matrix"]] <- rate_matrix
    class(state) <- "state"

    state

}
