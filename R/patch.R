##' Constructor for patch class
##'
##' Build a patch with seir compartments and patch-specific birth and death rates
##'
##' @param s_patch Number of suspectible individuals in patch
##' @param e_patch Number of exposed individuals in patch
##' @param i_patch Number of infected individuals in patch
##' @param r_patch Number of recovered individuals in patch
##' @param birth_rate Birth rate for this patch. Number of births per unit time
##' @param death_rate Death rate for this patch. Number of deaths per unit time in this patch
##' @return an object of class patch
##' @author Sangeeta Bhatia
make_patch <- function(s_patch,
                       e_patch,
                       i_patch,
                       r_patch,
                       birth_rate,
                       death_rate) {
    args <- list(
        s_patch, e_patch, i_patch, r_patch, birth_rate, death_rate
    )

    nonnumeric <- unlist(lapply(args, is.numeric))

    if (! any(nonnumeric)) {
        stop(
            "Error when trying to make a patch. At least one argument must be numeric",
            call. = FALSE
        )
    }


    out <- list(
        susceptible = s_patch,
        exposed = e_patch,
        infected = i_patch,
        recovered = r_patch,
        birth_rate = birth_rate,
        death_rate = death_rate
    )
    class(out) <- "patch"
    out
}
