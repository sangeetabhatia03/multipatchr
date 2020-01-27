##' Constructor for patch class
##'
##' Build a patch with seir compartments and patch-specific birth and death rates
##'
##' @param s_patch Number of susceptible individuals in patch
##' @param e_patch Number of exposed individuals in patch
##' @param i_patch Number of infected individuals in patch
##' @param r_patch Number of recovered individuals in patch
##' @param birth_rate Birth rate for this patch. Number of births per unit time
##' @param death_rate Death rate for this patch. Number of deaths per unit time in this patch
##' @param transmission_rate effective contact rate
##' @param infection_rate rate of moving from compartment E to I
##' @param recovery_rate rate of moving from compartment I to R.
##' @return a list of class patch which has the following items
##' ** susceptible ** Number of susceptibles in this patch
##' ** exposed ** Number of exposed individuals in this patch
##' ** infected ** Number of infected individuals in this patch
##' ** recovered ** Number of recovered individuals in this patch
##' ** birth_rate ** Patch-specific birth rate
##' ** death_rate ** Patch-specific death rate
##' ** transmission_rate ** Effective contact rate for this patch
##' ** Infection rate ** Rate at which individuals move from E to I
##' ** Recovery rate ** Rate at which individuals move from  I to R.
##' @export
##' @author Sangeeta Bhatia
patch <- function(s_patch,
                  e_patch,
                  i_patch,
                  r_patch,
                  birth_rate,
                  death_rate,
                  transmission_rate,
                  infection_rate,
                  recovery_rate) {
    args <- list(
        s_patch, e_patch, i_patch, r_patch, birth_rate, death_rate,
        transmission_rate, infection_rate, recovery_rate
    )

    nonnumeric <- unlist(lapply(args, is.numeric))

    if (! any(nonnumeric)) {
        stop(
            "Error when trying to make a patch.
             At least one argument must be numeric",
            call. = FALSE
        )
    }

    ## Number of individuals in any compartment must be
    ## an integer, else rbinom will generate NAs.
    are_integers <- lapply(
        c(s_patch, e_patch, i_patch, r_patch),
        function(x) all.equal(x, as.ineteger(x))
    )

    are_positive <- c(s_patch, e_patch, i_patch, r_patch) >= 0

    stopifnot(all(are_positive))

    if (! all(are_integers)) {

        warning(
            "Number of individuals in each compartment should be a
             positive integer. Rounding up."
        )

        s_patch <- ceiling(s_patch)
        e_patch <- ceiling(e_patch)
        i_patch <- ceiling(i_patch)
        r_patch <- ceiling(r_patch)

    }

    out <- list(
        susceptible = s_patch,
        exposed = e_patch,
        infected = i_patch,
        recovered = r_patch,
        birth_rate = birth_rate,
        death_rate = death_rate,
        transmission_rate = transmission_rate,
        infection_rate = infection_rate,
        recovery_rate = recovery_rate
    )
    class(out) <- "patch"
    out
}
