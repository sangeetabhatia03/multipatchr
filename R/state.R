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
##' @param transmission_rates Vector of patch-specific transmission
##' rate
##' @param infection_rates Vector of patch-specific infection rate
##' @param recovery_rates Vector of patch-specific recovery rate
##' @param movement_rate M rate of movement between
##' patches is defined the rate matrix M where
##' m[i, j] is the rate of movement from i to j
##' conditional on moving out of i per unit time.
##' Must be n X n matrix where is the number of
##' patches i.e, max of lengths of all input
##' vectors s_patches, e_patches, i_patches,
##' r_patches, birth_rates, and death_rates.
##' @return List with class `state`. Each element
##' of `state` is a list of class `patch`.
##' In addition, `state` contains a `movement_rate`
##' which is a non-negative matrix of rates of
##' movement between patches.
##' @seealso [patch()]
##' @export
##' @author Sangeeta Bhatia
state <- function(s_patches,
                  e_patches,
                  i_patches,
                  r_patches,
                  birth_rates,
                  death_rates,
                  transmission_rates,
                  infection_rates,
                  recovery_rates,
                  movement_rate) {
  args <- list(
    s_patches = s_patches,
    e_patches = e_patches,
    i_patches = i_patches,
    r_patches = r_patches,
    birth_rates = birth_rates,
    death_rates = death_rates,
    transmission_rates = transmission_rates,
    infection_rates = infection_rates,
    recovery_rates = recovery_rates
  )

  nonnumeric <- unlist(lapply(args, is.numeric))

  if (!any(nonnumeric)) {
    stop(
      "when trying to make a state.
             At least one argument must be numeric",
      call. = FALSE
    )
  }

  if (!is.matrix(movement_rate) || any(movement_rate < 0)) {
    stop(
      "when trying to make a state.
             movement_rate should be a matrix of non-negative rates.",
      call. = FALSE
    )
  }

  n <- unlist(lapply(args, length))
  ## Check if all vectors have the same length
  ## If not, give a warning
  if (max(n) != min(n)) {
    warning(
      "Not all input vectors have the same length.
             Shorter vectors will be recycled."
    )
  }

  n_patches <- max(n)
  args <- lapply(args, rep, length.out = n_patches)

  state <- vector(
    mode = "list", length = 3
  )

  names(state) <- c("patches", "movement_rate", "n_individuals")

  state[["patches"]] <- vector(
    mode = "list", length = n_patches
  )

  ## Number of individuals in each patch in each compartment
  ## at this time. Columns are the 4 compartments.
  n_individuals <- matrix(NA, nrow = n_patches, ncol = 4)
  for (idx in seq_len(n_patches)) {
    n_individuals[idx, 1] <- s_patches
    n_individuals[idx, 2] <- e_patches
    n_individuals[idx, 3] <- i_patches
    n_individuals[idx, 4] <- r_patches
  }

  for (idx in seq_len(n_patches)) {
    patch_args <- lapply(args, "[[", idx)
    state[["patches"]][[idx]] <- patch(
      s_patch = patch_args$s_patches,
      e_patch = patch_args$e_patches,
      i_patch = patch_args$i_patches,
      r_patch = patch_args$r_patches,
      birth_rate = patch_args$birth_rates,
      death_rate = patch_args$death_rates,
      transmission_rate = patch_args$transmission_rate,
      infection_rate = patch_args$infection_rate,
      recovery_rate = patch_args$recovery_rate
    )
  }

  if (
    nrow(movement_rate) != n_patches ||
      ncol(movement_rate) != n_patches
  ) {
    stop(
      "when trying to make a state.
             movement_rate should be a ", n_patches,
      " X ", n_patches,
      " matrix",
      call. = FALSE
    )
  }

  state[["movement_rate"]] <- movement_rate
  class(state) <- "state"

  state
}
