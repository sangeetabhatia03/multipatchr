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
                  movement_rate
                  ) {

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

    if (! is.matrix(movement_rate) || any(movement_rate < 0)) {
        stop(
            "when trying to make a state.
             movement_rate should be a matrix of non-negative rates.",
            call. = FALSE
        )
    }

    n <- unlist(lapply(args, length))
    ## Check if all vectors have the same length
    ## If not, give warning a
    if (max(n) !=  min(n)) {
        warning(
            "Not all input vectors have the same length.
             Shorter vectors will be recycled."
        )
    }

    n_patches <- max(n)
    args <- lapply(args, rep, length.out = n_patches)

    state <- vector(
        mode = "list", length = 2
    )

    names(state) <- c("patches", "movement_rate")

    state[["patches"]] <- vector(
        mode = "list", length = n_patches
    )

    for (idx in seq_len(n_patches)) {

        patch_args <- lapply(args, '[[', idx)
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



##' Constructor for state class
##'
##' Build a state, which is a collection of
##' patches where each patch has SEIR
##' compartments, with the addition of asymptomatic,
##' pre-symptomatic and symptomatic compartments,
##' and patch-specific birth and
##' death rates.
##'
##' @param s_patches Vector of number of
##' susceptibles in each patch,
##' @param e_patches Vector of number of
##' exposed in each patch
##' @param i_a_patches Vector of number of
##' asymptomatic infected in each patch
##' @param i_p_patches Vector of number of
##' pre-symptomatic infected in each patch
##' @param i_s_patches Vector of number of
##' symptomatic infected in each patch
##' @param r_patches Vector of number of
##' recovered in each patch
##' @param birth_rates patch-specific birth rates
##' @param death_rates patch-specific death rates
##' @param transmission_rates Vector of patch-specific transmission
##' rate
##' @param asymptomatic_infectiousness Vector of proportions that 
##' give the infectivity of asymptomatic cases relative to
##' symptomatic cases.
##' @param presymptomatic_infectiousness Vector of proportions that 
##' give the infectivity of presymptomatic cases relative to
##' symptomatic cases.
##' @param prop_symptomatic Vector of patch-specific proportion
##' of cases that are symptomatic
##' @param infection_rates Vector of patch-specific infection rate
##' @param symptom_rates Vector of patch-specific rate at which
##' pre-symptomatic infections develop symptoms
##' @param recovery_rates_asym Vector of patch-specific recovery rate
##' for asymptomatic infections
##' @param recovery_rates_sym Vector of patch-specific recovery rate
##' for symptomatic infections
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
##' @author Jack Wardle
state_symptoms <- function(s_patches,
                  e_patches,
                  i_a_patches,
                  i_p_patches,
                  i_s_patches,
                  r_patches,
                  birth_rates,
                  death_rates,
                  transmission_rates,
                  asymptomatic_infectiousness,
                  presymptomatic_infectiousness,
                  prop_symptomatic,
                  infection_rates,
                  symptom_rates,
                  recovery_rates_asym,
                  recovery_rates_sym,
                  movement_rate
) {
  
  args <- list(
    s_patches = s_patches,
    e_patches = e_patches,
    i_a_patches = i_a_patches,
    i_p_patches = i_p_patches,
    i_s_patches = i_s_patches,
    r_patches = r_patches,
    birth_rates = birth_rates,
    death_rates = death_rates,
    transmission_rates = transmission_rates,
    asymptomatic_infectiousness = asymptomatic_infectiousness,
    presymptomatic_infectiousness = presymptomatic_infectiousness,
    prop_symptomatic = prop_symptomatic,
    infection_rates = infection_rates,
    symptom_rates = symptom_rates,
    recovery_rates_asym = recovery_rates_asym,
    recovery_rates_sym = recovery_rates_sym
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
  
  if (! is.matrix(movement_rate) || any(movement_rate < 0)) {
    stop(
      "when trying to make a state.
             movement_rate should be a matrix of non-negative rates.",
      call. = FALSE
    )
  }
  
  n <- unlist(lapply(args, length))
  ## Check if all vectors have the same length
  ## If not, give warning a
  if (max(n) !=  min(n)) {
    warning(
      "Not all input vectors have the same length.
             Shorter vectors will be recycled."
    )
  }
  
  n_patches <- max(n)
  args <- lapply(args, rep, length.out = n_patches)
  
  state <- vector(
    mode = "list", length = 2
  )
  
  names(state) <- c("patches", "movement_rate")
  
  state[["patches"]] <- vector(
    mode = "list", length = n_patches
  )
  
  for (idx in seq_len(n_patches)) {
    
    patch_args <- lapply(args, '[[', idx)
    state[["patches"]][[idx]] <- patch_symptoms(
      s_patch = patch_args$s_patches,
      e_patch = patch_args$e_patches,
      i_a_patch = patch_args$i_a_patches,
      i_p_patch = patch_args$i_p_patches,
      i_s_patch = patch_args$i_s_patches,
      r_patch = patch_args$r_patches,
      birth_rate = patch_args$birth_rates,
      death_rate = patch_args$death_rates,
      transmission_rate = patch_args$transmission_rates,
      asymptomatic_infectiousness = patch_args$asymptomatic_infectiousness,
      presymptomatic_infectiousness = patch_args$presymptomatic_infectiousness,
      prop_symptomatic = patch_args$prop_symptomatic,
      infection_rate = patch_args$infection_rates,
      symptom_rate = patch_args$symptom_rates,
      recovery_rate_asym = patch_args$recovery_rates_asym,
      recovery_rate_sym = patch_args$recovery_rates_sym
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


##' Constructor for state class
##'
##' Build a state, which is a collection of
##' patches where each patch has SEIR
##' compartments, with the addition of asymptomatic,
##' pre-symptomatic and symptomatic compartments,
##' and patch-specific birth and
##' death rates.
##'
##' @param s_patches Vector of number of
##' susceptibles in each patch,
##' @param e_patches Vector of number of
##' exposed in each patch
##' @param i_a_patches Vector of number of
##' asymptomatic infected in each patch
##' @param i_p_patches Vector of number of
##' pre-symptomatic infected in each patch
##' @param i_s_patches Vector of number of
##' symptomatic infected in each patch
##' @param e_diag_patches Vector of number of
##' diagnosed exposed in each patch
##' @param i_a_diag_patches Vector of number of
##' diagnosed asymptomatic infected in each patch
##' @param i_p_diag_patches Vector of number of
##' diagnosed pre-symptomatic infected in each patch
##' @param i_s_diag_patches Vector of number of
##' diagnosed symptomatic infected in each patch
##' @param r_patches Vector of number of
##' recovered in each patch
##' @param r_diag_patches Vector of number of
##' diagnosed recovered in each patch
##' @param s_false_patches Vector of number of 
##' susceptibles in patch with false positive 
##' test result.
##' @param r_false_patches Vector of number of 
##' recovered in patch with false positive 
##' test result.
##' @param birth_rates patch-specific birth rates
##' @param death_rates patch-specific death rates
##' @param transmission_rates Vector of patch-specific transmission
##' rate
##' @param asymptomatic_infectiousness Vector of proportions that 
##' give the infectivity of asymptomatic cases relative to
##' symptomatic cases.
##' @param presymptomatic_infectiousness Vector of proportions that 
##' give the infectivity of presymptomatic cases relative to
##' symptomatic cases.
##' @param prop_symptomatic Vector of patch-specific proportion
##' of cases that are symptomatic
##' @param infection_rates Vector of patch-specific infection rate
##' @param symptom_rates Vector of patch-specific rate at which
##' pre-symptomatic infections develop symptoms
##' @param recovery_rates_asym Vector of patch-specific recovery rate
##' for asymptomatic infections
##' @param recovery_rates_sym Vector of patch-specific recovery rate
##' for symptomatic infections
##' @param testing_rates Vector of patch-specific testing rates
##' @param test_sensitivity Vector of patch-specific test sensitivity
##' @param test_specificity Vector of patch-specific test specificity
##' @param isolation_period Vector of patch-specific mean duration of 
##' isolation following positive test
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
##' @author Jack Wardle
state_symptoms_testing <- function(s_patches,
                           e_patches,
                           i_a_patches,
                           i_p_patches,
                           i_s_patches,
                           r_patches,
                           e_diag_patches,
                           i_a_diag_patches,
                           i_p_diag_patches,
                           i_s_diag_patches,
                           r_diag_patches,
                           s_false_patches,
                           r_false_patches,
                           birth_rates,
                           death_rates,
                           transmission_rates,
                           asymptomatic_infectiousness,
                           presymptomatic_infectiousness,
                           prop_symptomatic,
                           infection_rates,
                           symptom_rates,
                           recovery_rates_asym,
                           recovery_rates_sym,
                           testing_rates,
                           test_sensitivity,
                           test_specificity,
                           isolation_period,
                           movement_rate
) {
  
  args <- list(
    s_patches = s_patches,
    e_patches = e_patches,
    i_a_patches = i_a_patches,
    i_p_patches = i_p_patches,
    i_s_patches = i_s_patches,
    r_patches = r_patches,
    e_diag_patches = e_diag_patches,
    i_a_diag_patches = i_a_diag_patches,
    i_p_diag_patches = i_p_diag_patches,
    i_s_diag_patches = i_s_diag_patches,
    r_diag_patches = r_diag_patches,
    s_false_patches = s_false_patches,
    r_false_patches = r_false_patches,
    birth_rates = birth_rates,
    death_rates = death_rates,
    transmission_rates = transmission_rates,
    asymptomatic_infectiousness = asymptomatic_infectiousness,
    presymptomatic_infectiousness = presymptomatic_infectiousness,
    prop_symptomatic = prop_symptomatic,
    infection_rates = infection_rates,
    symptom_rates = symptom_rates,
    recovery_rates_asym = recovery_rates_asym,
    recovery_rates_sym = recovery_rates_sym,
    testing_rates = testing_rates,
    test_sensitivity = test_sensitivity,
    test_specificity = test_specificity,
    isolation_period = isolation_period
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
  
  if (! is.matrix(movement_rate) || any(movement_rate < 0)) {
    stop(
      "when trying to make a state.
             movement_rate should be a matrix of non-negative rates.",
      call. = FALSE
    )
  }
  
  n <- unlist(lapply(args, length))
  ## Check if all vectors have the same length
  ## If not, give warning a
  if (max(n) !=  min(n)) {
    warning(
      "Not all input vectors have the same length.
             Shorter vectors will be recycled."
    )
  }
  
  n_patches <- max(n)
  args <- lapply(args, rep, length.out = n_patches)
  
  state <- vector(
    mode = "list", length = 2
  )
  
  names(state) <- c("patches", "movement_rate")
  
  state[["patches"]] <- vector(
    mode = "list", length = n_patches
  )
  
  for (idx in seq_len(n_patches)) {
    
    patch_args <- lapply(args, '[[', idx)
    state[["patches"]][[idx]] <- patch_symptoms_testing(
      s_patch = patch_args$s_patches,
      e_patch = patch_args$e_patches,
      i_a_patch = patch_args$i_a_patches,
      i_p_patch = patch_args$i_p_patches,
      i_s_patch = patch_args$i_s_patches,
      r_patch = patch_args$r_patches,
      e_diag_patch = patch_args$e_diag_patches,
      i_a_diag_patch = patch_args$i_a_diag_patches,
      i_p_diag_patch = patch_args$i_p_diag_patches,
      i_s_diag_patch = patch_args$i_s_diag_patches,
      r_diag_patch = patch_args$r_diag_patches,
      s_false_patch = patch_args$s_false_patches,
      r_false_patch = patch_args$r_false_patches,
      birth_rate = patch_args$birth_rates,
      death_rate = patch_args$death_rates,
      transmission_rate = patch_args$transmission_rates,
      asymptomatic_infectiousness = patch_args$asymptomatic_infectiousness,
      presymptomatic_infectiousness = patch_args$presymptomatic_infectiousness,
      prop_symptomatic = patch_args$prop_symptomatic,
      infection_rate = patch_args$infection_rates,
      symptom_rate = patch_args$symptom_rates,
      recovery_rate_asym = patch_args$recovery_rates_asym,
      recovery_rate_sym = patch_args$recovery_rates_sym,
      testing_rate = patch_args$testing_rates,
      test_sensitivity = patch_args$test_sensitivity,
      test_specificity = patch_args$test_specificity,
      isolation_period = patch_args$isolation_period
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



state_symptoms_testing_copy <- function(s_patches,
                                   e_patches,
                                   i_a_patches,
                                   i_p_patches,
                                   i_s_patches,
                                   r_patches,
                                   e_diag_patches,
                                   i_a_diag_patches,
                                   i_p_diag_patches,
                                   i_s_diag_patches,
                                   r_diag_patches,
                                   birth_rates,
                                   death_rates,
                                   transmission_rates,
                                   asymptomatic_infectiousness,
                                   presymptomatic_infectiousness,
                                   prop_symptomatic,
                                   infection_rates,
                                   symptom_rates,
                                   recovery_rates_asym,
                                   recovery_rates_sym,
                                   testing_rates,
                                   movement_rate
) {
  
  args <- list(
    s_patches = s_patches,
    e_patches = e_patches,
    i_a_patches = i_a_patches,
    i_p_patches = i_p_patches,
    i_s_patches = i_s_patches,
    r_patches = r_patches,
    e_diag_patches = e_diag_patches,
    i_a_diag_patches = i_a_diag_patches,
    i_p_diag_patches = i_p_diag_patches,
    i_s_diag_patches = i_s_diag_patches,
    r_diag_patches = r_diag_patches,
    birth_rates = birth_rates,
    death_rates = death_rates,
    transmission_rates = transmission_rates,
    asymptomatic_infectiousness = asymptomatic_infectiousness,
    presymptomatic_infectiousness = presymptomatic_infectiousness,
    prop_symptomatic = prop_symptomatic,
    infection_rates = infection_rates,
    symptom_rates = symptom_rates,
    recovery_rates_asym = recovery_rates_asym,
    recovery_rates_sym = recovery_rates_sym,
    testing_rates = testing_rates
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
  
  if (! is.matrix(movement_rate) || any(movement_rate < 0)) {
    stop(
      "when trying to make a state.
             movement_rate should be a matrix of non-negative rates.",
      call. = FALSE
    )
  }
  
  n <- unlist(lapply(args, length))
  ## Check if all vectors have the same length
  ## If not, give warning a
  if (max(n) !=  min(n)) {
    warning(
      "Not all input vectors have the same length.
             Shorter vectors will be recycled."
    )
  }
  
  n_patches <- max(n)
  args <- lapply(args, rep, length.out = n_patches)
  
  state <- vector(
    mode = "list", length = 2
  )
  
  names(state) <- c("patches", "movement_rate")
  
  state[["patches"]] <- vector(
    mode = "list", length = n_patches
  )
  
  for (idx in seq_len(n_patches)) {
    
    patch_args <- lapply(args, '[[', idx)
    state[["patches"]][[idx]] <- patch_symptoms_testing(
      s_patch = patch_args$s_patches,
      e_patch = patch_args$e_patches,
      i_a_patch = patch_args$i_a_patches,
      i_p_patch = patch_args$i_p_patches,
      i_s_patch = patch_args$i_s_patches,
      r_patch = patch_args$r_patches,
      e_diag_patch = patch_args$e_diag_patches,
      i_a_diag_patch = patch_args$i_a_diag_patches,
      i_p_diag_patch = patch_args$i_p_diag_patches,
      i_s_diag_patch = patch_args$i_s_diag_patches,
      r_diag_patch = patch_args$r_diag_patches,
      birth_rate = patch_args$birth_rates,
      death_rate = patch_args$death_rates,
      transmission_rate = patch_args$transmission_rates,
      asymptomatic_infectiousness = patch_args$asymptomatic_infectiousness,
      presymptomatic_infectiousness = patch_args$presymptomatic_infectiousness,
      prop_symptomatic = patch_args$prop_symptomatic,
      infection_rate = patch_args$infection_rates,
      symptom_rate = patch_args$symptom_rates,
      recovery_rate_asym = patch_args$recovery_rates_asym,
      recovery_rate_sym = patch_args$recovery_rates_sym,
      testing_rate = patch_args$testing_rates
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