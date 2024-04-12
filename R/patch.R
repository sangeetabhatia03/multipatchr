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
    are_integers <- sapply(
        c(s_patch, e_patch, i_patch, r_patch),
        function(x) all.equal(x, as.integer(x))
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


##' Constructor for patch class
##'
##' Build a patch with seir compartments, including (a/pre-)symptomatic compartments,
##' and patch-specific birth and death rates
##'
##' @param s_patch Number of susceptible individuals in patch
##' @param e_patch Number of exposed individuals in patch
##' @param i_a_patch Number of asymptomatic infected individuals in patch
##' @param i_p_patch Number of presymptomatic infected individuals in patch
##' @param i_s_patch Number of symptomatic infected individuals in patch
##' @param r_patch Number of recovered individuals in patch
##' @param birth_rate Birth rate for this patch. Number of births per unit time
##' @param death_rate Death rate for this patch. Number of deaths per unit time in this patch
##' @param transmission_rate effective contact rate
##' @param asymptomatic_infectiousness Proportion giving 
##' the infectivity of asymptomatic cases relative to
##' symptomatic cases.
##' @param presymptomatic_infectiousness Proportion giving 
##' the infectivity of presymptomatic cases relative to
##' symptomatic cases.
##' @param prop_symptomatic Proportion of cases that are symptomatic.
##' @param infection_rate rate of moving from compartment E to I
##' @param symptom_rate Rate of moving from I_presymptomatic to I_symptomatic.
##' @param recovery_rate_asym Rate of moving from I_asymtpomatic to recovered
##' @param recovery_rate_sym Rate of moving from I_symtpomatic to recovered
##' @return a list of class patch which has the following items
##' ** susceptible ** Number of susceptibles in this patch
##' ** exposed ** Number of exposed individuals in this patch
##' ** infected_asymptomatic ** Number of asymptomatic infected individuals in this patch
##' ** infected_presymptomatic ** Number of presymptomatic infected individuals in this patch
##' ** infected_symptomatic ** Number of symptomatic infected individuals in this patch
##' ** recovered ** Number of recovered individuals in this patch
##' ** birth_rate ** Patch-specific birth rate
##' ** death_rate ** Patch-specific death rate
##' ** transmission_rate ** Effective contact rate for this patch
##' ** asymptomatic_infectiousness ** Relative infectivity of asymptomatic infected individuals
##' ** presymptomatic_infectiousness ** Relative infectivity of presymptomatic infected individuals
##' ** prop_symptomatic ** Proportion of infections with symptoms
##' ** Infection rate ** Rate at which individuals move from E to I_a or I_p
##' ** Symptom rate ** Rate at which individuals move from I_p to I_s
##' ** Asymptomatic recovery rate ** Rate at which individuals move from  I_a to R.
##' ** Symptomatic recovery rate ** Rate at which individuals move from  I_s to R.
##' 
##' @export
##' @author Jack Wardle
patch_symptoms <- function(s_patch,
                  e_patch,
                  i_a_patch,
                  i_p_patch,
                  i_s_patch,
                  r_patch,
                  birth_rate,
                  death_rate,
                  transmission_rate,
                  asymptomatic_infectiousness,
                  presymptomatic_infectiousness,
                  prop_symptomatic,
                  infection_rate,
                  symptom_rate,
                  recovery_rate_asym,
                  recovery_rate_sym) {
  args <- list(
    s_patch, e_patch, i_a_patch, i_p_patch, i_s_patch, r_patch,
    birth_rate, death_rate, transmission_rate, asymptomatic_infectiousness,
    presymptomatic_infectiousness, prop_symptomatic, infection_rate, symptom_rate,
    recovery_rate_asym, recovery_rate_sym
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
  are_integers <- sapply(
    c(s_patch, e_patch, i_a_patch, i_p_patch, i_s_patch, r_patch),
    function(x) all.equal(x, as.integer(x))
  )
  
  are_positive <- c(s_patch, e_patch,
                    i_a_patch, i_p_patch, i_s_patch, r_patch) >= 0
  
  stopifnot(all(are_positive))
  
  if (! all(are_integers)) {
    
    warning(
      "Number of individuals in each compartment should be a
             positive integer. Rounding up."
    )
    
    s_patch <- ceiling(s_patch)
    e_patch <- ceiling(e_patch)
    i_a_patch <- ceiling(i_a_patch)
    i_p_patch <- ceiling(i_p_patch)
    i_s_patch <- ceiling(i_s_patch)
    r_patch <- ceiling(r_patch)
    
  }
  
  out <- list(
    susceptible = s_patch,
    exposed = e_patch,
    infected_asymptomatic = i_a_patch,
    infected_presymptomatic = i_p_patch,
    infected_symptomatic = i_s_patch,
    recovered = r_patch,
    birth_rate = birth_rate,
    death_rate = death_rate,
    transmission_rate = transmission_rate,
    asymptomatic_infectiousness = asymptomatic_infectiousness,
    presymptomatic_infectiousness = presymptomatic_infectiousness,
    prop_symptomatic = prop_symptomatic,
    infection_rate = infection_rate,
    symptom_rate = symptom_rate,
    recovery_rate_asym = recovery_rate_asym,
    recovery_rate_sym = recovery_rate_sym
  )
  class(out) <- "patch"
  out
}


##' Constructor for patch class
##'
##' Build a patch with seir compartments, including (a/pre-)symptomatic compartments,
##' compartments for diagnosed individuals, and patch-specific birth and death rates
##'
##' @param s_patch Number of susceptible individuals in patch
##' @param e_patch Number of exposed individuals in patch
##' @param i_a_patch Number of asymptomatic infected individuals in patch
##' @param i_p_patch Number of presymptomatic infected individuals in patch
##' @param i_s_patch Number of symptomatic infected individuals in patch
##' @param r_patch Number of recovered individuals in patch
##' @param e_diag_patch Number of diagnosed exposed individuals in patch
##' @param i_a_diag_patch Number of diagnosed asymptomatic infected individuals in patch
##' @param i_p_diag_patch Number of diagnosed presymptomatic infected individuals in patch
##' @param i_s_diag_patch Number of diagnosed symptomatic infected individuals in patch
##' @param r_diag_patch Number of diagnosed recovered individuals in patch
##' @param birth_rate Birth rate for this patch. Number of births per unit time
##' @param death_rate Death rate for this patch. Number of deaths per unit time in this patch
##' @param transmission_rate effective contact rate
##' @param asymptomatic_infectiousness Proportion giving 
##' the infectivity of asymptomatic cases relative to
##' symptomatic cases.
##' @param presymptomatic_infectiousness Proportion giving 
##' the infectivity of presymptomatic cases relative to
##' symptomatic cases.
##' @param prop_symptomatic Proportion of cases that are symptomatic.
##' @param infection_rate rate of moving from compartment E to I
##' @param symptom_rate Rate of moving from I_presymptomatic to I_symptomatic.
##' @param recovery_rate_asym Rate of moving from I_asymptomatic to recovered
##' @param recovery_rate_sym Rate of moving from I_symptomatic to recovered
##' @param testing_rate Proportion of people that are tested
##' @return a list of class patch which has the following items
##' ** susceptible ** Number of susceptibles in this patch
##' ** exposed ** Number of exposed individuals in this patch
##' ** exposed_diagnosed ** Number of diagnosed exposed individuals in this patch
##' ** infected_asymptomatic ** Number of asymptomatic infected individuals in this patch
##' ** infected_asymptomatic_diagnosed ** Number of diagnosed asymptomatic infected individuals in this patch
##' ** infected_presymptomatic ** Number of presymptomatic infected individuals in this patch
##' ** infected_presymptomatic_diagnosed ** Number of diagnosed presymptomatic infected individuals in this patch
##' ** infected_symptomatic ** Number of symptomatic infected individuals in this patch
##' ** infected_symptomatic_diagnosed ** Number of diagnosed symptomatic infected individuals in this patch
##' ** recovered ** Number of recovered individuals in this patch
##' ** recovered_diagnosed ** Number of diagnosed recovered individuals in this patch
##' ** birth_rate ** Patch-specific birth rate
##' ** death_rate ** Patch-specific death rate
##' ** transmission_rate ** Effective contact rate for this patch
##' ** asymptomatic_infectiousness ** Relative infectivity of asymptomatic infected individuals
##' ** presymptomatic_infectiousness ** Relative infectivity of presymptomatic infected individuals
##' ** prop_symptomatic ** Proportion of infections with symptoms
##' ** Infection rate ** Rate at which individuals move from E to I_a or I_p
##' ** Symptom rate ** Rate at which individuals move from I_p to I_s
##' ** Asymptomatic recovery rate ** Rate at which individuals move from  I_a to R.
##' ** Symptomatic recovery rate ** Rate at which individuals move from  I_s to R.
##' ** Testing rate ** Proportion of individuals that are tested.
##' 
##' @export
##' @author Jack Wardle
patch_symptoms_testing <- function(s_patch,
                                   e_patch,
                                   i_a_patch,
                                   i_p_patch,
                                   i_s_patch,
                                   r_patch,
                                   e_diag_patch,
                                   i_a_diag_patch,
                                   i_p_diag_patch,
                                   i_s_diag_patch,
                                   r_diag_patch,
                                   birth_rate,
                                   death_rate,
                                   transmission_rate,
                                   asymptomatic_infectiousness,
                                   presymptomatic_infectiousness,
                                   prop_symptomatic,
                                   infection_rate,
                                   symptom_rate,
                                   recovery_rate_asym,
                                   recovery_rate_sym,
                                   testing_rate) {
  args <- list(
    s_patch, e_patch, i_a_patch, i_p_patch, i_s_patch, r_patch,
    e_diag_patch, i_a_diag_patch, i_p_diag_patch, i_s_diag_patch, r_diag_patch,
    birth_rate, death_rate, transmission_rate,
    asymptomatic_infectiousness, presymptomatic_infectiousness,
    prop_symptomatic, infection_rate, symptom_rate,
    recovery_rate_asym, recovery_rate_sym, testing_rate
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
  are_integers <- sapply(
    c(s_patch, e_patch, i_a_patch, i_p_patch, i_s_patch, r_patch,
      e_diag_patch, i_a_diag_patch, i_p_diag_patch, i_s_diag_patch, r_diag_patch),
    function(x) all.equal(x, as.integer(x))
  )
  
  are_positive <- c(s_patch, e_patch,
                    i_a_patch, i_p_patch, i_s_patch, r_patch,
                    e_diag_patch, i_a_diag_patch, i_p_diag_patch, i_s_diag_patch,
                    r_diag_patch) >= 0
  
  stopifnot(all(are_positive))
  
  if (! all(are_integers)) {
    
    warning(
      "Number of individuals in each compartment should be a
             positive integer. Rounding up."
    )
    
    s_patch <- ceiling(s_patch)
    e_patch <- ceiling(e_patch)
    i_a_patch <- ceiling(i_a_patch)
    i_p_patch <- ceiling(i_p_patch)
    i_s_patch <- ceiling(i_s_patch)
    r_patch <- ceiling(r_patch)
    e_diag_patch <- ceiling(e_diag_patch)
    i_a_diag_patch <- ceiling(i_a_diag_patch)
    i_p_diag_patch <- ceiling(i_p_diag_patch)
    i_s_diag_patch <- ceiling(i_s_diag_patch)
    r_diag_patch <- ceiling(r_diag_patch)
    
  }
  
  out <- list(
    susceptible = s_patch,
    exposed = e_patch,
    exposed_diagnosed = e_diag_patch,
    infected_asymptomatic = i_a_patch,
    infected_asymptomatic_diagnosed = i_a_diag_patch,
    infected_presymptomatic = i_p_patch,
    infected_presymptomatic_diagnosed = i_p_diag_patch,
    infected_symptomatic = i_s_patch,
    infected_symptomatic_diagnosed = i_s_diag_patch,
    recovered = r_patch,
    recovered_diagnosed = r_diag_patch,
    birth_rate = birth_rate,
    death_rate = death_rate,
    transmission_rate = transmission_rate,
    asymptomatic_infectiousness = asymptomatic_infectiousness,
    presymptomatic_infectiousness = presymptomatic_infectiousness,
    prop_symptomatic = prop_symptomatic,
    infection_rate = infection_rate,
    symptom_rate = symptom_rate,
    recovery_rate_asym = recovery_rate_asym,
    recovery_rate_sym = recovery_rate_sym,
    testing_rate = testing_rate
  )
  class(out) <- "patch"
  out
}