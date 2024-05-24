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

# Function to allow mixing between all the sub-patches in ksa
update_ksa_patch <- function(patch, dt, patch_exposure_rate) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  exposure_rate <- patch_exposure_rate
  
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

# Modified functions for a model w/ symptom compartments
update_patch_symptoms <- function(patch, dt, screening = FALSE) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  # First apply transitions to diagnosed compartments if applicable
  if (screening) {
    patch <- update_patch_screening(patch, dt)
  }
  
  # Now apply transitions to undiagnosed compartments
  
  # conditional statement means that we can handle patches that become empty
  exposure_rate <- ifelse(patch$susceptible +
                            patch$exposed +
                            patch$infected_asymptomatic +
                            patch$infected_presymptomatic +
                            patch$infected_symptomatic +
                            patch$recovered > 0,
                          patch$transmission_rate * (
                            patch$infected_symptomatic +
                              (patch$presymptomatic_infectiousness * patch$infected_presymptomatic) + 
                              (patch$asymptomatic_infectiousness * patch$infected_asymptomatic)) /
                            (patch$susceptible +
                               patch$exposed +
                               patch$infected_asymptomatic +
                               patch$infected_presymptomatic +
                               patch$infected_symptomatic +
                               patch$recovered
                            ),
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
  
  newly_infected_presymptomatic <- round(patch$prop_symptomatic * newly_infected)
  newly_infected_asymptomatic <- newly_infected - newly_infected_presymptomatic
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic, patch$recovery_rate_asym, dt
  )
  
  patch$infected_asymptomatic <- patch$infected_asymptomatic -
    newly_recovered_asymptomatic -
    deaths(patch$infected_asymptomatic, patch$death_rate, dt) +
    newly_infected_asymptomatic
  
  newly_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic <- patch$infected_presymptomatic -
    newly_symptomatic -
    deaths(patch$infected_presymptomatic, patch$death_rate, dt) +
    newly_infected_presymptomatic
  
  newly_recovered_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic <- patch$infected_symptomatic -
    newly_recovered_symptomatic -
    deaths(patch$infected_symptomatic, patch$death_rate, dt) +
    newly_symptomatic
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered_asymptomatic +
    newly_recovered_symptomatic
  
  patch
}


update_ksa_patch_symptoms <- function(patch, dt, patch_exposure_rate,
                                      finished_isolating_s,
                                      finished_isolating_r,
                                      finished_isolating_infected,
                                      old_state,
                                      screening = FALSE) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  # First apply transitions to diagnosed compartments if applicable
  # if (screening) { # masking this part. We now use a numerical approximation to determine end isolation states
  # patch <- update_patch_screening(patch, dt)
  # }
  
  # Now apply transitions to undiagnosed compartments
  
  exposure_rate <- patch_exposure_rate
  
  newly_exposed <- to_next_compartment(
    patch$susceptible, exposure_rate, dt
  )
  
  if (patch$isolation_period == 0) {
    
    warning(
      "Isolation period is set at 0 days. All false positive susceptibles and
      recovereds will be returned to S and R."
    )
  }
  
  # newly_released_s_false <- to_next_compartment(
  #   patch$susceptible_false_positive, 1/patch$isolation_period, dt
  #   )
  
  newly_released_s_false <- finished_isolating_s

  patch$susceptible <- patch$susceptible -
    newly_exposed -
    deaths(patch$susceptible, patch$death_rate, dt) +
    births(patch$susceptible, patch$birth_rate, dt) +
    newly_released_s_false
  
  patch$susceptible_false_positive <- patch$susceptible_false_positive -
    newly_released_s_false
  
  newly_infected <-  to_next_compartment(
    patch$exposed, patch$infection_rate, dt
  )
  
  newly_infected_presymptomatic <- stats::rbinom(1, newly_infected, patch$prop_symptomatic)
  newly_infected_asymptomatic <- newly_infected - newly_infected_presymptomatic
  
  patch$exposed <- patch$exposed -
    newly_infected -
    deaths(patch$exposed, patch$death_rate, dt) +
    newly_exposed
  
  newly_recovered_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic, patch$recovery_rate_asym, dt
  )
  
  patch$infected_asymptomatic <- patch$infected_asymptomatic -
    newly_recovered_asymptomatic -
    deaths(patch$infected_asymptomatic, patch$death_rate, dt) +
    newly_infected_asymptomatic
  
  newly_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic <- patch$infected_presymptomatic -
    newly_symptomatic -
    deaths(patch$infected_presymptomatic, patch$death_rate, dt) +
    newly_infected_presymptomatic
  
  newly_recovered_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic <- patch$infected_symptomatic -
    newly_recovered_symptomatic -
    deaths(patch$infected_symptomatic, patch$death_rate, dt) +
    newly_symptomatic
  
  # newly_released_r_false <- to_next_compartment(
  #   patch$recovered_false_positive, 1/patch$isolation_period, dt
  # )
  newly_released_r_false <- finished_isolating_r
  
  patch$recovered <- patch$recovered -
    deaths(patch$recovered, patch$death_rate, dt) +
    newly_recovered_asymptomatic +
    newly_recovered_symptomatic +
    newly_released_r_false
  
  patch$recovered_false_positive <- patch$recovered_false_positive -
    newly_released_r_false
  
  
  if (!is.null(old_state)) {
   
    # Define the numbers leaving each of the isolation compartments
    newly_released_exposed <- finished_isolating_infected["exposed_diagnosed"]
    newly_released_infected_asymptomatic <- finished_isolating_infected["infected_asymptomatic_diagnosed"]
    newly_released_infected_presymptomatic <- finished_isolating_infected["infected_presymptomatic_diagnosed"]
    newly_released_infected_symptomatic <- finished_isolating_infected["infected_symptomatic_diagnosed"]
    newly_released_recovered <- finished_isolating_infected["recovered_diagnosed"]
    
    patch$exposed <- patch$exposed + newly_released_exposed
    patch$infected_asymptomatic <- patch$infected_asymptomatic + newly_released_infected_asymptomatic
    patch$infected_presymptomatic <- patch$infected_presymptomatic + newly_released_infected_presymptomatic
    patch$infected_symptomatic <- patch$infected_symptomatic + newly_released_infected_symptomatic
    patch$recovered <- patch$recovered + newly_released_recovered
    
  }
    
  patch
}

update_patch_screening <- function(patch, dt) {
  
  if (! inherits(patch, "patch")) {
    stop(
      "Error in updating patch. Argument not of class patch.",
      call. = FALSE
    )
  }
  
  newly_infected_diagnosed <-  to_next_compartment(
    patch$exposed_diagnosed, patch$infection_rate, dt
  )
  
  newly_infected_diag_presymptomatic <- stats::rbinom(1, newly_infected_diagnosed, patch$prop_symptomatic)
  newly_infected_diag_asymptomatic <- newly_infected_diagnosed - newly_infected_diag_presymptomatic
  
  patch$exposed_diagnosed <- patch$exposed_diagnosed -
    newly_infected_diagnosed -
    deaths(patch$exposed_diagnosed, patch$death_rate, dt)
  
  newly_recovered_diag_asymptomatic <-  to_next_compartment(
    patch$infected_asymptomatic_diagnosed, patch$recovery_rate_asym, dt
  ) 
  
  patch$infected_asymptomatic_diagnosed <- patch$infected_asymptomatic_diagnosed -
    newly_recovered_diag_asymptomatic -
    deaths(patch$infected_asymptomatic_diagnosed, patch$death_rate, dt) +
    newly_infected_diag_asymptomatic
  
  newly_diag_symptomatic <- to_next_compartment(
    patch$infected_presymptomatic_diagnosed, patch$symptom_rate, dt
  )
  
  patch$infected_presymptomatic_diagnosed <- patch$infected_presymptomatic_diagnosed -
    newly_diag_symptomatic -
    deaths(patch$infected_presymptomatic_diagnosed, patch$death_rate, dt) +
    newly_infected_diag_presymptomatic
  
  newly_recovered_diag_symptomatic <-  to_next_compartment(
    patch$infected_symptomatic_diagnosed, patch$recovery_rate_sym, dt
  )
  
  patch$infected_symptomatic_diagnosed <- patch$infected_symptomatic_diagnosed -
    newly_recovered_diag_symptomatic -
    deaths(patch$infected_symptomatic_diagnosed, patch$death_rate, dt) +
    newly_diag_symptomatic
  
  patch$recovered_diagnosed <- patch$recovered_diagnosed -
    deaths(patch$recovered_diagnosed, patch$death_rate, dt) +
    newly_recovered_diag_asymptomatic +
    newly_recovered_diag_symptomatic
  
  # newly_released_s_false <- to_next_compartment(
  #   patch$susceptible_false_positive, 1/patch$isolation_period, dt
  #   )
  # 
  # patch$susceptible
  
  patch
}