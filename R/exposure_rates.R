# Function for computing the pilgrim exposure rate
compute_ksa_exposure_rate_original <- function(state, ksa_index) {
  
  # Find the sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # How many symptomatic, pre-symptomatic, and asymptomatic
  # infections in KSA at this time, across all sub-patches
  ksa_infections_symptomatic <- ksa_patches %>%
    purrr::map_dbl(~ .x$infected_symptomatic) %>% 
    purrr::reduce(`+`)
  
  ksa_infections_presymptomatic <- ksa_patches %>% 
    purrr::map_dbl(~ .x$infected_presymptomatic) %>% 
    purrr::reduce(`+`)
  
  ksa_infections_asymptomatic <- ksa_patches %>% 
    purrr::map_dbl(~ .x$infected_asymptomatic) %>% 
    purrr::reduce(`+`)
  
  # How many people in KSA at this time, across all sub-patches
  ksa_total <- ksa_patches %>%
    purrr::map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                         pluck(.x, "infected_asymptomatic"), pluck(.x, "infected_presymptomatic"),
                         pluck(.x, "infected_symptomatic"), pluck(.x, "recovered"))) %>%
    purrr::reduce(`+`)
  
  # What is ksa exposure rate, calculated across all sub-patches
  ksa_transmission_rate <- ksa_patches[[1]]$transmission_rate_pilgrims  # same for all so can extract form patch 1 only
  asymptomatic_infectivity <- ksa_patches[[1]]$asymptomatic_infectiousness
  presymptomatic_infectivity <- ksa_patches[[1]]$presymptomatic_infectiousness
  
  ksa_exposure_rate <- ifelse(ksa_total > 0,
                              ksa_transmission_rate * (
                                ksa_infections_symptomatic + 
                                  presymptomatic_infectivity * ksa_infections_presymptomatic +
                                  asymptomatic_infectivity * ksa_infections_asymptomatic
                              )  / ksa_total,
                              0)
  
}


compute_ksa_exposure_rate <- function(state, ksa_index, atrisk_index) {
  
  # Find the pilgrim sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # Get the at-risk non-pilgrim sub-patches in KSA using pre-defined atrisk_index
  atrisk_patches <- state[["patches"]][atrisk_index]
  
  # This returns I/N for i) pilgrim population and ii) at-risk non-pilgrims,
  # taking into account the relative infectiousness of symptomatic, presymptomatic,
  # and asymptomatic cases.
  
  pilgrim_relative_infections <- get_relative_number_infectious_cases(ksa_patches)
  atrisk_relative_infections <- get_relative_number_infectious_cases(atrisk_patches)
  
  mixing_population_size <- pilgrim_relative_infections$patch_group_total +
    atrisk_relative_infections$patch_group_total
  
  # Extract transmission rates between different groups
  pilgrim2pilgrim_transmission_rate <- ksa_patches[[1]]$transmission_rate_pilgrims  # same for all so can extract form patch 1 only
  atrisk2pilgrim_transmission_rate <- ksa_patches[[1]]$transmission_rate_pilgrims_and_atrisk
  
  # Calculate the components of the overall exposure rate
  pilgrim2pilgrim_exposure_rate <- ifelse(pilgrim_relative_infections$patch_group_total > 0,
                                          pilgrim2pilgrim_transmission_rate *
                                            pilgrim_relative_infections$relative_infectious_cases / 
                                            mixing_population_size,
                                          0)
  
  atrisk2pilgrim_exposure_rate <- ifelse(atrisk_relative_infections$patch_group_total > 0,
                                         atrisk2pilgrim_transmission_rate *
                                           atrisk_relative_infections$relative_infectious_cases / 
                                           mixing_population_size,
                                         0)
  
  # Sum to give the overall exposure rate
  overall_pilgrim_exposure_rate <- pilgrim2pilgrim_exposure_rate + atrisk2pilgrim_exposure_rate
  
}

# Function for computing the non pilgrim at risk rate
compute_atrisk_exposure_rate <- function(state, ksa_index, atrisk_index) {
  
  # Find the pilgrim sub-patches in KSA using pre-defined ksa_index
  ksa_patches <- state[["patches"]][ksa_index]
  
  # Get the at-risk non-pilgrim sub-patches in KSA using pre-defined atrisk_index
  atrisk_patches <- state[["patches"]][atrisk_index]
  
  # This returns I/N for i) pilgrim population and ii) at-risk non-pilgrims,
  # taking into account the relative infectiousness of symptomatic, presymptomatic,
  # and asymptomatic cases.
  
  pilgrim_relative_infections <- get_relative_number_infectious_cases(ksa_patches)
  atrisk_relative_infections <- get_relative_number_infectious_cases(atrisk_patches)
  
  mixing_population_size <- pilgrim_relative_infections$patch_group_total +
    atrisk_relative_infections$patch_group_total
  
  # Extract transmission rates between different groups
  atrisk2atrisk_transmission_rate <- ksa_patches[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  pilgrim2atrisk_transmission_rate <- ksa_patches[[1]]$transmission_rate_pilgrims_and_atrisk
  
  # Calculate the components of the overall exposure rate
  atrisk2atrisk_exposure_rate <- ifelse(atrisk_relative_infections$patch_group_total > 0,
                                        atrisk2atrisk_transmission_rate *
                                          atrisk_relative_infections$relative_infectious_cases / 
                                          mixing_population_size,
                                          0)
  
  pilgrim2atrisk_exposure_rate <- ifelse(pilgrim_relative_infections$patch_group_total > 0,
                                         pilgrim2atrisk_transmission_rate *
                                           pilgrim_relative_infections$relative_infectious_cases / 
                                           mixing_population_size,
                                         0)
  
  # Sum to give the overall exposure rate
  overall_atrisk_exposure_rate <- atrisk2atrisk_exposure_rate + pilgrim2atrisk_exposure_rate
  
}

get_relative_number_infectious_cases <- function(patch_group) {
  
  # How many symptomatic, pre-symptomatic, and asymptomatic
  # infections in patch group at this time, across all sub-patches of group
  infections_symptomatic <- patch_group %>%
    purrr::map_dbl(~ .x$infected_symptomatic) %>% 
    purrr::reduce(`+`)
  
  infections_presymptomatic <- patch_group %>% 
    purrr::map_dbl(~ .x$infected_presymptomatic) %>% 
    purrr::reduce(`+`)
  
  infections_asymptomatic <- patch_group %>% 
    purrr::map_dbl(~ .x$infected_asymptomatic) %>% 
    purrr::reduce(`+`)
  
  # How many ppl in patch group at this time
  patch_group_total <- patch_group %>%
    purrr::map_dbl(~ sum(pluck(.x, "susceptible"), pluck(.x, "exposed"),
                         pluck(.x, "infected_asymptomatic"), pluck(.x, "infected_presymptomatic"),
                         pluck(.x, "infected_symptomatic"), pluck(.x, "recovered"))) %>%
    purrr::reduce(`+`)
  
  
  # What is exposure rate, calculated across all sub-patches
  # transmission_rate <- patch_group[[1]]$transmission_rate  # same for all so can extract form patch 1 only
  asymptomatic_infectivity <- patch_group[[1]]$asymptomatic_infectiousness
  presymptomatic_infectivity <- patch_group[[1]]$presymptomatic_infectiousness
  
  relative_infectious_cases <- infections_symptomatic + 
    presymptomatic_infectivity * infections_presymptomatic +
    asymptomatic_infectivity * infections_asymptomatic
  
  out <- list(relative_infectious_cases = relative_infectious_cases,
              patch_group_total = patch_group_total)
  out
  
}


