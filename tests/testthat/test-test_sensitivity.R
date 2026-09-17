context("state-dependent test sensitivity")

test_that("validate_test_sensitivity recycles a single scalar to all states", {

  out <- validate_test_sensitivity(0.9)

  expect_equal(
    out,
    c(exposed = 0.9, infected_presymptomatic = 0.9,
      infected_asymptomatic = 0.9, infected_symptomatic = 0.9)
  )
})

test_that("validate_test_sensitivity accepts a complete named vector", {

  input <- c(infected_symptomatic = 0.9, exposed = 0.3,
            infected_presymptomatic = 0.6, infected_asymptomatic = 0.6)

  out <- validate_test_sensitivity(input)

  # order is standardised regardless of the order supplied
  expect_equal(
    out,
    c(exposed = 0.3, infected_presymptomatic = 0.6,
      infected_asymptomatic = 0.6, infected_symptomatic = 0.9)
  )
})

test_that("validate_test_sensitivity accepts a complete named list", {

  input <- list(exposed = 0.3, infected_presymptomatic = 0.6,
               infected_asymptomatic = 0.6, infected_symptomatic = 0.9)

  out <- validate_test_sensitivity(input)

  expect_equal(
    out,
    c(exposed = 0.3, infected_presymptomatic = 0.6,
      infected_asymptomatic = 0.6, infected_symptomatic = 0.9)
  )
})

test_that("validate_test_sensitivity errors on unrecognised names", {

  input <- c(exposed = 0.3, infected_presymptomatic = 0.6,
            infected_asymptomatic = 0.6, not_a_state = 0.9)

  expect_error(validate_test_sensitivity(input), "unrecognised names")
})

test_that("validate_test_sensitivity errors on missing states", {

  input <- c(exposed = 0.3, infected_presymptomatic = 0.6)

  expect_error(validate_test_sensitivity(input), "missing entries")
})

test_that("validate_test_sensitivity errors on out-of-range values", {

  expect_error(validate_test_sensitivity(1.1), "between 0 and 1")
  expect_error(validate_test_sensitivity(-0.1), "between 0 and 1")
})

test_that("patch_symptoms_testing stores a scalar test_sensitivity as a full named vector", {

  p <- patch_symptoms_testing(
    s_patch = 100, e_patch = 0, i_a_patch = 0, i_p_patch = 0, i_s_patch = 0,
    r_patch = 0, diag_patch = 0, s_false_patch = 0, r_false_patch = 0,
    birth_rate = 0, death_rate = 0, transmission_rate = 0.2,
    transmission_rate_pilgrims = 0.2, transmission_rate_pilgrims_and_atrisk = 0.2,
    asymptomatic_infectiousness = 0.5, presymptomatic_infectiousness = 1,
    prop_symptomatic = 0.6, infection_rate = 0.2, symptom_rate = 0.5,
    recovery_rate_asym = 0.1, recovery_rate_sym = 0.1, testing_rate = 0.5,
    test_sensitivity = 0.9, test_specificity = 0.9, isolation_period = 10
  )

  expect_equal(
    p$test_sensitivity,
    c(exposed = 0.9, infected_presymptomatic = 0.9,
      infected_asymptomatic = 0.9, infected_symptomatic = 0.9)
  )
})

test_that("state_symptoms_testing applies a state-dependent test_sensitivity to every patch", {

  state_dependent_sensitivity <- c(
    exposed = 0.3, infected_presymptomatic = 0.6,
    infected_asymptomatic = 0.6, infected_symptomatic = 0.9
  )

  s <- state_symptoms_testing(
    s_patches = c(100, 100), e_patches = c(0, 0), i_a_patches = c(0, 0),
    i_p_patches = c(0, 0), i_s_patches = c(0, 0), r_patches = c(0, 0),
    diag_patches = c(0, 0), s_false_patches = c(0, 0), r_false_patches = c(0, 0),
    birth_rates = c(0, 0), death_rates = c(0, 0), transmission_rates = c(0.2, 0.2),
    transmission_rates_pilgrims = c(0.2, 0.2),
    transmission_rates_pilgrims_and_atrisk = c(0.2, 0.2),
    asymptomatic_infectiousness = c(0.5, 0.5), presymptomatic_infectiousness = c(1, 1),
    prop_symptomatic = c(0.6, 0.6), infection_rates = c(0.2, 0.2), symptom_rates = c(0.5, 0.5),
    recovery_rates_asym = c(0.1, 0.1), recovery_rates_sym = c(0.1, 0.1), testing_rates = c(0.5, 0.5),
    test_sensitivity = state_dependent_sensitivity, test_specificity = c(0.9, 0.9),
    isolation_period = c(10, 10), movement_rate = matrix(c(1, 0, 0, 1), nrow = 2)
  )

  expect_equal(s$patches[[1]]$test_sensitivity, state_dependent_sensitivity)
  expect_equal(s$patches[[2]]$test_sensitivity, state_dependent_sensitivity)
})
