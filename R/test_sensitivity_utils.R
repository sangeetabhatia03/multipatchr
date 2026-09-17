# Helper to allow test sensitivity to be either a single value (applied
# uniformly, as before) or state-dependent (a named vector/list giving a
# distinct sensitivity for each infection state).

test_sensitivity_states <- c("exposed", "infected_presymptomatic",
                             "infected_asymptomatic", "infected_symptomatic")

validate_test_sensitivity <- function(test_sensitivity) {

  if (is.list(test_sensitivity)) {
    test_sensitivity <- unlist(test_sensitivity)
  }

  # A single unnamed value is recycled across all infection states,
  # preserving the previous (state-independent) behaviour.
  if (length(test_sensitivity) == 1 && is.null(names(test_sensitivity))) {
    test_sensitivity <- stats::setNames(
      rep(test_sensitivity, length(test_sensitivity_states)),
      test_sensitivity_states
    )
  }

  if (is.null(names(test_sensitivity)) || any(names(test_sensitivity) == "")) {
    stop(
      "test_sensitivity must be either a single number or a named vector/",
      "list with entries for: ", paste(test_sensitivity_states, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(names(test_sensitivity), test_sensitivity_states)
  if (length(unknown_names) > 0) {
    stop(
      "test_sensitivity contains unrecognised names: ",
      paste(unknown_names, collapse = ", "),
      ". Valid names are: ", paste(test_sensitivity_states, collapse = ", "),
      call. = FALSE
    )
  }

  missing_names <- setdiff(test_sensitivity_states, names(test_sensitivity))
  if (length(missing_names) > 0) {
    stop(
      "test_sensitivity is missing entries for: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(test_sensitivity < 0) || any(test_sensitivity > 1)) {
    stop("test_sensitivity values must be between 0 and 1.", call. = FALSE)
  }

  test_sensitivity[test_sensitivity_states]

}
