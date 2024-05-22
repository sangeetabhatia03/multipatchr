set_finished_isolators <- function(old_state,
                                   finished_isolating_numbers,
                                   idx) {
  
  if (is.null(old_state)) {
    
    finished_isolating <- 0
    
  } else {
    
    finished_isolating <- finished_isolating_numbers[idx]
    
  }
  
  finished_isolating
  
}