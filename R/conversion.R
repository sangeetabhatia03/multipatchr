##' Conversion of state object
##'
##' Convert state object to a data.frame
##'
##' @param state object of class `state`
##'
##' @param row.names names to be assigned to patches. If NULL,
##' names will be numeric from 1 to number of patches. Note that these
##' will actually be column values in the data.frame, rather than
##' row.names.
##' @param long if true, data frame is in long format with the columns
##' patch, variable, and value. Column patch contains the patch name
##' assigned through row.names (1:number of patches, if row.names
##' is NULL. the column variable are the elements of patch (
##' susceptible, exposed, infected, recovred, transmission_rate,
##' recovery_rate, birth_rate, death_rate
##' ))
##' @return data.frame
##' @author Sangeeta Bhatia
##' @export
as.data.frame.state <- function(state, row.names = NULL, ..., long = TRUE) {

    stopifnot(inherits(state, "state"))


    patches <- state[["patches"]]
    n_patches <- length(patches)

    if (is.null(row.names)) row.names <- seq_len(n_patches)

    names(patches) <- row.names

    out <- data.frame(
        patch = character(),
        variable = character(),
        value = numeric(),
        stringsAsFactors = FALSE
    )


    for (patch_idx in row.names) {

        patch <- patches[[patch_idx]]
        variables <- names(patch)
        for (var in variables) {
            df <- data.frame(
                patch = patch_idx,
                variable = var,
                value = patch[[var]]
            )

            out <- rbind(out, df)
        }
    }
    ## TODO: Affix movement_rate to this data frame.
    ## TODO: implement as.data.frame.patch and use it in the loop.
    ## TODO: Implement long = FALSE
   out

}
