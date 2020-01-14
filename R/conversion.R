##' Test for class state
##'
##' @param x object
##' @param ... optional arguments not used
##' @return TRUE of x if of class state, FALSE
##' otherwise.
##' @author Sangeeta Bhatia
##' @export
is.state <- function(x, ...) {

    inherits(x, "state")
}

##' Test for class patch
##'
##'
##' @param x object
##' @param ... optional arguments not used
##' @return TRUE if x is of class patch, FALSE otherwise
##' @author Sangeeta Bhatia
##' @export
is.patch <- function(x, ...) {

    inherits(x, "patch")
}



##' Convert a patch object to a data frame
##'
##' ..
##' @param x object of class `patch`
##' @param row.names name of patch. Note that this will be a column in
##' the returned data.frame. If NULL, name is set to 1. If a vector of
##' length greater than 1 is given, only the first name is used.
##'
##' @param ... any other arguments. not used
##' @param long whether the data.frame should be long (1 row for each
##' variable in the object of class `patch`), or wide (a single row
##' where each variable in the object of class patch `patch` is a
##' column.)
##' @return data.frame
##' @author Sangeeta Bhatia
as.data.frame.patch <- function(x, row.names = NULL, ..., long = TRUE) {

    stopifnot(is.patch(x))

    if (is.null(row.names)) row.names <- 1

    out <- data.frame(
        patch = character(),
        variable = character(),
        value = numeric(),
        stringsAsFactors = FALSE
    )

    variables <- names(x)
    for (var in variables) {
        df <- data.frame(
            patch = row.names[1],
            variable = var,
            value = patch[[var]]
        )
        out <- rbind(out, df)
     }

    out
}

##' Conversion of state object
##'
##' Convert state object to a data.frame
##'
##' @param x object of class `state`
##'
##' @param row.names names to be assigned to patches. If NULL,
##' names will be numeric from 1 to number of patches. Note that these
##' will actually be column values in the data.frame, rather than
##' row.names.
##' @param ...  unused arguments
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
as.data.frame.state <- function(x, row.names = NULL, ..., long = TRUE) {

    stopifnot(inherits(x, "state"))


    patches <- x[["patches"]]
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
        df <- as.data.frame(patch, row.names = patch_idx)
        out <- rbind(out, df)
    }
    ## TODO: Affix movement_rate to this data frame.
    ## TODO: Implement long = FALSE
   out

}
