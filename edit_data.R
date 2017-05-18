add_mass_column <-
    function(dataframe, rnames = rownames(dataframe)) {
        # Adds a column containing mz to the dataframe.
        dataframe = as.data.frame(dataframe)
        dataframe$m = str_replace(string = rnames,
                                  pattern = ".*_",
                                  replacement = "")
        return(dataframe)
    }

ParetoScale <- function(x, column = 2) {
    # Returns the pareto scaled version of the each column/row of the matrix.
    return(apply(
        X = x,
        MARGIN = column,
        FUN = function(j)
            (j - mean(j)) / sqrt(sd(j))
    ))
}

add_euclidean_distance <-
    function(dataframe, x = dataframe[, 1], y = dataframe[, 2]) {
        # Calculates Euclidean distance between x and y and adds it to a
        # dataframe column named 'ed'. If no value is given, assumes 1st and 2nd
        # columns are 'x' and 'y' respectively. Returns the entire dataframe.
        dataframe = as.data.frame(dataframe)
        dataframe$ed = sqrt(x ^ 2 + y ^ 2)
        return(dataframe)
    }
