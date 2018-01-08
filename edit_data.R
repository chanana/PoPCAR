add_mass_column <-
    function(dataframe, rnames = rownames(dataframe)) {
        # Adds a column containing mz to the dataframe.
        dataframe = as.data.frame(dataframe)
        dataframe$m = str_replace(string = rnames,
                                  pattern = ".*_",
                                  replacement = "")
        dataframe$m = sapply(dataframe$m, as.numeric)
        return(dataframe)
    }
add_rt_column <-
    function(dataframe, rnames = rownames(dataframe)) {
        # Adds a column containing rt to the dataframe.
        dataframe = as.data.frame(dataframe)
        dataframe$rt = str_replace(string = rnames,
                                  pattern = "_.*",
                                  replacement = "")
        dataframe$rt = sapply(dataframe$rt, as.numeric)
        return(dataframe)
    }
pareto_scale <- function(matrix, column = 2) {
    # Returns the pareto scaled version of the each column/row of the matrix.
    return(apply(
        X = matrix,
        MARGIN = column,
        FUN = function(x)
            (x - mean(x)) / sqrt(sd(x))
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

fd.sort <- function(x) {
    sorted = sort(x, decreasing = T, index.return = T)
    return(sorted)
}
