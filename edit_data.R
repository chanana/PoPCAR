add_mass_column <- function(dataframe, rnames = rownames(dataframe)) {
    # Adds a column containing mz to the dataframe.
    dataframe = as.data.frame(dataframe)
    dataframe$m = str_replace(
        string = rnames,
        pattern = ".*_",
        replacement = ""
    )
    return(dataframe)
}
