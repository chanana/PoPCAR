# Read data from spectral intensity table. Find row and column names and return
# the dataframe with proper row/column names.
get_spectral_file <- function(filename) {
    # Reads the bucket table file.
    library(data.table)
    spectral_file = as.data.frame(fread(
        input = filename,
        header = T,
        showProgress = T
    ))
    return(spectral_file)
}

get_spectral_data <- function(spectral_file) {
    # Returns only the data without the first column since this column ends up
    # containing the rownames.
    spectral_data = as.data.frame(spectral_file[, 2:dim(spectral_file)[2]])
    return(spectral_data)
}

# The following functions could be combined. We'll see...
cleanup_rows <- function(dataframe, dirty_rownames, pattern = "[A-Z]+\\d+") {
    # Assumes the default naming scheme is <letters><numbers>. It looks for that
    # combination in the first column of the dataframe produced by the function
    # 'get_spectral_data'. This must be done
    library(stringr)
    dataframe = as.data.frame(dataframe)
    rownames(dataframe) = str_extract(
        string = dirty_rownames,
        pattern = pattern
        )
    return(dataframe)
}

cleanup_cols <- function(dataframe, t = "min") {
    # Makes the column names of the form 'RT_M/Z'.
    library(stringr)
    if (t == "s") {
        pattern = "(\\d+.\\d+)s : (\\d+.\\d+)m\\/z"
    }
    colnames(dataframe) = str_replace(
        string = colnames(dataframe),
        pattern = pattern,
        replacement = "\\1_\\2"
    )
    return(dataframe)
}
