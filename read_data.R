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

replicate_remover <- function(spectral_data, replicates) {
    # Removes the replicates by averaging them into one line (row).
    rows = nrow(spectral_data)
    cols = ncol(spectral_data)
    temp_matrix = matrix(
        data = rep(NA, rows * cols / replicates),
        nrow = rows / replicates,
        ncol = cols
    )
    for (i in 1:(rows / replicates)) {
        temp_matrix[i, ] =
            colMeans(spectral_data[(replicates * i - replicates + 1):(replicates * i),])
    }
    colnames(temp_matrix) = colnames(spectral_data)
    return(as.data.frame(temp_matrix))
}

get_spectral_data <- function(spectral_file, replicates = 1) {
    # Returns only the data without the first column since this column has the
    # rownames.
    spectral_data = as.data.frame(spectral_file[, 2:dim(spectral_file)[2]])
    if (replicates > 1) {
        spectral_data = replicate_remover(spectral_data, replicates)
    }
    return(spectral_data)
}

# The following functions could be combined. We'll see...
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


cleanup_rows <- function(dataframe,
                         dirty_rownames,
                         pattern = "[A-Z]+\\d+",
                         replicates = 1) {
    # Cleans up rownames. Allows for replicates.
    library(stringr)
    dataframe = as.data.frame(dataframe)
    if (replicates > 1) {
        repsequence = c(seq(from = 1,
                            to = length(dirty_rownames),
                            by = replicates))
        rownames(dataframe) = str_extract(
            string = dirty_rownames[repsequence],
            pattern = pattern)
    }
    else {
    rownames(dataframe) = str_extract(string = dirty_rownames,
                                      pattern = pattern)
    }
    return(dataframe)
}
