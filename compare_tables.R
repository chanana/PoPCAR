setwd("~/Lab/bucket_tables/")
source("~/Programming/R/PoPCAR/read_data.R")
source("~/Programming/R/PoPCAR/edit_data.R")
b4 = get_spectral_file("B499.txt")
b5 = get_spectral_file("B500.txt")
b45 = get_spectral_file("B499+B500.txt")

b4.d = get_spectral_data(b4)
b5.d = get_spectral_data(b5)
b45.d = get_spectral_data(b45)

b4.d = cleanup_cols(b4.d)
b5.d = cleanup_cols(b5.d)
b45.d = cleanup_cols(b45.d)

b4.d = cleanup_rows(b4.d, b4[, 1])
b5.d = cleanup_rows(b5.d, b5[, 1])
b45.d = cleanup_rows(b45.d, b45[, 1])

get_data_row <- function(parent_dataframe, name = NA, position = NA) {
    if (!is.na(name)) {
        return(parent_dataframe[name, ])
    }
    if (!is.na(position)) {
        return(parent_dataframe[position, ])
    }
    return(-1) # if neither is supplied, return -1
}

remove_zeroes <- function(dataframe, byrows = F) {
    # Removes zero columns or rows (as specified)
    if (!byrows) {
        # if the data is a row, change it to a column so it's faster
        dataframe = t(dataframe)
    }
    dataframe[which(dataframe == 0)] = NA
    dataframe = na.omit(dataframe)
    return(as.data.frame(dataframe))
}

get_merged_dataframe <- function(dataframe_a, dataframe_b,
                          name = NA, position_a = NA, position_b = NA) {
    # takes input of two dataframes 'a' and 'b' and merges them into one based
    # on the specified row either by name or by position. Assumes the rownmaes
    # are unique (key) and data is in a wide format.
    a = get_data_row(parent_dataframe = dataframe_a, name, position_a)
    b = get_data_row(parent_dataframe = dataframe_b, name, position_b)
    if (a == -1 || b == -1) {
        return(-1)
    }
    a = remove_zeroes(a)
    b = remove_zeroes(b)
    # will merge a and b and put things in a but not in b and vice versa
    # also into the merged dataframe. It assumes the key as the rownames and
    # that the dataframe is now in the long format.
    merged_dataframe = merge(
        x = a,
        y = b,
        by = 0,
        all.x = T,
        all.y = T,
        sort = F
    )
    return(as.data.frame(merged_dataframe))
}

# b4dup = t(b45.d[1, ])
# b5dup = t(b45.d[2, ])
#
# b4dup[which(b4dup == 0)] = NA
# b4dup = na.omit(b4dup)
#
# b4_4dup = merge(
#     x = as.data.frame(t(b4.d)),
#     y = as.data.frame(b4dup),
#     by = 0,
#     all.x = T,
#     all.y = T,
#     sort = F
# )

e = add_mass_column(t(b4.d), rnames = rownames(b4.d))

in_single_but_not_double = b4_4dup[is.na(b4_4dup$B499.y), ]
in_double_but_not_single = b4_4dup[is.na(b4_4dup$B499.x), ]

in_single_but_not_double$m = str_replace(
    string = in_single_but_not_double$Row.names,
    pattern = ".*_",
    replacement = ""
)
in_double_but_not_single$m = str_replace(
    string = in_double_but_not_single$Row.names,
    pattern = ".*_",
    replacement = ""
)

in_single_but_not_double = in_single_but_not_double[order(as.numeric(in_single_but_not_double$m), decreasing = T), ]
in_double_but_not_single = in_double_but_not_single[order(as.numeric(in_double_but_not_single$m), decreasing = T), ]

plot(as.numeric(in_single_but_not_double$m), pch = 16, ylab = "MZ", main = "'Different' buckets by mass")
points(
    as.numeric(in_double_but_not_single$m),
    col = "yellow",
    pch = 16,
    cex = 0.4
)

b4_4dup$m = str_replace(
    string = b4_4dup$Row.names,
    pattern = ".*_",
    replacement = ""
)
b4_4dup = b4_4dup[order(as.numeric(b4_4dup$m), decreasing = T),]
d = b4_4dup
d$B499.x[is.na(d$B499.x)] = 0
d$B499.y[is.na(d$B499.y)] = 0
d$B499.x[d$B499.x != 0] = 1
d$B499.y[d$B499.y != 0] = 2

e = as.numeric(in_double_but_not_single$m) - as.numeric(in_single_but_not_double$m)
png(filename = "OneVTwo-B499.png", width = 12, height = 9, units = "in", res = 300)
plot(
    x = e,
    # y = d$m,
    pch = 16,
    col = "blue",
    ylab = "Difference in MZ Value",
    xlab = "Index (Decreasing MZ >>>)",
    main = "Differences in MZ for the buckets found 'different'"
)
dev.off()
