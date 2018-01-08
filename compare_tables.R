setwd("~/Lab/bucket_tables/")
source("~/Programming/R/PoPCAR/read_data.R")
source("~/Programming/R/PoPCAR/edit_data.R")

df = get_spectral_file(filename = "BT_with_Adducts.txt")
d = get_spectral_data(spectral_file = df)
d = cleanup_cols(dataframe = d)
d = cleanup_rows(dataframe = d, dirty_rownames = df[, 1], pattern = "\\d+-[0-9]")

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

get_merged_dataframe <- function(dataframe_a,
                                 dataframe_b,
                                 dataframe_c,
                                 name = NA,
                                 position_a = NA,
                                 position_b = NA,
                                 position_c = NA) {
    # takes input of two dataframes 'a' and 'b' and merges them into one based
    # on the specified row either by name or by position. Assumes the rownmaes
    # are unique (key) and data is in a wide format.
    a = get_data_row(parent_dataframe = dataframe_a, name, position_a)
    b = get_data_row(parent_dataframe = dataframe_b, name, position_b)
    c = get_data_row(parent_dataframe = dataframe_c, name, position_c)
    if (a == -1 || b == -1 || c == -1) {
        return(-1)
    }
    a = remove_zeroes(a)
    b = remove_zeroes(b)
    c = remove_zeroes(c)
    # Add mass column for acting as the key
    a = add_mass_column(a)
    b = add_mass_column(b)
    c = add_mass_column(c)
    # will merge a and b and put things in a but not in b and vice versa
    # also into the merged dataframe. It assumes the key as the rownames and
    # that the dataframe is now in the long format.
    merged_dataframe = Reduce(function(df1, df2) merge(x = df1,
                                                       y = df2,
                                                       by = "m",
                                                       all.x = T,
                                                       all.y = T,
                                                       sort = F),
                              list(a, b, c))
    return(as.data.frame(merged_dataframe))
}
d.merged = get_merged_dataframe(dataframe_a = d,
                                dataframe_b = d,
                                dataframe_c = d,
                                position_a = 1,
                                position_b = 2,
                                position_c = 3)

d.merged = d.merged[order(d.merged$m, decreasing = T), ]

visualize_venn <- function(merged_dataframe, draw = T) {
    md = merged_dataframe
    names = colnames(md)[2:4]
    colnames(md)[2:4] = c('a', 'b', 'c')
    md$a[md$a != 0] = 1
    md$b[md$b != 0] = 2
    md$c[md$c != 0] = 4
    md[is.na(md)] = 0
    md$sum = as.factor(apply(X = md[,2:4], MARGIN = 1, FUN = sum))
    areas = table(md$sum)
    areas.list = c(
        sum(areas[c(1, 3, 5, 7)]),
        sum(areas[c(3, 2, 6, 7)]),
        sum(areas[c(4, 5, 6, 7)]),
        sum(areas[c(3, 7)]),
        sum(areas[c(6, 7)]),
        sum(areas[c(5, 7)]),
        sum(areas[7])
    )
    # areas.list = areas.list / sum(areas.list)
    if (draw) {
        library(VennDiagram)
        draw.triple.venn(
            area1 = areas.list[1],
            area2 = areas.list[2],
            area3 = areas.list[3],
            n12 = areas.list[4],
            n23 = areas.list[5],
            n13 = areas.list[6],
            n123 = areas.list[7],
            category = names,
            fill = c("#FF0000", "#00FF00", "#0000FF")
            # alpha = rep(0.5, 3) #Looks better without these
        )
    }
}
visualize_venn(d.merged)

get_differences <- function(merged_dataframe){
    md = merged_dataframe
    names = colnames(md)[2:4]
    colnames(md)[2:4] = c('a', 'b', 'c')
    md$a[md$a != 0] = 1
    md$b[md$b != 0] = 2
    md$c[md$c != 0] = 4
    md[is.na(md)] = 0
    md$sum = as.factor(apply(X = md[,2:4], MARGIN = 1, FUN = sum))
    a = merged_dataframe[md$sum == 1, ]
    b = merged_dataframe[md$sum == 3, ]
    c = merged_dataframe[md$sum == 2, ]
    d = merged_dataframe[md$sum == 5, ]
    e = merged_dataframe[md$sum == 4, ]
    f = merged_dataframe[md$sum == 6, ]
    g = merged_dataframe[md$sum == 7, ]
    plot(a[,1], col = "red", cex = 0.5, pch = 15)
    points(b[,1], col = "blue", cex = 0.5, pch = 16)
    points(c[,1], col = "green", cex = 0.5, pch = 17)
    points(d[,1], col = "grey", cex = 0.5, pch = 18)
    points(e[,1], col = "orange", cex = 0.5, pch = 19)
    points(f[,1], col = "purple", cex = 0.5, pch = 20)
    points(g[,1], col = "black", cex = 0.5, pch = 1)
    legend("topright",
           legend = c("1", "1+2", "2", "1+3", "3", "2+3", "1+2+3"),
           pch = c(seq(15, 20), 1),
           col = c("red", "blue", "green", "grey", "orange", "purple", "black"))
    # get masses in 'x' and not 'y' as x_not_y
    a_not_b = md[is.na(md$b), ]
    b_not_a = md[is.na(md$a), ]
    # sort in decreasing order of mass
    a_not_b = a_not_b[order(a_not_b$m, decreasing = T), ]
    b_not_a = b_not_a[order(b_not_a$m, decreasing = T), ]
    # make a list of the two dataframes
    return(list('df' = list('a' = a_not_b, 'b' = b_not_a),
                'number' = list('a' = dim(a_not_b)[1], 'b' = dim(b_not_a)[1])))
}
d.diff = get_differences(d.merged)

visualize_differences <- function(merged_dataframe) {
    md = get_differences(merged_dataframe) #md is now a list
    # Figure out which of a and b are larger and set the xlim to that +25
    m = c(dim(md$df$a)[1], dim(md$df$b)[1])
    # larger = which(m == max(m))
    # smaller = which(m == min(m))
    plot(md$df$b$m,
         col = "red",
         pch = 16,
         ylab = "m/z",
         main = "'Different' buckets by decreasing mass"
         )
    points(
        md$df$a$m,
        col = "blue",
        pch = 16,
        cex = 0.75
    )
    legend(
        x = "topright",
        legend = c(md$number$b, md$number$a),
        col = c("red", "blue"),
        pch = 16
    )
}
visualize_differences(d.merged)

view_difference_between_masses <- function(merged_dataframe){
    md = get_differences(merged_dataframe)
}
# Can't really do the following unless the buckets are exactly the same with minute differences.
# e = as.numeric(in_double_but_not_single$m) - as.numeric(in_single_but_not_double$m)
# png(filename = "OneVTwo-B499.png", width = 12, height = 9, units = "in", res = 300)
# plot(
#     x = e,
#     # y = d$m,
#     pch = 16,
#     col = "blue",
#     ylab = "Difference in MZ Value",
#     xlab = "Index (Decreasing MZ >>>)",
#     main = "Differences in MZ for the buckets found 'different'"
# )
# dev.off()
