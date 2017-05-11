# Libraries
library(data.table)
library(stringr)
library(xlsx)
library(colorspace)
# Definitions
setwd("~/Lab/Bucket Tables/")
# bucket_table_name = "A133-B531.txt"
# bucket_table_name = "smaller-dataset-bucket.txt"
bucket_table_name = "Hits-All3Media.txt"
colors = rainbow_hcl(4)
FLAG = 1 # Set to 1 if you want to include Antibase
if (FLAG == 1) {
    ppm = 2 # Can be changed to allow user to input desired accuracy
    PROTON = 1.6726231 / 1.6605402 # mass in kg to mass in u
    ELECTRON = (9.1093897 / 1.6605402) / 10000
    SODIUM = 22.989768
    SODIUM_Plus = SODIUM - ELECTRON
}

# Functions
# General Functions ####
fd.sort <- function(x) {
    sorted = sort(x, decreasing = T, index.return = T)
    return(sorted)
}
ParetoScale <- function(x) {
    return(apply(
        X = x,
        MARGIN = 2,
        FUN = function(x)
            (x - mean(x)) / sqrt(sd(x))
    ))
}
add_euclidean_distance <- function(x) {
    x = cbind(x, sqrt((x[, 1]) ^ 2 + (x[, 2]) ^ 2))
    colnames(x)[3] = "Euclidean Distance"
    return(x)
}
ParetoScale <- function(x) {
    return(apply(
        X = x,
        MARGIN = 2,
        FUN = function(x)
            (x - mean(x)) / sqrt(sd(x))
    ))
}
filter_zero <- function(x) {
    x = x[which(x != 0)]
    return(x)
}
lookup_function_triplicate <- function(x) {
    a = rownames(f.d)[x[1]]
    b = rownames(f.d)[x[2]]
    c = rownames(f.d)[x[3]]
    if (a == b & b == c) {
        result = x[1]
    }
    else {
        result = -1
    }
    return(result)
}
lookup_function_duplicate <- function(x) {
    c = rownames(f.d)[x[1]]
    d = rownames(f.d)[x[2]]
    if (c == d) {
        result = x[1]
    }
    else {
        result = -1
    }
    return(result)
}
different_function <-
    function(x) {
        # if there are non-zero elements, this returns the indices of those elements in the list
        n = which(x != 0)
        # if there is only one non-zero element, return it.
        if (length(n) == 1) {
            result = n
        }
        # if there are exactly two values, check if they are duplicates.
        if (length(n) == 2) {
            result = lookup_function_duplicate(n)
        }
        # if there are exactly three values,check if they are triplicates.
        if (length(n) == 3) {
            result = lookup_function_triplicate(n)
        }
        # if there are more than three values, return -1.
        if (length(n) > 3) {
            result = -1
        }
        return(result)
    }

filter_M <-
    function(m, a, p) {
        which(abs((m - a) * 1000000 / a) <= p) #it is a hit
    }
filter_H <-
    function(m, a, p) {
        which(abs((m - PROTON - a) * 1000000 / a) <= p) #it is a hit
    }
filter_Na <-
    function(m, a, p) {
        which(abs((m - SODIUM_Plus - a) * 1000000 / a) <= p) #it is a hit
    }
DrawFigures <-
    function(filename, x) {
        y = fd.highestPCs[1:2, x]
        tiff(
            filename = filename,
            height = 6,
            width = 12,
            units = "in",
            res = 300
        )
        par(mfrow = c(1, 2))
        plot(fd.pca$x[, y],
             pch = 16,
             col = "blue",
             main = "Scores")
        points(x = fd.pca$x[x, y[1]],
               y = fd.pca$x[x, y[2]],
               pch = 16,
               col = "red")
        # text(
        #     fd.pca$x[x, y],
        #     col = "orange",
        #     labels = rownames(fd.pca$x)[x],
        #     pos = 1,
        #     font = 2
        # )
        plot(fd.pca$rotation[, y],
             pch = 16,
             col = "blue",
             main = "Loadings")
        points(zd[, 4:5],
               col = "red",
               pch = 16)
        dev.off()
    }

# Data ####
f = as.data.frame(fread(input = bucket_table_name, header = T))
f.d = f[, 2:length(f)]  # Take only the numerical values

### Add the row and column names
#Specific pattern depending on naming convention
# pattern = "[A-Z]\\d+" #For both datasets
pattern = "[A-Z]{4}\\d+_[A-Z, 0-9]*"
rownames(f.d) = str_extract(string = f[, 1], pattern = pattern) # For larger dataset
# rownames(f.d)[4:dim(f.d)[1]] = str_extract(string = f[4:dim(f)[1] ,1],
#                             pattern = pattern) # For smaller dataset
# rownames(f.d)[1:3] = c("B499", "A165", "A2171") # For smaller dataset
# Generic pattern for column names assuming data is in min not seconds
pattern = "(\\d+.\\d+)min : (\\d+.\\d+)m\\/z"
colnames(f.d) = str_replace(string = colnames(f)[2:length(f)],
                            pattern = pattern,
                            replacement = "\\1_\\2")
### Apply Pareto Scaling
f.d.scaled.pareto = as.data.frame(ParetoScale(f.d))
### Run the PCA (scaling and centering has been done above)
fd.pca = prcomp(x = f.d.scaled.pareto,
                scale. = F,
                center = F)
## Figure containing Loadings plot with color and graph showing how quickly Euclidean Distance decreases in a typical loadings plot.
d = fd.pca$rotation[, 1:2]
d = as.data.frame(add_euclidean_distance(d))
d = d[order(d$`Euclidean Distance`, decreasing = T), ]
d$col = rep(c("red", "blue", "green", "orange"),
            times = c(2000, 2000, 2000, dim(d)[1] - 6000))
tiff(
    filename = "Colored-Loadings-smaller-dataset.tiff",
    width = 12,
    height = 6,
    units = "in",
    res = 300
)
par(mfrow = c(1, 2))
plot(
    x = d$PC1,
    y = d$PC2,
    pch = 16,
    col = d$col,
    xlab = "PC1",
    ylab = "PC2",
    main = "Loadings Plot"
)
plot(
    x = d$`Euclidean Distance`,
    pch = 16,
    col = d$col,
    xlab = "Index",
    ylab = "Euclidean Distance",
    main = "Euclidean Distance vs Index"
)
dev.off()

### Highest PCs
fd.highestPCs = as.data.frame(apply(
    X = abs(fd.pca$x),
    MARGIN = 1,
    FUN = function(i)
        fd.sort(i)
))
fd.highestPCs = fd.highestPCs[, which(str_detect(string = colnames(fd.highestPCs),
                                                 pattern = ".ix"))]
colnames(fd.highestPCs) = str_replace(
    string = colnames(fd.highestPCs),
    pattern = ".ix",
    replacement = ""
)
## To get pareto scaled data back scores %*% t(loadings)
#invisible(fd.pca$x %*% t(fd.pca$rotation)) #invisible hides output

### Get unique masses
fd.unique = apply(
    X = f.d,
    MARGIN = 2,
    FUN = function(x)
        different_function(x)
)
## Fish out the -1s and return everything else's column index
filter_uniquemass = as.integer(which(fd.unique != -1))
# If we apply the Unique masses filter to the bucket table,
# and then from that,
# we filter out the columns which are zero,
# we can get the unique masses for THAT row (strain). Note:the zero filtering step happens later,
# for now we simply take the masses which are unique to a strain from  the bucket table.
buckets_unique = f.d[, filter_uniquemass]
# Antibase

if (FLAG == 1) {
    # Read antibase
    antibase <-
        read.csv("~/Lab/antibase_tableout.csv",
                 stringsAsFactors = FALSE)
    # Take exact mass column
    antibase.exactmass = antibase$StructCalc
    # Cleanup: finds the first space and deletes it and
    # everything after it leaving only the exact mass.
    antibase.exactmass = str_replace(string = antibase.exactmass,
                                     pattern = " .*",
                                     replacement = "")
    # Make it of numeric type. This causes
    # an 'NAs introduced by coercion' warning. This is ok because we want
    # all non numeric values to be NA and can convert them to 0 later
    antibase.exactmass = as.numeric(antibase.exactmass)
    antibase.exactmass[which(is.na(antibase.exactmass))] = 0
    # Fill in blanks in the 'Names' column
    antibase$Name[which(antibase$Name == "")] = "NO_NAME"
}
# mainDir = "~/Lab/bucket_tables/"
# subDir = "outputDirectory/"
# ifelse(!dir.exists(file.path(mainDir, subDir)),
#        dir.create(file.path(mainDir, subDir)),
#        FALSE)
wb = createWorkbook()
saveWorkbook(wb, paste(bucket_table_name, ".xlsx", sep = ""))
# Final Program
# Takes row number 'r' from the modified bucket table and outputs the unique masses ordered by their euclidean distance on the loadings plot for that strain. Also outputs list of AntiBase matches and non -
#     matches into an Excel spreadsheet based on user set FLAG at the beginning of the script.
# Note z,
# d are temporary variables for the loop.
for (r in 1:(dim(f.d)[1])) {
    d = fd.pca$rotation[, fd.highestPCs[1:2, r]] #table containing PC planes of the strain in question
    d = as.data.frame(add_euclidean_distance(d)) #add euclidean distance
    d$m = str_replace(string = rownames(d),
                      pattern = ".*_",
                      replacement = "") # add column of masses to serve as a key
    z = as.data.frame(t(buckets_unique[r, ])) # masses for that strain
    z$rt = str_replace(string = rownames(z),
                       pattern = "_.*",
                       replacement = "") # add a time column
    z$m = str_replace(string = rownames(z),
                      pattern = ".*_",
                      replacement = "") # Add mass column to serve as a key
    z = as.data.frame(z[z != 0, ]) # remove zeroes - only unique strain masses left
    zd = merge(
        x = z,
        y = d,
        by.x = 'm',
        by.y = 'm',
        sort = F
    ) # merge on keys
    zd = zd[order(zd$`Euclidean Distance`, decreasing = T), ] # sort by Euclidean distance
    colnames(zd)[2] = "Intensity"
    zd = zd[, c(3, 1, 2, 4, 5, 6)]
    if (dim(zd)[1] >= 100) {
        zd = zd[1:100, ]
    }
    ## DRAW FIGURE!
    DrawFigures(paste(rownames(fd.pca$x)[r], ".tiff", sep = ""), r)
    if (FLAG == 1) {
        # Now search for this mass list in antibase, and for each hit,
        # report the name of the hit. Pretty sure I don't need the next three lines.
        # mzrt = read.table(text = rownames(zd),
        #                   sep = "_",
        #                   colClasses = "numeric")
        # z$rt = mzrt$V1
        # z$mz = mzrt$V2
        # Column for putting in (later) which masses are NOT in AntiBase
        zd$Unique = NA
        # Search in AntiBase
        match_mass = lapply(as.numeric(zd$m), function(x)
            filter_M(x, antibase.exactmass, ppm))
        match_hydrogen = lapply(as.numeric(zd$m), function(x)
            filter_H(x, antibase.exactmass, ppm))
        match_sodium = lapply(as.numeric(zd$m), function(x)
            filter_Na(x, antibase.exactmass, ppm))
        rows_for_mass = sum(unlist(lapply(match_mass, function(x)
            length(x))))
        rows_for_hydrogen = sum(unlist(lapply(match_hydrogen, function(x)
            length(x))))
        rows_for_sodium = sum(unlist(lapply(match_sodium, function(x)
            length(x))))
        #Finds how many matches were found for all of the masses across M, H, Na
        totalrowsneeded = sum(rows_for_mass,
                              rows_for_sodium,
                              rows_for_hydrogen)

        # Next, we want to filter out the masses from the 'match_mass' which are hits in AntiBase. It is a list of lists, where the upper level is the index in the 'match_mass' and lower number is the index in the antibase.exactmass list.
        # Excel sheet preparations
        wb = loadWorkbook(file = paste(bucket_table_name, ".xlsx", sep = ""))
        # make sheet = strain name
        sheet = createSheet(wb,
                            sheetName = rownames(buckets_unique)[r])
        # make headers for the row
        firstRow = matrix(
            data = c("RT_MZ",
                     "M Match",
                     "M--H Match",
                     "M--Na Match"),
            nrow = 1
        )
        # define a cell block
        cb = CellBlock(
            sheet = sheet,
            startRow = 1,
            startColumn = 1,
            noRows = 1,
            noColumns = 4,
            create = T
        )
        # put the headers there
        CB.setMatrixData(
            cellBlock = cb,
            x = firstRow,
            startRow = 1,
            startColumn = 1,
            showNA = T
        )
        # an index to keep track of where to put the next mass being matched
        previousEnd = 2
        for (i in 1:dim(zd)[1]) {
            cb = CellBlock(
                sheet = sheet,
                startRow = previousEnd,
                startColumn = 1,
                noRows = 1,
                noColumns = 1,
                create = T
            )
            CB.setMatrixData(
                cellBlock = cb,
                x = as.matrix(paste(zd$rt[i], "_", zd$m[i], sep = "")),
                startRow = 1,
                startColumn = 1,
                showNA = T
            )
            # define cell block for M, H, Na matches and add data to the block
            # For Mass matches
            if (length(match_mass[[i]]) != 0) {
                massdata = as.matrix(antibase$Name[match_mass[[i]]])
                numberofRowsM = length(match_mass[[i]])
                cb = CellBlock(
                    sheet = sheet,
                    startRow = previousEnd,
                    startColumn = 2,
                    noRows = numberofRowsM,
                    noColumns = 1,
                    create = T
                )
                CB.setMatrixData(
                    cellBlock = cb,
                    x = massdata,
                    startRow = 1,
                    startColumn = 1,
                    showNA = T
                ) # write the data
            } else {
                numberofRowsM = 0
            }
            # For Hydrogen Matches
            if (length(match_hydrogen[[i]]) != 0) {
                hydrogendata = as.matrix(antibase$Name[match_hydrogen[[i]]])
                numberofRowsH = length(match_hydrogen[[i]])
                cb = CellBlock(
                    sheet = sheet,
                    startRow = previousEnd,
                    startColumn = 3,
                    noRows = numberofRowsH,
                    noColumns = 1,
                    create = T
                )
                CB.setMatrixData(
                    cellBlock = cb,
                    x = hydrogendata,
                    startRow = 1,
                    startColumn = 1,
                    showNA = T
                ) # write the data
            } else {
                numberofRowsH = 0
            }
            # For Sodium Matches
            if (length(match_sodium[[i]]) != 0) {
                sodiumdata = as.matrix(antibase$Name[match_sodium[[i]]])
                numberofRowsNa = length(match_sodium[[i]])
                cb = CellBlock(
                    sheet = sheet,
                    startRow = previousEnd,
                    startColumn = 4,
                    noRows = length(match_sodium[[i]]),
                    noColumns = 1,
                    create = T
                )
                CB.setMatrixData(
                    cellBlock = cb,
                    x = sodiumdata,
                    startRow = 1,
                    startColumn = 1,
                    showNA = T
                ) # write the data
            } else {
                numberofRowsNa = 0
            }
            if (max(numberofRowsNa, numberofRowsH, numberofRowsM) != 0) {
                previousEnd = previousEnd +
                    (max(numberofRowsNa,
                         numberofRowsH,
                         numberofRowsM))
            } else {
                cb = CellBlock(
                    sheet = sheet,
                    startRow = previousEnd,
                    startColumn = 2,
                    noRows = 1,
                    noColumns = 1,
                    create = T
                ) # cell block of 1 cell
                CB.setMatrixData(
                    cellBlock = cb,
                    x = matrix(c("NO MATCHES FOUND"), nrow = 1),
                    startRow = 1,
                    startColumn = 1,
                    showNA = T
                ) # say no matches found
                zd$Unique[i] = 1
                previousEnd = previousEnd + 1
            }
        }
    }
    # Write the matrix 'zd' as a data frame to the sheet for the current strain. This is the list of masses sorted in order of euclidean distance on the loadings plot. The last column includes a '1' for unique and a blank for something found in AntiBase.
    if (FLAG == 0) {
        # Excel sheet preparations
        wb = loadWorkbook(file = paste(bucket_table_name, ".xlsx", sep = ""))
        # make sheet = strain name
        sheet = createSheet(wb,
                            sheetName = rownames(buckets_unique)[r])
    }
    addDataFrame(
        x = as.data.frame(zd),
        sheet = sheet,
        col.names = T,
        row.names = F,
        startRow = 1,
        startColumn = 6,
        showNA = F
    )
    saveWorkbook(wb, file = paste(bucket_table_name, ".xlsx", sep = ""))
}