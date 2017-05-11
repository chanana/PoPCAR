setwd("~/Programming/R/PoPCAR/")
source("read_data.R")
filename = "~/Lab/bucket_tables/12Strains.txt"
a = get_spectral_file(filename)
b = get_spectral_data(a)
c = cleanup_cols(b, "s")
d = cleanup_rows(dataframe = c, dirty_rownames = a[,1])
