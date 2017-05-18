#Draw Colored Loadings
colored_loadings <- function(PCA, PCx = 1, PCy = 2, npoints = 2000) {
    d = as.data.frame(PCA$rotation[, PCx:PCy])
    d = add_euclidean_distance(d)
    d = d[order(d$ed, decreasing = T), ]
    v = c(rep(npoints ,3), dim(d)[1] - (3 * npoints))
    d$col = rep(c("red", "blue", "green", "orange"), v)
    par(mfrow = c(1, 2))
    plot(
        x = d[, 1],
        y = d[, 2],
        pch = 16,
        col = d$col,
        xlab = colnames(d)[1],
        ylab = colnames(d)[2],
        main = "Loadings Plot"
    )
    plot(
        x = d$ed,
        pch = 16,
        col = d$col,
        xlab = "Index",
        ylab = "Euclidean Distance",
        main = "Euclidean Distance vs Index"
    )
}
