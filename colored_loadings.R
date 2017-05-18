#Draw Colored Loadings
colored_loadings <- function(PCA, PCx = 1, PCy = 2, npoints = 2000, ntimes = 3) {
    d = as.data.frame(PCA$rotation[, PCx:PCy])
    d = add_euclidean_distance(d)
    d = d[order(d$ed, decreasing = T), ]
    c = c(rainbow(ntimes + 1))
    v = c(rep(npoints, ntimes), dim(d)[1] - (ntimes * npoints))
    d$col = rep(x = c, times = v)
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
    legend(
        x = "topright",
        legend = v,
        fill = c,
        border = "black",
        # bg = "black",
        text.font = 2,
        x.intersp = 0.5,
        y.intersp = 0.5
        # cex = 0.75
    )
}
