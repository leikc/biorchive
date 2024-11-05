data(iris)

head(iris)
pca.res <- prcomp(scale(iris[,1:4]))
plot(pca.res$x)