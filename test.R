data(iris)

head(iris)
pca.res <- prcomp(scale(iris[,1:4]))
plot(pca.res$x)

library(matrix)
x <- matrix(1,1,1)

# Test message
message("Hello world!")

# Add another message and folder
message("Add folder")