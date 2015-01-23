ctrl <- genAlgControl(populationSize = 200, numGenerations = 30, minVariables = 5,
    maxVariables = 12, verbosity = 1)

evaluatorSRCV <- evaluatorPLS(numReplications = 2L, innerSegments = 7L, testSetSize = 0.4,
    numThreads = 2L)

evaluatorRDCV <- evaluatorPLS(numReplications = 2L, innerSegments = 7L, outerSegments = 4L,
    numThreads = 2L)

# Generate demo-data
set.seed(12345)
X <- matrix(rnorm(10000, sd = 1:5), ncol = 50, byrow = TRUE)
y <- drop(-1.2 + rowSums(X[, seq(1, 43, length = 8)]) + rnorm(nrow(X), 1.5));

resultSRCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorSRCV, seed = 123)
resultRDCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorRDCV, seed = 123)

subsets(resultSRCV, 1:5)
subsets(resultRDCV, 1:5)
