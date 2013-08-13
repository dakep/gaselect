ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, numSegments = 4L, numThreads = 2L)

result <- genAlg(y, x, control = ctrl, evaluator = evaluator)

subsets(result, names = TRUE, indices = 1:5) # best 5 variable subsets as a list of names

result@subsets[ , 1:5] # best 5 variable subsets as a logical matrix with the subsets in the columns