ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, innerSegments = 4L, numThreads = 2L)

\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)

fitness(result)
}
