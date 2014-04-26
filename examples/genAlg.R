ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluatorSRCV <- evaluatorPLS(numReplications = 5L, innerSegments = 10L, testSetSize = 0.4, numThreads = 2L)
evaluatorRDCV <- evaluatorPLS(numReplications = 5L, innerSegments = 7L, outerSegments = 4L, numThreads = 2L)

\dontrun{
resultSRCV <- genAlg(y, x, control = ctrl, evaluator = evaluatorSRCV)

resultRDCV <- genAlg(y, x, control = ctrl, evaluator = evaluatorRDCV)
}
