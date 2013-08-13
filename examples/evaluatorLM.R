ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorLM(statistic = "adjusted.r.squared", maxCor = 0.8)

result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
