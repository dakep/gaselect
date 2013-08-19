ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)

# Use the BIC of a linear model to evaluate the fitness of a variable subset
evalFun <- function(y, X) {
		return(BIC(lm(y ~ X)));
}

# Dummy function that returns the residuals standard deviation and not the SEP
sepFUN <- function(genAlg) {
		return(apply(genAlg@subsets, 2, function(subset) {
			m <- lm(genAlg@response ~ genAlg@covariates[, subset]);
			return(sd(m$residuals));
		}));
}
evaluator <- evaluatorUserFunction(FUN = evalFun, sepFUN = sepFUN)
\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
}
