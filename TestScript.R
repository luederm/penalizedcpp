require(microbenchmark)
require(penalizedcpp)
require(penalized)
data(nki70)

# Function which compares two penfit objects
compPenfitObj <- function(fit1, fit2, i = 0) {
  passed = T
  
  # penalized
  pen = all.equal(fit1@penalized, fit2@penalized, check.names = T, check.attributes = T)
  if (!isTRUE(pen)){ warning(paste("Difference detected comparing penalized slot:", pen)); passed = F}
  
  # unpenalized
  upen = all.equal(fit1@unpenalized, fit2@unpenalized, check.names = T, check.attributes = T)
  if (!isTRUE(upen)){ warning(paste("Difference detected comparing unpenalized slot:", upen)); passed = F}
  
  # residuals
  res = all.equal(fit1@residuals, fit2@residuals, check.names = T, check.attributes = T)
  if (!isTRUE(res)){ warning(paste("Difference detected comparing residuals slot:", res)); passed = F}
  
  # fitted
  fitted = all.equal(fit1@fitted, fit2@fitted, check.names = T, check.attributes = T)
  if (!isTRUE(fitted)){ warning(paste("Difference detected comparing fitted slot:", fitted)); passed = F}
  
  # lin.pred
  lp = all.equal(fit1@lin.pred, fit2@lin.pred, check.names = T, check.attributes = T)
  if (!isTRUE(lp)){ warning(paste("Difference detected comparing lin.pred slot:", lp)); passed = F}
  
  # loglik
  loglik = all.equal(fit1@loglik, fit2@loglik, check.names = T, check.attributes = T)
  if (!isTRUE(loglik)){ warning(paste("Difference detected comparing loglik slot:", loglik)); passed = F}
  
  # penalty
  penalty = all.equal(fit1@penalty, fit2@penalty, check.names = T, check.attributes = T)
  if (!isTRUE(penalty)){ warning(paste("Difference detected comparing penalty slot:", penalty)); passed = F}
  
  # iterations
  iter = all.equal(fit1@iterations, fit2@iterations, check.names = T, check.attributes = T)
  if (!isTRUE(iter)){ warning(paste("Difference detected comparing iterations slot:", iter)); passed = F}
  
  # converged
  converged = all.equal(fit1@converged, fit2@converged, check.names = T, check.attributes = T)
  if (!isTRUE(converged)){ warning(paste("Difference detected comparing converged slot:", converged)); passed = F}
  
  # model
  model = all.equal(fit1@model, fit2@model, check.names = T, check.attributes = T)
  if (!isTRUE(model)){ warning(paste("Difference detected comparing model slot:", model)); passed = F}
  
  # lambda1
  l1 = all.equal(fit1@lambda1, fit2@lambda1, check.names = T, check.attributes = T)
  if (!isTRUE(l1)){ warning(paste("Difference detected comparing lambda1 slot:", l1)); passed = F}
  
  # lambda2
  l2 = all.equal(fit1@lambda2, fit2@lambda2, check.names = T, check.attributes = T)
  if (!isTRUE(l2)){ warning(paste("Difference detected comparing lambda2 slot:", l2)); passed = F}
  
  # weights
  weights = all.equal(fit1@weights, fit2@weights, check.names = T, check.attributes = T)
  if (!isTRUE(weights)){ warning(paste("Difference detected comparing weights slot:", weights)); passed = F}
  
  return(passed)
}

# Function which compares output of opt/prof functions
compareCVOutput <- function(fit1, fit2) {
  # lambda
  lambda = all.equal(fit1$lambda, fit2$lambda, check.names = T, check.attributes = T)
  if (!isTRUE(lambda)) warning(paste("Difference detected comparing lambda element:", lambda))
  
  # cvl
  cvl = all.equal(fit1$cvl, fit2$cvl, check.names = T, check.attributes = T)
  if (!isTRUE(cvl)) warning(paste("Difference detected comparing cvl element:", cvl))
  
  # fold
  fold = all.equal(fit1$fold, fit2$fold, check.names = T, check.attributes = T)
  if (!isTRUE(fold)) warning(paste("Difference detected comparing cvl element:", fold))
  
  # predictions
  fold = all.equal(fit1$fold, fit2$fold, check.names = T, check.attributes = T)
  if (!isTRUE(fold)) warning(paste("Difference detected comparing cvl element:", fold))
  
  # fullfit
  if ( (is.list(fit1$fullfit) != is.list(fit2$fullfit)) | (length(fit1$fullfit) != length(fit2$fullfit)) )
  {
    warning("Mismatch between fullfit types or length")
    
  } else if (is.list(fit1$fullfit)) {
    for (i in 1:length(fit1$fullfit))
    {
      pass = compPenfitObj(fit1$fullfit[[i]], fit2$fullfit[[i]])
      if (!pass) print(paste("Failures occured at fullfit index", i))
    }
  } else {
    compPenfitObj(fit1$fullfit, fit2$fullfit)
  }
}

# ***** COMPARE RESULTS *****

# Cox model w/ lasso
result1 = penalized:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
              data = nki70, standardize=TRUE, lambda1 = 0.5)
result2 = penalizedcpp:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
              data = nki70, standardize=TRUE, lambda1 = 0.5)
compPenfitObj(result1, result2)

# Cox model w/ ridge
result1 = penalized:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                    data = nki70, standardize=TRUE, lambda2 = 0.5)
result2 = penalizedcpp:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                    data = nki70, standardize=TRUE, lambda2 = 0.5)
compPenfitObj(result1, result2)

# Cox model w/ l1 penalty optimization
set.seed(42);
result1 = penalized:::optL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::optL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Cox model w/ l1 profile
set.seed(42);
result1 = penalized:::profL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::profL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Cox model w/ both l1 and l2 penalty and n < p
result1 = penalized:::penalized(Surv(nki70[1:20,1], nki70[1:20,2]), penalized = nki70[1:20,8:77], data = nki70,
                        lambda1 = .001, lambda2 = .002, fusedl = F)
result2 = penalizedcpp:::penalized(Surv(nki70[1:20,1], nki70[1:20,2]), penalized = nki70[1:20,8:77], data = nki70,
                        lambda1 = .001, lambda2 = .002, fusedl = F)
compPenfitObj(result1, result2)

# Logistic regression w/ lasso
result1 = penalized:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda1=3)
result2 = penalizedcpp:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda1=3)
compPenfitObj(result1, result2)

# Logistic regression w/ ridge
result1 = penalized:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda2=3)
result2 = penalizedcpp:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda2=3)
compPenfitObj(result1, result2)

# Logistic regression w/ l1 penalty optimization
set.seed(42);
result1 = penalized:::optL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::optL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Logistic regression w/ l1 profile
set.seed(42);
result1 = penalized:::profL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::profL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Logistic regression model w/ both l1 and l2 penalty and n < p
result1 = penalized:::penalized(nki70[1:20,5], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
result2 = penalizedcpp:::penalized(nki70[1:20,5], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
compPenfitObj(result1, result2)

# Linear model w/ lasso
result1 = penalized:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda1=3)
result2 = penalizedcpp:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda1=3)
compPenfitObj(result1, result2)

# Linear model w/ ridge
result1 = penalized:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda2=3)
result2 = penalizedcpp:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda2=3)
compPenfitObj(result1, result2)

# Linear model w/ l1 penalty optimization
set.seed(42);
result1 = penalized:::optL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::optL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Linear model w/ l1 profile
set.seed(42);
result1 = penalized:::profL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5)
set.seed(42);
result2 = penalizedcpp:::profL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5)
compareCVOutput(result1, result2)

# Linear model w/ both l1 and l2 penalty and n < p
result1 = penalized:::penalized(nki70[1:20,1], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
result2 = penalizedcpp:::penalized(nki70[1:20,1], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
compPenfitObj(result1, result2)

# High dimensional analysis data set
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36000
#source("http://bioconductor.org/biocLite.R") 
#biocLite("affy")
library(affy)
setwd("/Users/matthewlueder/desktop/exprDat/")
eset <- justRMA()
exprs <- exprs(eset)
set.seed(1234)
CR_event <- sample( c(0,1,2), 76, replace=TRUE )
hdim <- cbind( nki70[1:76,1:2], CR_event, nki70[1:76,3:length(nki70)], t(exprs))

# H-Dim Cox model w/ lasso
result1 = penalized:::penalized(Surv(time, event), penalized = hdim[,9:54753],
                data = hdim, standardize=TRUE, lambda1 = 3)
result2 = penalizedcpp:::penalized(Surv(time, event), penalized = hdim[,9:54753],
                data = hdim, standardize=TRUE, lambda1 = 3)
compPenfitObj(result1, result2)

# H-Dim competing risks model w/ lasso
result1 = penalized:::penalized(Surv(time, CR_event == 1), penalized = hdim[,9:54753],
                    data = hdim, standardize=TRUE, lambda1 = 3)
result2 = penalizedcpp:::penalized(Surv(time, CR_event == 1), penalized = hdim[,9:54753],
                    data = hdim, standardize=TRUE, lambda1 = 3)
compPenfitObj(result1, result2)


# ***** COMPARE SPEED *****

# Cox model w/ lasso
microbenchmark(
  penalized:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                      data = nki70, standardize=TRUE, lambda1 = 0.5),
  penalizedcpp:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                      data = nki70, standardize=TRUE, lambda1 = 0.5)
)
#   min       lq      mean    median     uq      max     neval
# 554.3666 584.2260 699.3707 712.1897 745.1892 1304.951   100
# 224.7664 233.5997 264.8149 244.1310 253.3629  470.378   100

# Cox model w/ ridge
microbenchmark(
  penalized:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                    data = nki70, standardize=TRUE, lambda2 = 0.5),
  penalizedcpp:::penalized(Surv(time, event)~strata(ER), penalized = nki70[,8:76],
                    data = nki70, standardize=TRUE, lambda2 = 0.5)
)
#  min       lq     mean   median       uq      max neval
#30.26409 37.73726 40.98312 39.27207 40.28391 221.2033   100
#29.74334 35.88285 40.34717 36.99988 38.36700 214.6791   100


# Cox model w/ l1 penalty optimization
microbenchmark(
  penalized:::optL1(Surv(time, event)~strata(ER), penalized = nki70[,8:77], unpenalized = nki70[,77],
        data = nki70, standardize=TRUE, fold = 5),
  penalizedcpp:::optL1(Surv(time, event)~strata(ER), penalized = nki70[,8:77], unpenalized = nki70[,77],
        data = nki70, standardize=TRUE, fold = 5),
  times = 10
)
#       min       lq     mean   median       uq      max neval
# 1533.1250 1762.368 2410.892 1977.835 2845.067 4991.450    10
#  970.0893 1028.945 1324.902 1078.799 1116.879 3121.993    10

# Cox model w/ l1 profile
microbenchmark(
  penalized:::profL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5),
  penalizedcpp:::profL1(Surv(time, event)~strata(ER), penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5),
  times = 10
)

# Cox model w/ both l1 and l2 penalty and n < p
microbenchmark(
  penalized:::penalized(Surv(nki70[1:20,1], nki70[1:20,2]), penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F),
  penalizedcpp:::penalized(Surv(nki70[1:20,1], nki70[1:20,2]), penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
)

# Logistic regression w/ lasso
microbenchmark(
  penalized:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda1=3),
  penalizedcpp:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda1=3)
)

# Logistic regression w/ ridge
microbenchmark(
  penalized:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda2=3),
  penalizedcpp:::penalized(ER, penalized = nki70[,8:76], data=nki70, lambda2=3)
)

# Logistic regression w/ l1 penalty optimization
microbenchmark(
  penalized:::optL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5),
  penalizedcpp:::optL1(ER, penalized = nki70[,8:76], unpenalized = nki70[,77],
                data = nki70, standardize=TRUE, fold = 5),
  times = 20
)

# Linear model w/ lasso
microbenchmark(
  penalized:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda1=3),
  penalizedcpp:::penalized(time, penalized = nki70[,8:76], data=nki70, lambda1=3)
)

# Linear model w/ l1 profile
microbenchmark(
  penalized:::profL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5),
  penalizedcpp:::profL1(time, penalized = nki70[,8:76], unpenalized = nki70[,77],
                 data = nki70, standardize=TRUE, fold = 5),
  times = 10
)
  
# Linear model w/ both l1 and l2 penalty and n < p
microbenchmark(
  penalized:::penalized(nki70[1:20,1], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F),
  penalizedcpp:::penalized(nki70[1:20,1], penalized = nki70[1:20,8:77], data = nki70,
                    lambda1 = .001, lambda2 = .002, fusedl = F)
)

# H-Dim Cox model w/ lasso
microbenchmark(
  penalized:::penalized(Surv(time, event), penalized = hdim[,9:54753],
              data = hdim, standardize=TRUE, lambda1 = 3),
  penalizedcpp:::penalized(Surv(time, event), penalized = hdim[,9:54753],
              data = hdim, standardize=TRUE, lambda1 = 3),
  times = 5
)

# H-Dim competing risks model w/ lasso
microbenchmark(
  penalized:::penalized(Surv(time, CR_event == 1), penalized = hdim[,9:54753],
                    data = hdim, standardize=TRUE, lambda1 = 3),
  penalizedcpp:::penalized(Surv(time, CR_event == 1), penalized = hdim[,9:54753],
                    data = hdim, standardize=TRUE, lambda1 = 3),
  times = 5
)




