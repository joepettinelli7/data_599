# caret train() and rfe() use kernlab for SVM by default which was not
# sufficient with the radial basis function kernel because ksvm algorithm
# could not find SVM solution.
# http://topepo.github.io/caret/using-your-own-model-in-train.html

e1071_rbf_svm <- list(label = "Custom RBF SVM with e1071",
                       type = "Classification",
                       library = "e1071",
                       loop = NULL) 

# Parameters
e1071_rbf_svm$parameters <- data.frame(parameter = c("cost", "gamma"),
                                        class = rep("numeric", 2),
                                        label = c("Cost", "Gamma"))

# Grid
e1071_rbf_svm$grid <- function (x, y, len = NULL, search = "grid") {
  warning('Will need to implement!')
  browser()
}

# Loop
e1071_rbf_svm$loop <- NULL

# Fit
e1071_rbf_svm$fit <- function (x, y, wts, param, lev, last, classProbs, ...) {
  out <- e1071::svm(x = x, y = y, kernel = "radial", cost = param$cost,
                    gamma = param$gamma, probability = classProbs,
                    type = "C-classification", tol = 0.001, ...)
  return(out)
}

# Prediction
e1071_rbf_svm$predict <- function (modelFit, newdata, submodels = NULL) {
  
  svmPred <- function(obj, x) {
    pred <- predict(obj, x, type = "class") 
    return(pred)
  }
  
  out <- try(svmPred(modelFit, newdata), silent = TRUE)
  if (class(out)[1] == "try-error") {
    warning("e1071 prediction calculations failed; returning NAs")
    out <- rep(NA, nrow(newdata))
  }
  return(out)
}

# Probability
e1071_rbf_svm$prob <- function (modelFit, newdata, submodels = NULL) {
  out <- try(predict(modelFit, newdata, probability = TRUE), silent = TRUE)
  out <- attr(out, "probabilities")
  if (class(out)[1] != "try-error") {
    # If there are any negative probabilities, set to 0 then
    # normalize probabilities to sum to 1
    if (any(out < 0)) {
      out[out < 0] <- 0
      out <- t(apply(out, 1, function(x) x / sum(x)))
    }
    out <- out[, modelFit$levels, drop = FALSE]
  } else {
    num_lvls <- length(modelFit$levels)
    warning("e1071 class probability calculations failed; returning NAs")
    out <- matrix(NA, nrow = nrow(newdata) * num_lvls, ncol = num_lvls)
    colnames(out) <- modelFit$levels
  }
  return(out)
}

# Predictors
e1071_rbf_svm$predictors <- function (x, ...) {
  warning("Will need to implement!")
  browser()
}

# Levels needed for s4
e1071_rbf_svm$levels <- function (x) x$levels

# Sort
e1071_rbf_svm$sort <- function (x) {
  x[order(x$cost, -x$gamma), ]
}
