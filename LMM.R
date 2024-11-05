library(nlme);library(lme4)


LMMsetup <- function(form, dat, ref = list()) {
  # Construct fixed effects model matrix (X) from formula and data
  X <- model.matrix(form, dat)
  
  # Construct random effects model matrix (Z)
  Z_list <- lapply(ref, function(vars) model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "-1")), dat))
  Z <- do.call(cbind, Z_list)
  
  # Response variable
  y <- model.response(model.frame(form, dat))
  
  return(list(X = X, Z=Z,y=y))
}
LMMsetup(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
