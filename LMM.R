
# The goal is to implement a linear mixed-effects model (LMM) estimation procedure based on maximum likelihood. 
# This involves deriving estimates for fixed effects (Betas) and variance components (theta's) by minimizing the negative log-likelihood.
# The process utilizes the QR decomposition of the random effects design matrix (Z) for efficiency, avoiding large matrix inversions. 
# The key steps include:
# 1- Identifying Y (Response Variable) , X (Fixed Effects) and Z (Random Effects) this is done within the LMMsetup function.
# 2- Formulating the likelihood function with respect to theta this is done in the LMMprof function.
# 3- Using Cholesky decomposition for matrix computations, this is done in the solve_chol function. 
# 4- Iteratively optimizing the log-likelihood to obtain parameter estimates this is done in the lmm function.
# L-BFGS-B optimization method is chosen because it efficiently handles bound constraints, ensuring variances are non-negative.
# it uses gradient information for faster convergence. It works well for models with random effects and simpler models without random effects, making it a versatile choice for this problem.

# Load necessary libraries
library(nlme)
library(lme4)
library(MASS)
LMMsetup <- function(form, dat, ref = list()) {
  # This function constructs the fixed-effects design matrix (X), random-effects 
  # design matrix (Z), and response variable (y) for a given linear mixed-effects model
  
  # Inputs :
  # form: A formula specifying the fixed effects part of the model (y ~ x1 + x2).
  # dat: A data frame containing the variables in the model.
  # ref: A list specifying the random effects structure.
  
  # Returns : 
  # A list with the following components :
  # X:Fixed-effects design matrix
  # Z:Random-effects design matrix
  # y:Response variable vector
  # dimensions:A list of column dimensions for each random effect grouping in Z
  
  # Construct fixed effects model matrix (X) from formula and data
  X <- model.matrix(form, dat)
  
  # Construct random effects model matrix (Z)
  # Z_list <- lapply(ref, function(vars) model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "-1")), dat))
  Z_list <- lapply(ref, function(vars) {
    # Create the model matrix
    model_mat <- model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "-1")), dat)
    
    # Get the dimensions of the matrix
    attr(model_mat, "columns") <- dim(model_mat)[2]
    return(model_mat)
  })
  
  Z <- do.call(cbind, Z_list)
  dimensions<-lapply(Z_list, function(Z) attr(Z, "columns"))
  
  # Response variable
  y <- model.response(model.frame(form, dat))
  
  return(list(X = X, Z=Z,y=y, dimensions=dimensions))
}
solve_chol <- function(L, b) {
  # This function solves a linear system of equations (Ax = b) using 
  # The Cholesky decomposition of A, where A is symmetric positive-definite. 
  # It is an efficient method for solving such systems without explicitly computing the inverse of A.
  
  # Inputs :
  # L: A upper triangular matrix from the Cholesky decomposition of A such that A=L^TL.
  # b: A numeric vector representing the right-hand side of the linear system.
  
  # Returns : 
  # A numeric vector x, the solution to the linear system Ax=b.
  
  return (backsolve(L, forwardsolve(t(L),b)))
}
LMMprof <- function(theta, setup) {
  
  #"This function computes the negative log-likelihood for a linear mixed model 
  # given a set of parameters (theta) and model matrices. It handles cases 
  # with or without random effects by formulating the likelihood function and using 
  # efficient matrix computations such as QR decomposition and Cholesky decomposition.
  
  #Inputs : 
  # theta : A numeric vector of parameters to be estimated:
  # The first element corresponds to the log of the residual variance.
  # Subsequent elements correspond to the log-transformed variances of the random effects.
  
  # setup : A list of model components, as created by the `LMMsetup` function:
  
  # Returns : 
  # A scalar representing the negative log-likelihood for the given parameters. 
  # The estimated fixed effects (beta) are returned as an attribute "beta_hat". "
  
  # browser()
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  qr_decomp<-setup$qr_decomp
  R<-setup$R
  if(is.null(Z)){
    sigma <- exp(theta[1])
    n <- length(y)
    p<-ncol(R)
    beta_hat<-backsolve(R, qr.qty(qr_decomp, y)[1:p])
    residual <- y - (X %*% beta_hat)
    # Minus log-likelihood calculation
    minus_log_likelihood <- 0.5 * ((t(residual)%*%residual/sigma^2)+log(sigma^(2*n)))
    
    # print("minus_Log-likelihood:")
    # print(minus_log_likelihood)
    attr(minus_log_likelihood, "beta_hat") <- beta_hat
    return(minus_log_likelihood)
    
    
  }else{
    Z_cols<-setup$dimensions
    # Dimensions
    n <- length(y)
    p <- ncol(Z)
    
    # Extract sigma and random effects variances from theta
    sigma <- exp(theta[1])#I think we don't need to expone
    psi_diag <- exp(2 * theta[-1])
    # browser()
    # Construct the covariance matrix Psi_b for the random effects
    Psi_b <- diag(rep(psi_diag, Z_cols), ncol(Z), ncol(Z))
    small_block<-R%*%Psi_b%*%t(R)+diag(1, nrow = p, ncol = p)*(sigma^2)
    # Cholesky decomposition
    small_block_chol <- chol(small_block)
    
    small_block_Inv<-solve_chol(small_block_chol, diag(x = 1, nrow = p, ncol = p))
    # Step 1: Compute Q^T y and split it into two parts
    QTy <- qr.qty(qr_decomp, y) # Compute Q^T y
    QTy_1 <- QTy[1:p]           # First p rows
    QTy_2 <- QTy[(p + 1):n]     # Remaining (n - p) rows
    
    # Step 2: Compute W * y without forming W
    Wy_1 <- small_block_Inv%*%QTy_1 # Apply small block inverse to QTy_1
    Wy_2 <- QTy_2/sigma^(2)                     # Scale QTy_2 by 1/sigma^2
    WQTy <- c(Wy_1, Wy_2)                         # Combine to form Wy
    
    # Step 3: Compute X^T W X and X^T W y without splitting Q^T X
    QTX <- qr.qty(qr_decomp, X)  # Compute Q^T X without splitting
    
    # Apply the weights to QTX
    WX_1 <- small_block_Inv %*% QTX[1:p, , drop = FALSE]  # Apply the small block inverse to the first p rows
    WX_2 <- QTX[(p + 1):n, , drop = FALSE]/sigma^(2)              # Scale the remaining rows by 1 / sigma^2
    WQTX <- rbind(WX_1, WX_2)                                          # Combine weighted parts without further splitting
    
    # Compute X^T W X
    XWX <- t(X) %*% qr.qy(qr_decomp, WQTX)           
    
    # Compute X^T W y
    XWy<-t(X) %*% qr.qy(qr_decomp, WQTy)
    XWX_chol <- chol(XWX)
    beta_hat<-solve_chol(L=XWX_chol, b=XWy)
    residual <- y - (X %*% beta_hat)
    
    # Step 4: Estimate beta using the Cholesky decomposition of X^T W X
    XWX_chol <- chol(XWX)
    beta_hat <- solve_chol(XWX_chol, XWy)
    
    # Compute residuals and negative log-likelihood
    residual <- y - (X %*% beta_hat)
    QTy_res <- qr.qty(qr_decomp, residual)
    QTy_res_1 <- QTy_res[1:p]
    QTy_res_2 <- QTy_res[(p + 1):n]
    
    #W_res_1 <- solve_chol(small_block_chol, QTy_res_1)
    W_res_1 <- small_block_Inv %*% QTy_res_1
    W_res_2 <- QTy_res_2 / sigma^2
    W_res <- c(W_res_1, W_res_2)
    
    # Minus log-likelihood calculation
    minus_log_likelihood <- 0.5 * (t(residual)%*%qr.qy(qr_decomp, W_res) + 2*sum(log(diag(small_block_chol))) + (n - p) * log(sigma^2))
    
    # print("minus_Log-likelihood:")
    # print(minus_log_likelihood)
    attr(minus_log_likelihood, "beta_hat") <- beta_hat
    return(minus_log_likelihood)
  }
}
# “Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”,“Brent
lmm <- function(form, dat, ref = list()) {
  # This function estimates the parameters of a Linear Mixed-Effects Model (LMM) using maximum likelihood.
  # It prepares model matrices for fixed and random effects, passes those data to LMMprof for optimization, and calculates
  # variance components and fixed effects by minimizing the negative log-likelihood.
  
  # Inputs :
  # form: A formula specifying the fixed effects part of the model (y ~ x1 + x2).
  # dat: A data frame containing the variables in the model.
  # ref: A list specifying the random effects structure.
  
  # Returns :
  # A list containing :
  # beta: Estimated fixed effect coefficients. 
  # theta: Optimized variance components (residual variance and random effects variances)
  # Step 1: Setup model matrices and data
  setup <- LMMsetup(form, dat, ref)
  # QR Decomposition of Z
  if(is.null(setup$Z)){
    setup$qr_decomp <- qr(setup$X)
    setup$R <- qr.R(setup$qr_decomp)
    
    
  }else{
    setup$qr_decomp <- qr(setup$Z)
    setup$R <- qr.R(setup$qr_decomp)
  }
  
  
  # Step 2: Initial guesses for theta: log(sigma) and log standard deviations for random effects
  # theta_init<-rnorm(length(ref) + 1)
  theta_init <- rep(0, length(ref) + 1)  # Starting with a small positive value
  # lower_bounds <- rep(-Inf, length(ref) + 1)  # Example: ensuring all theta > -2
  # upper_bounds <- rep(Inf, length(ref) + 1)  # No upper bounds, or set specifically if needed
  # Step 3: Optimize negative log-likelihood using `optim`
  opt <- optim(theta_init, LMMprof, setup = setup, method = "L-BFGS-B", control = list(fnscale = 1))
  # opt <- optim(theta_init, LMMprof, setup = setup, lower = lower_bounds,
  #              upper = upper_bounds, method = "L-BFGS-B", control = list())
  final_cost_value <- LMMprof(theta = opt$par, setup = setup)
  
  # Accessing the attribute "beta"
  beta_hat <- attr(final_cost_value, "beta_hat")
  
  
  return(list(beta = beta_hat, theta = opt$par))
}
# Load the Machines dataset and use `lmm` function
data("Machines", package = "nlme")
result <- lmm(score ~ Machine, dat = Machines, ref = list("Worker", c("Worker", "Machine")))
# Compare to lme4 results
lmer_model <- lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), data = Machines, REML = FALSE)
summary(lmer_model)
print(exp(result$theta))
print(result$beta)

#compare to lm results
result <- lmm(score ~ Machine, dat = Machines, ref = list())
lm_model<-lm(score ~ Machine, data = Machines)
predictions <- predict(lm_model, Machines)
sqrt(sum((Machines$score-predictions)^2)/dim(Machines)[1])
summary(lm_model)
print(exp(result$theta))
print(result$beta)
