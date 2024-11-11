# Load necessary libraries
library(debug)
library(nlme)
library(lme4)
library(MASS)
LMMsetup <- function(form, dat, ref = list()) {
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
setup <- LMMsetup(score ~ Machine, dat = Machines, ref = list("Worker", c("Worker", "Machine")))
solve_chol <- function(L, b) {
  return (backsolve(L, forwardsolve(t(L),b)))
}
LMMprof <- function(theta, setup) {
  # browser()
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  Z_cols<-setup$dimensions
  qr_decomp<-setup$qr_decomp
  R<-setup$R
  # Dimensions
  n <- length(y)
  p <- ncol(Z)
  
  # Extract sigma and random effects variances from theta
  sigma <- exp(theta[1])#I think we don't need to expone
  psi_diag <- exp(2 * theta[-1])
  # browser()
  # Construct the covariance matrix Psi_b for the random effects
  Psi_b <- diag(rep(psi_diag, Z_cols), ncol(Z), ncol(Z))
  small_block<-(R%*%Psi_b%*%t(R))+(diag(1, nrow = p, ncol = p)*(sigma^2))
  # Cholesky decomposition
  small_block_chol <- chol(small_block)
  
  small_block_Inv<-solve_chol(small_block_chol, diag(x = 1, nrow = p, ncol = p))
  I_n_p <- diag(n - p)/(sigma^(2))
  
  # Construct the full W matrix
  W_middle <- rbind(
    cbind(small_block_Inv, matrix(0, nrow = p, ncol = n - p)),
    cbind(matrix(0, nrow = n - p, ncol = p), I_n_p)
  )
  
  # Calculate W using Q and the block matrix
  
  # W3 <- qr.qy(qr_decomp, t(qr.qy(qr_decomp, t(W_middle))))
  # W2<-Q %*% W_middle %*% t(Q)
  # Compute log likelihood terms
  XWX<-t(X)%*%qr.qy(qr_decomp, W_middle)%*%qr.qty(qr_decomp, X)
  XWy<-t(X)%*%qr.qy(qr_decomp, W_middle)%*%qr.qty(qr_decomp, y)
  XWX_chol <- chol(XWX)
  beta_hat<-solve_chol(L=XWX_chol, b=XWy)
  residual <- y - (X%*%beta_hat)
  minus_log_likelihood<- 0.5*(t(residual)%*%qr.qy(qr_decomp, W_middle)%*%qr.qty(qr_decomp, residual)+2*sum(log(diag(small_block_chol)))+ ((n - p) * log(sigma^2)))
  # log_likelihood <- -0.5 * (t(residual) %*% W %*% residual + sum(log(diag(chol(small_block)))) + (n - p) * log(sigma^2))
  
  print("minus_Log-likelihood:")
  print(minus_log_likelihood)
  attr(minus_log_likelihood, "beta_hat") <- beta_hat
  return(minus_log_likelihood)  # Return negative log-likelihood for minimization
}
lmm <- function(form, dat, ref = list()) {
  
  # Step 1: Setup model matrices and data
  setup <- LMMsetup(form, dat, ref)
  # QR Decomposition of Z
  setup$qr_decomp <- qr(setup$Z)
  setup$R <- qr.R(setup$qr_decomp)
  
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
# mtrace(lmm, FALSE)
# mtrace(LMMprof, FALSE)
# Load the Machines dataset and use `lmm` function
data("Machines", package = "nlme")
result <- lmm(score ~ Machine, dat = Machines, ref = list("Worker", c("Worker", "Machine")))
print(exp(result$theta))
print(result$beta)

# Compare to lme4 results
lmer_model <- lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), data = Machines, REML = FALSE)
summary(lmer_model)



