library(nlme);library(lme4)
library(debug)

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
solve_chol <- function(L, b) {
  return (backsolve(L, forwardsolve(t(L),b)))
}
result<-LMMsetup(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
# Step 2: Initial guesses for theta: log(sigma) and log standard deviations for random effects
theta_init <- rep(0, length(list("Worker",c("Worker","Machine"))) + 1)  # Starting with a small positive value

LMMprof <- function(theta, setup) {
  # browser()
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  
    # QR Decomposition of Z
  qr_decomp <- qr(Z)
  Q <- qr.Q(qr_decomp, complete = TRUE)
  R <- qr.R(qr_decomp)
  
  # Dimensions
  n <- length(y)
  p <- ncol(Z)
  
  # Extract sigma and random effects variances from theta
  sigma <- exp(theta[1])#I think we don't need to expone
  psi_diag <- exp(2 * theta[-1])
  # browser()
  # Construct the covariance matrix Psi_b for the random effects
  Psi_b <- diag(psi_diag, ncol(Z), ncol(Z))
  small_block<-R%*%Psi_b%*%t(R)+diag(1, nrow = p, ncol = p)*(sigma^2)
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
  residual <- y - X %*% beta_hat
  minus_log_likelihood<- 0.5*t(residual)%*%qr.qy(qr_decomp, W_middle)%*%qr.qty(qr_decomp, residual)+sum(log(diag(small_block_Inv)))+ ((n - p) * log(sigma^2)/2)
  # log_likelihood <- -0.5 * (t(residual) %*% W %*% residual + sum(log(diag(chol(small_block)))) + (n - p) * log(sigma^2))
  
  print("minus_Log-likelihood:")
  print(minus_log_likelihood)
  attr(minus_log_likelihood, "beta_hat") <- beta_hat
  return(minus_log_likelihood)  # Return negative log-likelihood for minimization
}
mtrace(LMMprof, FALSE)
mtrace(lmm)


x<-LMMprof(theta = theta_init, setup = result)
attr(x, 'beta_hat')

lmm <- function(form, dat, ref = list()) {
  # Step 1: Setup model matrices and data
  setup <- LMMsetup(form, dat, ref)
  
  # Step 2: Initial guesses for theta: log(sigma) and log standard deviations for random effects
  theta_init <- rep(0, length(ref) + 1)  # Starting with a small positive value
  lower_bounds <- rep(log(.001)/2, length(ref) + 1)  # Example: ensuring all theta > -2
  upper_bounds <- rep(Inf, length(ref) + 1)  # No upper bounds, or set specifically if needed
  
  # Step 3: Optimize negative log-likelihood using `optim`
  # opt <- optim(theta_init, LMMprof, setup = setup, method = "L-BFGS-B", control = list(fnscale = 1))
  opt <- optim(theta_init, LMMprof, setup = setup, lower = lower_bounds,
               upper = upper_bounds, method = "L-BFGS-B", control = list())
  
  # Extract optimal theta and compute beta estimate
  theta_opt <- opt$par
  sigma_opt <- exp(theta_opt[1])
  psi_diag_opt <- exp(2 * theta_opt[-1])
  
  # # Recompute beta estimate using optimized theta
  # #fix this one hac
  # W <- chol2inv(chol(qr.R(qr(setup$Z)) %*% diag(psi_diag_opt, ncol(setup$Z), ncol(setup$Z)) %*% t(qr.R(qr(setup$Z))) + sigma_opt^2 * diag(length(setup$y))))
  # beta_hat <- solve(t(setup$X) %*% W %*% setup$X, t(setup$X) %*% W %*% setup$y)
  
  return(list(beta = beta_hat, theta = theta_opt))
}


# Load necessary libraries
library(nlme)
library(lme4)
library(MASS)
# Load the Machines dataset and use `lmm` function
data("Machines", package = "nlme")
result <- lmm(score ~ Machine, dat = Machines, ref = list("Worker", c("Worker", "Machine")))

# Compare to lme4 results
lmer_model <- lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), data = Machines, REML = FALSE)
summary(lmer_model)
print(result)



#.............................
# 𝑍
# Z:
#   
#   Suppose 
# 𝑍
# Z has dimensions 
# 𝑛
# ×
# 𝑝
# n×p, where:
#   𝑛
# n is the number of observations (rows).
# 𝑝
# p is the number of random effects (columns).
# So, 
# 𝑍
# Z is an 
# 𝑛
# ×
# 𝑝
# n×p matrix.
# Ψ
# 𝜃
# Ψ 
# θ
# 
# :
#   
#   Ψ
# 𝜃
# Ψ 
# θ
# 
# is a diagonal matrix representing the variances of the random effects.
# It has dimensions 
# 𝑝
# ×
# 𝑝
# p×p, where:
#   𝑝
# p is the number of random effects.
# Ψ
# 𝜃
# Ψ 
# θ
# 
# typically contains the variance parameters associated with each random effect along its diagonal.

# This final 
# 𝑛
# ×
# 𝑛
# n×n matrix is used as 
# 𝑊
# W, 