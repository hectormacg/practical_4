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
dim(result$Z)
qrZ<-qr(result$Z)
dim(qr.R(qrZ))
LMMprof<-function(X, Z, y){
          qrZ<-qr(Z)
  
}
result<-LMMsetup(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
dim(resul$Z)
unique(Machines$Worker)

LMMprof <- function(theta, setup) {
  X <- setup$X
  Z <- setup$Z
  y <- setup$y
  
  # Extract sigma and random effects variances from theta
  sigma <- exp(theta[1])#I think we don't need to expone
  psi_diag <- exp(2 * theta[-1])
  
  # Construct the covariance matrix Psi_b for the random effects
  Psi_b <- diag(psi_diag, ncol(Z), ncol(Z))
  
  # QR Decomposition of Z
  qr_decomp <- qr(Z)
  Q <- qr.Q(qr_decomp, complete = TRUE)
  R <- qr.R(qr_decomp)
  
  # Dimensions
  n <- length(y)
  p <- ncol(Z)
  
  # Subset R to match dimensions of Psi_b
  R_subset <- R[1:p, ]  # Ensure R has compatible dimensions with Psi_b
  
  # Calculate small_block with regularization
  regularization <- 1e-6  # Small positive value to ensure positive definiteness
  small_block <- R_subset %*% Psi_b %*% t(R_subset) + sigma^2 * diag(p) + regularization * diag(p)
  
  # # Check eigenvalues of small_block for positive definiteness
  # eigenvalues <- eigen(small_block, symmetric = TRUE, only.values = TRUE)$values
  # print("Current theta values:")
  # print(theta)
  # print("Eigenvalues of small_block:")
  # print(sort(eigenvalues))
   
  # If positive definite, proceed with Cholesky inversion
  small_block_inv <- chol2inv(chol(small_block))
  
  # Construct the full block matrix
  block_matrix <- rbind(
    cbind(small_block_inv, matrix(0, p, n - p)),
    cbind(matrix(0, n - p, p), diag(1 / sigma^2, n - p))
  )
  
  # Calculate W using Q and the block matrix
  W <- Q %*% block_matrix %*% t(Q)
  
  # Compute log likelihood terms
  beta_hat <- solve(t(X) %*% W %*% X, t(X) %*% W %*% y)
  residual <- y - X %*% beta_hat
  log_likelihood <- -0.5 * (t(residual) %*% W %*% residual + sum(log(diag(chol(small_block)))) + (n - p) * log(sigma^2))
  
  print("Log-likelihood:")
  print(log_likelihood)
  
  return(-log_likelihood)  # Return negative log-likelihood for minimization
}
lmm <- function(form, dat, ref = list()) {
  # Step 1: Setup model matrices and data
  setup <- LMMsetup(form, dat, ref)
  
  # Step 2: Initial guesses for theta: log(sigma) and log standard deviations for random effects
  theta_init <- rep(0, length(ref) + 1)  # Starting with a small positive value
  
  # Step 3: Optimize negative log-likelihood using `optim`
  opt <- optim(theta_init, LMMprof, setup = setup, method = "L-BFGS-B", control = list(fnscale = 1))
  
  # Extract optimal theta and compute beta estimate
  theta_opt <- opt$par
  sigma_opt <- exp(theta_opt[1])
  psi_diag_opt <- exp(2 * theta_opt[-1])
  
  # Recompute beta estimate using optimized theta
  #fix this one hac
  W <- chol2inv(chol(qr.R(qr(setup$Z)) %*% diag(psi_diag_opt, ncol(setup$Z), ncol(setup$Z)) %*% t(qr.R(qr(setup$Z))) + sigma_opt^2 * diag(length(setup$y))))
  beta_hat <- solve(t(setup$X) %*% W %*% setup$X, t(setup$X) %*% W %*% setup$y)
  
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