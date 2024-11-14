# The goal is to implement a linear mixed-effects model (LMM) estimation procedure based on maximum likelihood. 
# This involves deriving estimates for fixed effects (Betas) and log standar deviation components (theta's) by minimizing the negative log-likelihood.
# The process utilizes the QR decomposition of the random effects design matrix (Z) for efficiency, avoiding large matrix inversions. 
# The key steps include:
# 1- Identifying Y (Response Variable) , X (Fixed Effects) and Z (Random Effects) this is done within the LMMsetup function.
# 2- Formulating the likelihood function with respect to theta and finding the Betas this is done in the LMMprof function.
# 3- Using Cholesky decomposition for matrix computations, this is done in the solve_chol function. 
# 4- Iteratively optimizing the log-likelihood to obtain parameter estimates this is done in the lmm function.
# L-BFGS-B optimization method is chosen because it efficiently handles bound constraints, ensuring variances are non-negative.
# it uses gradient information for faster convergence. It works well for models with random effects and simpler models without random effects, making it a versatile choice for this problem.

# Load necessary libraries

LMMsetup <- function(form, dat, ref = list()) {
  # This function constructs the fixed-effects design matrix (X), random-effects 
  # design matrix (Z), and response variable (y) for a given linear mixed-effects model
  
  # Inputs :
  # form: A formula specifying the fixed effects part of the model e.g (y ~ x1 + x2).
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
  # browser()
  # Unique list from previous step
  unique_elements <- unique(unlist(ref))
  unique_list <- as.list(unique_elements)
  
  #convert the variables into factor variables if they are not
  # Iterate over each column name in the unique list
  for (col in unique_list) {
    # Check if the column is a factor
    if (!is.factor(dat[[col]])) {
         warning(paste("Column", col, "is not a factor variable"))
      }
    }

  # Construct random effects model matrix (Z)
  Z_list <- lapply(ref, function(vars) {
    # Create the model matrix
    model_mat <- model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "-1")), dat)
    
    # Get the dimensions of the matrix
    attr(model_mat, "columns") <- ncol(model_mat)
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
WBlock_mult_A <-function(small_block_Inv, A, sigma){
  # This function is auxiliary and helps to compute efficiently
  # the multiplication between a matrix 'Wblock' and 'A' (Wblock%*%A) where Wblock is of the form:
  
  # small_block_Inv    0
  #       0         I/(sigma^2)
  
  # where I is the identity matrix of n-p x n-p 
  
  # Inputs :
  # small_block_Inv: First diagonal block of pxp
  # sigma: constant dividing the second diagonal block (Identity matrix).
  
  # Returns : 
  # A matrix 'WBlock_by_A', the result of multiplying 'Wblock' and 'A' (Wblock %*% A).
  
  if (!is.matrix(A)) A<-as.matrix(A) #convert 'A' into a matrix if is a vector
  p<-ncol(small_block_Inv) 
  n<-nrow(A)
  WBlock_by_A_1_to_p <- small_block_Inv %*% A[1:p, , drop = FALSE]  # Multiply the small block inverse to the first p rows
  WBlock_by_A_p_1_to_n <- A[(p + 1):n, , drop = FALSE]/sigma^(2)    # Scale the remaining rows by 1 / sigma^2
  WBlock_by_A <- rbind(WBlock_by_A_1_to_p, WBlock_by_A_p_1_to_n)    # Combine weighted parts without further splitting
  return(WBlock_by_A)
}
LMMprof <- function(theta, setup) {
  
  #This function computes the negative -log-likelihood and Betas for a linear mixed model 
  # given a set of parameters (theta) and model matrices. It handles cases 
  # with or without random effects by formulating the likelihood function and using 
  # efficient matrix computations such as QR decomposition and Cholesky decomposition.
  # when no random effects provided then it computes the Betas and -log-likehood of a linear model
  # for linear mixed model:
  #   -log-likelihood = (1/2)*[(y-X*B_hat)^(t) * (W) * (y-X*B_hat) + log(det(R*Psi_b*R^(T)+I_p/sigma^(2))) + (n-p)log(sigma^2)]
  #    B_hat = (X^(t)*W*X)^(-1)*X^(T)*W*y
  #
  #             |( R*Psi_b*R^(t) + I_p*sigma^(2))^(-1)          0          |
  #   W =  Q *  |                0                         I_n_p/sigma^(2) | * Q^(t)
  #   Z= Q * | R |   
  #          | 0 |
  
  # for linear model:
  #    -log-likelihood = (1/2)*[(y-X*B_hat)^(t)*(y-X*B_hat)/sigma^(2) + log(det(I * sigma^2))]
  #    B_hat = R^(-1) * Q_f * y 
  #    X= Q * | R |   
  #           | 0 |
  
  
  
  #Inputs : 
  # theta : A numeric vector of parameters to be estimated:
  # The first element corresponds to the log of the residual variance.
  # Subsequent elements correspond to the log-transformed variances of the random effects.
  
  # setup : A named list of model components:
  #         X: Design Matrix 
  #         Z: Random effects Matrix
  #         y: Response variable Matrix
  #         qr_decomp:  qr decomposition of Z if random effects, qr decomposition of X if not random effects 
  #         R: qr.R(setup$qr_decomp) R matrix of the qr decomposition 
  
  # Returns : 
  # A scalar representing the negative log-likelihood for the given parameters. 
  # The estimated fixed effects (beta) are returned as an attribute "beta_hat". 
  
  #Extras de data
  X <- setup$X #Design Matrix
  Z <- setup$Z #Random effects Matrix
  y <- setup$y #Response variable Matrix
  qr_decomp<-setup$qr_decomp #qr_decom of Z if random effects
  R<-setup$R
  # Dimensions
  n <- length(y)
  p<-ncol(R)
  
  if(is.null(Z)){#if not random effects it estimates B and compute - log likelihood for a linear regression
    sigma <- exp(theta[1])
    beta_hat<-backsolve(R, qr.qty(qr_decomp, y)[1:p])# finds the betas (B) such that R %*% B =  Q_F %*% y
    residual <- y - (X %*% beta_hat) #Prediction error 
    # Minus log-likelihood calculation
    minus_log_likelihood <- 0.5 * ((t(residual)%*%residual/sigma^2)+log(sigma^(2*n)))
    attr(minus_log_likelihood, "beta_hat") <- beta_hat# creates an attribute with the betas 
    return(minus_log_likelihood)
    
    
  }else{ #if random effects it estimates beta and compute - log likelihood for random effects regression
    #Extract the length of each block in Z
    Z_cols<-setup$dimensions
    # Extract sigma and random effects variances from theta
    sigma <- exp(theta[1]) 
    psi_diag <- exp(2 * theta[-1])
    Psi_b <- diag(rep(psi_diag, Z_cols), ncol(Z), ncol(Z)) # Construct the covariance matrix Psi_b for the random effects
    
    # compute the first diagonal block of W ( R*Psi_b*R^(t) + I_p*sigma^(2))^(-1)
    small_block<-R%*%Psi_b%*%t(R)+diag(1, nrow = p, ncol = p)*(sigma^2)
    small_block_chol <- chol(small_block)# Cholesky decomposition of the first diagonal block
    small_block_Inv<-solve_chol(small_block_chol, diag(x = 1, nrow = p, ncol = p)) # Inverse the first diagonal block
    
    # Compute X^t*W*y
    Qty <- qr.qty(qr_decomp, y) # Compute Q^T y
    Wblock_Qty<-WBlock_mult_A(small_block_Inv=small_block_Inv, A=Qty, sigma=sigma)
    XtWy<-t(X) %*% qr.qy(qr_decomp, Wblock_Qty)
    
    # Compute X^t*W*X
    QtX <- qr.qty(qr_decomp, X)  # Compute Q^T X 
    Wblock_QtX<-WBlock_mult_A(small_block_Inv=small_block_Inv, A=QtX, sigma=sigma)
    XtWX <- t(X) %*% qr.qy(qr_decomp, Wblock_QtX)
    
    ######################################################
    # Compute Betas, residuals and negative log-likelihood 
    
    #B_hat computation solving XtWX*B=XtWy for B using de Cholesky decomposition
    XtWX_chol <- chol(XtWX) #Compute the Cholesky decomposition
    beta_hat<-solve_chol(L=XtWX_chol, b=XtWy) #solve for B 
    
    # residual
    residual <- y - (X %*% beta_hat)#compute the residual
    
    # computes efficiently res^(t)*W*rest  i.e (y-X*B_hat)^(t) * (W) * (y-X*B_hat)
    Qt_res <- qr.qty(qr_decomp, residual)
    Wblock_Qt_res<-WBlock_mult_A(small_block_Inv=small_block_Inv, A=Qt_res, sigma=sigma)
    restWrest<-t(residual)%*%qr.qy(qr_decomp, Wblock_Qt_res)
    
    #computes log(det(R*Psi_b*R^(T)+I_p/sigma^(2))) using the chol decomposition of A=R*Psi_b*R^(T)+I_p/sigma^(2)
    # A=S^(t)*S then log(det(A))= 2*sum(diag(S))
    det_small_block=2*sum(log(diag(small_block_chol)))
    
  
    # negative log-likelihood calculation
    minus_log_likelihood <- 0.5 * (restWrest + det_small_block + (n - p) * log(sigma^2))
    attr(minus_log_likelihood, "beta_hat") <- beta_hat
    return(minus_log_likelihood)
  }
}
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
  
  #If not random effects then apply the QR decomposition to X in order to compute a Linear regression
  if(is.null(setup$Z)){
    setup$qr_decomp <- qr(setup$X)
    setup$R <- qr.R(setup$qr_decomp)
    
    
  }else{#If there are random effects then apply the QR decomposition to Z 
    setup$qr_decomp <- qr(setup$Z)
    setup$R <- qr.R(setup$qr_decomp)
  }
  # Step 2: Initial guesses for theta: log(sigma) and log standard deviations for random effects
  theta_init <- rep(0, length(ref) + 1)  # Starting with a small positive value
  # Step 3: Optimize negative log-likelihood using `optim`
  opt <- optim(theta_init, LMMprof, setup = setup, method = "L-BFGS-B", control = list(fnscale = 1))
  #call the LMMprof with the thetas tha minimize the log likelihood 
  final_cost_value <- LMMprof(theta = opt$par, setup = setup)
  # Accessing the attribute "beta"
  beta_hat <- attr(final_cost_value, "beta_hat")
  return(list(beta = beta_hat, theta = opt$par))
}

