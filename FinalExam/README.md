1. Quantile regression model
----------------------------

    library(quantreg)

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    set.seed(12345)
    eps <- rnorm(50,0,2)
    x <- rchisq(50,df=2)-2
    y <- 2+x+eps

2. Maximum Likelihood
---------------------

    mle_func <- function(y,x){
      ones <- rep(1,50)
      x_mat <- matrix(c(ones,x),50,2)
      theta_mle <- solve(t(x_mat)%*%x_mat)%*%t(x_mat)%*%y
      err_var <- var(eps)
      inv_mat <- solve(t(x_mat)%*%x_mat)
      s <- sqrt((1/48)*sum((y-theta_mle[1]-theta_mle[2]*x_mat[,2])^2))
      ci_t0 <- theta_mle[1] + s*(sqrt(inv_mat[1,1]))*qt(c(.025, .975), df=48)
      ci_t1 <- theta_mle[2] + s*(sqrt(inv_mat[2,2]))*qt(c(.025, .975), df=48)
      
      return(list(parameters = theta_mle, ci = matrix(c(ci_t0,ci_t1),2,2)))
    }

    result <- mle_func(y,x)

    print("Fitted Parameters")

    ## [1] "Fitted Parameters"

    result$parameters

    ##           [,1]
    ## [1,] 2.3963213
    ## [2,] 0.8533437

    print("Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1")

    ## [1] "Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1"

    result$ci

    ##          [,1]      [,2]
    ## [1,] 1.765060 0.5111011
    ## [2,] 3.027583 1.1955863

3. Bootstrap
------------

    boot_t <- function(y,x){
      result <- rq(y~x,tau=0.5)
      residuals <- as.numeric(result$residuals)
      y_hat <- as.numeric(result$fitted.values)

      B = 2000
      n <- length(y_hat)
      theta_boot <- matrix(NA,B,2)
      
      for (i in 1:B){
        iz <- sample(1:n, n, replace=TRUE)
        yy <- y_hat+residuals[iz]
        result_boot <- rq(yy~x,tau=0.5)
        theta_boot[i,] <- result_boot$coefficients
      }
      
      alpha <- 0.05
      CI.t0 <- as.numeric(quantile(theta_boot[,1], c(alpha / 2, 1 - alpha / 2)))
      CI.t1 <- as.numeric(quantile(theta_boot[,2], c(alpha / 2, 1 - alpha / 2)))
      
      return(list(parameters = colMeans(theta_boot), ci = matrix(c(CI.t0,CI.t1),2,2)))
    }

    result <- boot_t(y,x)

    print("Fitted Parameters")

    ## [1] "Fitted Parameters"

    result$parameters

    ## [1] 2.7454520 0.7179819

    print("Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1")

    ## [1] "Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1"

    result$ci

    ##          [,1]     [,2]
    ## [1,] 1.955487 0.320771
    ## [2,] 3.293057 1.091556

4. Bayesian
-----------

    log.post <- function(t0,t1,X,y) {
      n <- length(X)
      return((-0.5*sum((y-t0-t1*X)*sign(y-t0-t1*X))))
    }

    rprop1 <- function(t0) rnorm(1,t0,1.1)
    rprop2 <- function(t1) rnorm(1,t1,0.55)

    mh.in.gibbs <- function(X,y,t0,t1,N){
      samples <- matrix(0,N,2)
      samples[1,] = c(t0,t1)
      ct1 <- ct2 <- 0
      for(i in 2:N){
        
        # sampling t0
        log.post.t0 <- function(theta0) log.post(theta0, samples[i-1,2], X, y)
        t0.star <- rprop1(samples[i-1,1])
        u <- runif(1)
        if(u <= min(1,exp(log.post.t0(t0.star)-log.post.t0(samples[i-1,1])))){
          samples[i,1] <- t0.star
          ct1<-ct1+1
        }else {
          samples[i,1] <- samples[i-1,1]
        }

        #sampling t1
        log.post.t1 <- function(theta1) log.post(samples[i,1], theta1, X, y)
        t1.star <- rprop2(samples[i-1,2])
        u <- runif(1)
        if(u <= min(1,exp(log.post.t1(t1.star)-log.post.t1(samples[i-1,2])))){
          samples[i,2] <- t1.star
          ct2<-ct2+1
        }else {
          samples[i,2] <- samples[i-1,2]
        }
      }
      
      alpha <- 0.05
      CI_t0 <- as.numeric(quantile(samples[,1], c(alpha / 2, 1 - alpha / 2)))
      CI_t1 <- as.numeric(quantile(samples[,2], c(alpha / 2, 1 - alpha / 2)))
      
      return(list(counts = c(ct1, ct2), samples = samples, ci = matrix(c(CI_t0,CI_t1),2,2)))
    }

    N<-10000
    result <- mh.in.gibbs(x, y, 1, 1, N)

    print("Acceptance rates")

    ## [1] "Acceptance rates"

    result$counts/N

    ## [1] 0.3279 0.3356

    print("Fitted Parameters")

    ## [1] "Fitted Parameters"

    colMeans(result$samples)

    ## [1] 2.7909689 0.7538453

    print("Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1")

    ## [1] "Confidence Intervals: Row 1: 2.5%, Row 2: 97.5%. Col 1: theta 0, Col 2: theta 1"

    result$ci

    ##          [,1]      [,2]
    ## [1,] 2.098357 0.4254873
    ## [2,] 3.375906 1.1217407

5. Monte Carlo Experiment
-------------------------

    set.seed(12345)
    x <- rchisq(50,df=2)-2
    set.seed(Sys.time())

    M = 100
    coverage_stat <- matrix(NA,M,6)
    length_stat <- matrix(NA,M,6)
    for (i in 1:M){
      theta_stat <- matrix(NA,6,2)
      eps <- rnorm(50,0,2)
      yy <- 2+x+eps
      result_mle <- mle_func(yy,x)
      theta_stat[1:2,] <- result_mle$ci
      result_boot <- boot_t(yy,x)
      theta_stat[3:4,] <- result_boot$ci
      result_gibbs <- mh.in.gibbs(x,yy, 1, 1, N)
      theta_stat[5:6,] <- result_gibbs$ci
      
      for (j in c(1,3,5)){
        coverage_stat[i,j] <- ifelse((theta_stat[j,1]<=2 & 2 <=theta_stat[j+1,1]),1,0)
        coverage_stat[i,j+1] <- ifelse((theta_stat[j,2]<=1 & 1<=theta_stat[j+1,2]),1,0)
      }
      
      for (l in c(1,3,5)){
        length_stat[i,l] <- theta_stat[l+1,1]-theta_stat[l,1]
        length_stat[i,l+1] <- theta_stat[l+1,2]-theta_stat[l,2]
      }
      
    }

    result_MCMC <- data.frame(t(colMeans(coverage_stat)))
    colnames(result_MCMC) <- c("theta_0_MLE","theta_1_MLE",
                               "theta_0_Boot","theta_1_Boot",
                               "theta_0_Bayesian","theta_1_Bayesian")

    avg_length <- data.frame(t(colMeans(length_stat)))
    colnames(avg_length) <- c("theta_0_MLE_AverageLength","theta_1_MLE_AverageLength",
                               "theta_0_Boot_AverageLength","theta_1_Boot_AverageLength",
                               "theta_0_Bayesian_AverageLength","theta_1_Bayesian_AverageLength")

    print(result_MCMC)

    ##   theta_0_MLE theta_1_MLE theta_0_Boot theta_1_Boot theta_0_Bayesian
    ## 1        0.95        0.97         0.94         0.92             0.93
    ##   theta_1_Bayesian
    ## 1             0.99

    print(avg_length)

    ##   theta_0_MLE_AverageLength theta_1_MLE_AverageLength
    ## 1                   1.16372                 0.6271224
    ##   theta_0_Boot_AverageLength theta_1_Boot_AverageLength
    ## 1                   1.362357                  0.7431997
    ##   theta_0_Bayesian_AverageLength theta_1_Bayesian_AverageLength
    ## 1                       1.272106                      0.6989108

Accroding to the results, MLE method is the most accurate and the most
of various methods achieve a coverage probability of 95%.
