Topics
------

Poisson regression model, Newton’s Method update for the MLE estimator,
Newton’s Method with the log link, Weibull distribution, MLEs using the
fixed point iteration, Frechet distirbution, “rank-one update”
quasi-Newton Method, “BFGS” method

Problems 1:
-----------

![](images/p1.png)

### 1.1

![](images/p1_1.png)

    library(glm2)
    data(crabs)

    crabs$Dark <- ifelse(crabs$Dark == "no",0,1)
    crabs$GoodSpine <- ifelse(crabs$GoodSpine == "no",0,1)

    poisson_loglikelihood <- function(theta, data){
      
      y <- data[1]
      x1 <- data[2]
      x2 <- data[3]
      x3 <- data[4]
      theta_1 <- theta[1]
      theta_2 <- theta[2]
      theta_3 <- theta[3]
      
      return( 
        sum(y*(x1*theta_1+x2*theta_2+x3*theta_3)-exp(x1*theta_1+x2*theta_2+x3*theta_3)-log(factorial(y)))
      )
    }

    gradient_poisson_loglikelihood <- function(theta, data){
      
      y <- data[1]
      x1 <- data[2]
      x2 <- data[3]
      x3 <- data[4]
      theta_1 <- theta[1]
      theta_2 <- theta[2]
      theta_3 <- theta[3]
      
      return( 
        
        matrix( 
          c(
            sum((y*x1)-(x1*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
            sum((y*x2)-(x2*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
            sum((y*x3)-(x3*exp(x1*theta_1+x2*theta_2+x3*theta_3)))
          )
          , 3,1))
    }

    hessian_poisson_loglikelihood <- function(theta, data){
      
      y <- data[1]
      x1 <- data[2]
      x2 <- data[3]
      x3 <- data[4]
      theta_1 <- theta[1]
      theta_2 <- theta[2]
      theta_3 <- theta[3]
      
      return( 
        
        matrix( c(  
          
          sum(-((x1^2)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x1^x2)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x1^x3)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x2*x1)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x2^2)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x2*x3)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x3*x1)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x3*x2)*exp(x1*theta_1+x2*theta_2+x3*theta_3))),
          sum(-((x3^2)*exp(x1*theta_1+x2*theta_2+x3*theta_3)))
          
          
        ) , 3, 3 )
        
      )
    }

    inv_grad_loglikelihood <- function(theta, data) {
      
      return(
        solve(hessian_poisson_loglikelihood(theta, data))
      )
    }

    newtons_method <- function(init_par, fun, grad, inv_hess, tol, max_iter, data){
      
      rel_tol <- 2*tol
      par_old <- init_par
      iter_count <- 0
      iter_par <- matrix(0,max_iter+1,length(init_par))
      iter_par[1,] <- t(init_par)   
      
      while(rel_tol > tol & iter_count < max_iter){
        
        par_new <- par_old - inv_hess(par_old, data)%*%grad(par_old, data)
        rel_tol <- max(abs(par_new - par_old)/par_old)
        par_old <- par_new
        iter_count <- iter_count + 1
        iter_par[iter_count+1,] <- t(par_new)
        
      }
      
      return(list(solution = par_new, fun_solution = fun(t(par_new), data), 
                  final_tol = rel_tol, num_iters = iter_count, all_iters = iter_par[1:iter_count,]))
      
    }


    NM_poisson <- newtons_method(c(0.1,0.1,0.1), poisson_loglikelihood, 
                                 gradient_poisson_loglikelihood, inv_grad_loglikelihood, 0.001, 10, crabs)

    NM_poisson

    ## $solution
    ##             [,1]
    ## [1,]  0.04665633
    ## [2,] -0.42138494
    ## [3,] -0.02511220
    ## 
    ## $fun_solution
    ## [1] -470.2652
    ## 
    ## $final_tol
    ## [1] 0.0005894555
    ## 
    ## $num_iters
    ## [1] 5
    ## 
    ## $all_iters
    ##            [,1]        [,2]        [,3]
    ## [1,] 0.10000000  0.10000000  0.10000000
    ## [2,] 0.06902628  0.08809015  0.19995135
    ## [3,] 0.05032395 -0.06592755  0.17245920
    ## [4,] 0.04627537 -0.31262213  0.03621928
    ## [5,] 0.04668385 -0.41741060 -0.02456145

    #calculate without intercept
    glm(Satellites~Width+Dark+GoodSpine-1, family=poisson(link="log"), data = crabs) #for comparison

    ## 
    ## Call:  glm(formula = Satellites ~ Width + Dark + GoodSpine - 1, family = poisson(link = "log"), 
    ##     data = crabs)
    ## 
    ## Coefficients:
    ##     Width       Dark  GoodSpine  
    ##   0.04666   -0.42152   -0.02523  
    ## 
    ## Degrees of Freedom: 173 Total (i.e. Null);  170 Residual
    ## Null Deviance:       1051 
    ## Residual Deviance: 585.2     AIC: 946.5

Problems 2:
-----------

![](images/p2.png)

### 2.1

![](images/p2_1.png)

    library(evd)
    set.seed(12345)
    wb <- rweibull(50, shape=5, scale = 7)

    beta_theta_mle <- function(beta_theta, x){
      
      n <- length(x)
      
      return(
        c(((sum((x^beta_theta[1])*log(x))/sum((x^beta_theta[1]))) - (1/n)*sum(log(x)))^-1,((1/n)*sum((x^beta_theta[2])))^(1/beta_theta[2]))
        
      )
    }

    b_t_last <- c(.1,.1)
    b_t_next <- c(1,1)

    while(abs(b_t_next[1] - b_t_last[1])>0.01 | abs(b_t_next[2] - b_t_last[2])>0.01){
      b_t_last[1] <- b_t_next[1]
      b_t_last[2] <- b_t_next[2]
      b_t_next <- beta_theta_mle(b_t_last,wb)
      
    }
    print("estimated beta and theta are:")

    ## [1] "estimated beta and theta are:"

    print(b_t_next)

    ## [1] 4.634475 7.485158

Problems 3:
-----------

![](images/p3.png)

### a). Use the Newton's Method algorithm (Secant Method as required by professor) from the lecture notes to find the MLEs.

    library(evd)
    set.seed(12345)
    x = rfrechet(100, shape=3, scale=2)

    frechet_loglikelihood <- function(k_theta, data){
      
      a <- k_theta[1]
      s <- k_theta[2]
      n <- length(data)
      
      return( 
        n*log(a)-n*log(s)-(1+a)*sum(log(data))+(n*(1+a)*log(s))-sum((s/data)^a)
      )
    }

    secant_method <- function(init_par, fun, tol, max_iter, data){
      
      eps <- 0.01
      rel_tol <- 2*tol
      par_old <- init_par
      iter_count <- 0
      iter_par <- matrix(0,max_iter+1,length(init_par))
      iter_par[1,] <- t(init_par)
      dim_par <- length(init_par)

      approx_grad <- rep(0, dim_par)
      approx_hess <- matrix(0, dim_par, dim_par)

      
      while(rel_tol > tol & iter_count < max_iter){
        
        for(i in 1:dim_par){
          del <- rep(0,dim_par)
          del[i] <- eps
          approx_grad[i] <- (0.5/eps)*(fun(t(par_old+matrix(del,dim_par,1)),data) - fun(t(par_old), data))
        }

        for(i in 1:dim_par){
          del1 <- rep(0,dim_par)
          del1[i] <- eps
          for(j in 1:dim_par){
            del2 <- rep(0,dim_par)
            del2[j] <- eps

            if(i == j){

              approx_hess[i,j] <- (1/(eps^2))*(fun(t(par_old
                                                     +matrix(del1,dim_par,1)),data)
                                               -2*fun(t(par_old),data) 
                                               +fun(t(par_old-matrix(del1,dim_par,1)),data))    
            }else {
              approx_hess[i,j] <- (0.25/(eps^2))*(fun(t(par_old+matrix(del1,dim_par,1)
                                                        +matrix(del2,dim_par,1)),data)
                                                  -fun(t(par_old+matrix(del1,dim_par,1)
                                                         -matrix(del2,dim_par,1)),data)
                                                  -fun(t(par_old-matrix(del1,dim_par,1)
                                                         +matrix(del2,dim_par,1)),data)
                                                  +fun(t(par_old-matrix(del1,dim_par,1)
                                                         -matrix(del2,dim_par,1)),data))
            }
          }
        }
        
        
        par_new <- par_old - solve(approx_hess)%*%approx_grad 
        rel_tol <- max(abs(par_new - par_old)/par_old)
        par_old <- par_new
        iter_count <- iter_count + 1
        iter_par[iter_count+1,] <- t(par_new)
        
      }
      
      return(list(solution = par_new, fun_solution = fun(t(par_new), data), 
                  final_tol = rel_tol, num_iters = iter_count, all_iters = iter_par[1:iter_count,]))
    }


    se_frechet <- secant_method(matrix(c(1,1),2,1), frechet_loglikelihood, 0.0000001, 5000, x)
    se_frechet

    ## $solution
    ##          [,1]
    ## [1,] 2.829017
    ## [2,] 2.013720
    ## 
    ## $fun_solution
    ## [1] -142.486
    ## 
    ## $final_tol
    ## [1] 7.824775e-08
    ## 
    ## $num_iters
    ## [1] 24
    ## 
    ## $all_iters
    ##           [,1]     [,2]
    ##  [1,] 1.000000 1.000000
    ##  [2,] 1.883128 2.059671
    ##  [3,] 2.233805 2.075944
    ##  [4,] 2.483284 2.051749
    ##  [5,] 2.640110 2.034872
    ##  [6,] 2.729628 2.025061
    ##  [7,] 2.777831 2.019698
    ##  [8,] 2.802955 2.016845
    ##  [9,] 2.815825 2.015347
    ## [10,] 2.822360 2.014565
    ## [11,] 2.825662 2.014158
    ## [12,] 2.827328 2.013947
    ## [13,] 2.828167 2.013838
    ## [14,] 2.828589 2.013781
    ## [15,] 2.828802 2.013751
    ## [16,] 2.828909 2.013736
    ## [17,] 2.828962 2.013728
    ## [18,] 2.828989 2.013724
    ## [19,] 2.829003 2.013722
    ## [20,] 2.829010 2.013721
    ## [21,] 2.829013 2.013721
    ## [22,] 2.829015 2.013720
    ## [23,] 2.829016 2.013720
    ## [24,] 2.829016 2.013720

Problems 3:
-----------

### b). Write a "rank-one update" quasi-Newton Method in R and use it to find the MLEs.

    library(evd)
    set.seed(12345)
    r_frechet = rfrechet(100, shape=3, scale=2)
    frechet_loglikelihood <- function(k_theta,x){
      
      a <- k_theta[1]
      s <- k_theta[2]
      n <- 100
      
      return( 
        n*log(a)-n*log(s)-(1+a)*sum(log(x))+(n*(1+a)*log(s))-sum((s/x)^a)
      )
    }


    gradient_frechet_loglikelihood <- function(k_theta,x){
      
      a <- k_theta[1]
      s <- k_theta[2]
      n <- 100
      
      return( 
        matrix(
          
          c( (n/a)-sum(log(x))+(n*log(s))-sum(((s/x)^a)*log(s/x)), 
             
             (-n/s)+(n*(1+a)/s)-(a/s)*sum((s/x)^a)
             
          ),2,1))
    }


    rank1_newtons_method <- function(init_par, fun, grad, tol, max_iter, data){
      
      rel_tol <- 2*tol
      par_old <- init_par
      iter_count <- 0
      iter_par <- matrix(0,max_iter+1,length(init_par))
      iter_par[1,] <- t(init_par)
      mu = 0.01
      H <- diag(2)
      
      while(rel_tol > tol & iter_count < max_iter){
        
        par_new <- par_old - mu*H%*%grad(par_old, data)
        
        if (fun(par_new,data) <= fun(par_old,data)){
          
          mu <- mu/2
        }
        
        z <- par_new-par_old
        y <- grad(par_new, data) - grad(par_old, data)
        add_term <- ((z-H%*%y)%*%(t(z-H%*%y)))/as.vector((t(z-H%*%y))%*%y)
        H <- H + add_term
        
        rel_tol <- max(abs(par_new - par_old)/par_old)
        par_old <- par_new
        iter_count <- iter_count + 1
        iter_par[iter_count+1,] <- t(par_new)
        
      }
      
      return(list(solution = par_new, fun_solution = fun(t(par_new), data),
                  final_tol = rel_tol, num_iters = iter_count))
      
    } #removed all_iters

    rank1_NM_rfrechet <- rank1_newtons_method(c(1,1), frechet_loglikelihood, 
                                              gradient_frechet_loglikelihood, 0.0000001, 5000, r_frechet)

    rank1_NM_rfrechet

    ## $solution
    ##          [,1]
    ## [1,] 2.829473
    ## [2,] 2.018691
    ## 
    ## $fun_solution
    ## [1] -142.4836
    ## 
    ## $final_tol
    ## [1] 9.975519e-08
    ## 
    ## $num_iters
    ## [1] 4694

Problems 3:
-----------

### c). Use R's optim function with the "BFGS" method (this is a rank-two update) to find the MLEs.

    library(evd)
    set.seed(12345)
    x = rfrechet(100, shape=3, scale=2)
    frechet_loglikelihood <- function(k_theta){
      
      a <- k_theta[1]
      s <- k_theta[2]
      n <- 100
      
      return( 
        -(n*log(a)-n*log(s)-(1+a)*sum(log(x))+(n*(1+a)*log(s))-sum((s/x)^a))
      )
    }


    gradient_frechet_loglikelihood <- function(k_theta){
      
      a <- k_theta[1]
      s <- k_theta[2]
      n <- 100
      
      return( 
        
        c( -((n/a)-sum(log(x))+(n*log(s))-sum(((s/x)^a)*log(s/x))), 
                    
           -((-n/s)+(n*(1+a)/s)-(a/s)*sum((s/x)^a))
        
      ))
    }

    rank2 <- optim(c(1, 1), fn=frechet_loglikelihood, gr = gradient_frechet_loglikelihood, 
                 method = "BFGS",control=list(trace = 1, reltol = 0.0000001))

    ## initial  value 224.477942

    ## Warning in log(s): NaNs produced

    ## Warning in log(s): NaNs produced

    ## Warning in log(s): NaNs produced

    ## Warning in log(s): NaNs produced

    ## Warning in log(s): NaNs produced

    ## Warning in log(s): NaNs produced

    ## iter  10 value 142.483718
    ## final  value 142.483581 
    ## converged

    rank2

    ## $par
    ## [1] 2.829554 2.018701
    ## 
    ## $value
    ## [1] 142.4836
    ## 
    ## $counts
    ## function gradient 
    ##       27       11 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

Problems 3:
-----------

### d). Compare the results. Are there any significant differences in number of iterations for convergence?

Under the same tolerance level at 0.0000001, both the Secant Method and
Rank-two Update Quasi-Newton Method appear to have the small number of
iterations (24 and 11) to achieve convergence. However, it takes around
4600 iterations for our version of Rank-one Update Quasi-Newton Method
to make convergence. I think the reason why this happens is that in
order to avoid overshooting and exploding gradients problems, we need a
very small step size mu at 0.01 and it will be updated to be the half of
itself for some iteration which new parameters does not increase the
likelihood of estimator, so that the mu will become even smaller and it
will need more iterations to make convergence.
