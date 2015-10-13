#-----------------------------------------------------------------------------#
#                                                                             # 
#                  Some Models for Process Variation                          #
#                                                                             #
# Reference:                                                                  #
# Wetherill, G.B. and Brown, D.W. (1991) ``Statistical Process Control'',     #
#   New York, Chapman and Hall, Chapter 3                                     #
#-----------------------------------------------------------------------------#

qcc.demo.menu <- function()
{

  title <- paste("#-----------------------------------------#",
                 "# Some models for process variation.      #",
                 "# Make a selection (or 0 to exit):        #",
                 "#-----------------------------------------#", 
                 sep="\n")
  choice <- 100
  while(choice != 0)
       { choice <- menu(c("Simple random variation", "Extra variation",
                          "Autocorrelation", "Recurring cycles",
                          "Trends", "Mixture", "Sudden jumps"),
                        title = title)
          switch(choice,
                 qcc.demo.random(),
                 qcc.demo.extravar(),
                 qcc.demo.autocor(),
                 qcc.demo.cycles(),
                 qcc.demo.trends(),
                 qcc.demo.mixture(),
                 qcc.demo.jumps()) 
       }
  invisible()
}

pause <- function () 
{ cat("Pause. Press <Enter> to continue...")
  readline()
  invisible()
}

qcc.demo.random <- function()
{
  cat("# 1) Simple random variation \n")
  cat("#    x_ij = mu + sigma_W * epsilon_ij \n")
  mu <- 100
  sigma.W <- 10
  epsilon <- rnorm(500)
  x <- matrix(mu + sigma.W*epsilon, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

qcc.demo.extravar <- function()
{
  cat("# 2) Between and Within group variation: extra variation \n")
  cat("#    x_ij = mu + sigma_B * u_i + sigma_W * epsilon_ij \n")
  mu <- 100
  sigma.W <- 10
  sigma.B <- 5
  epsilon <- rnorm(500)
  u <- as.vector(sapply(rnorm(50), rep, 10))
  x <- mu + sigma.B*u + sigma.W*epsilon
  x <- matrix(x, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

qcc.demo.autocor <- function()
{
  cat("# 3) Simple autocorrelation model \n")
  cat("#    x_ij = mu + W_i + sigma_W * epsilon_ij \n")
  cat("#     W_i = rho * W_i-1 + sigma_B * u_i =  \n")
  cat("#         = sigma_B * u_i + rho * sigma_B * u_i-1 + rho^2 * sigma_B * u_i-2 + .... \n")
  cat("#    (W_0=0) \n")
  mu <- 100
  rho <- 0.6
  sigma.W <- 10
  sigma.B <- 5
  epsilon <- rnorm(500)
  u <- rnorm(500)
  W <- rep(0,100)
  for (i in 2:length(W))
      W[i] <- rho*W[i-1] + sigma.B*u[i]
  x <- mu + sigma.B*u + sigma.W*epsilon
  x <- matrix(x, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

qcc.demo.cycles <- function()
{
  cat("# 4) Recurring cycles \n")
  cat("#    Assume we have 3 working turns of 8 hours each for each \n")
  cat("#    working day, so 8*3=24 points in time, and at each point \n")
  cat("#    we sample 5 units. \n")
  cat("#            x_ij = mu + W_i + sigma_W * epsilon_ij  \n")
  cat("#    where W_i (i=1,...,8) is the cycle \n")
  mu <- 100
  sigma.W <- 10
  epsilon <- rnorm(120, sd=0.3)
  W <- c(-4, 0, 1, 2, 4, 2, 0, -2) # assumed workers cycle
  W <- rep(rep(W, rep(5,8)), 3)
  x <- mu + W + sigma.W*epsilon
  x <- matrix(x, ncol=5, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

qcc.demo.trends <- function()
{
  cat("# 4) Trends \n")
  cat("#            x_ij = mu + W_i + sigma_W * epsilon_ij \n")
  cat("#    where W_i = 0.2 * i  \n")
  mu <- 100
  sigma.W <- 10
  epsilon <- rnorm(500)
  W <- rep(0.2*1:100, rep(5,100))
  x <- mu + W + sigma.W*epsilon
  x <- matrix(x, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

qcc.demo.mixture <- function()
{
  cat("# 5) Mixture \n")
  cat("#          x_ij = mu1 * p + mu2 * (1-p) + sigma_W * epsilon_ij \n")
  cat("#    where p = Pr(Process #1) \n")

  mu1 <- 90
  mu2 <- 110
  sigma.W <- 10
  epsilon <- rnorm(500)
  p <- rbinom(50, 1, 0.5)
  mu <- mu1*p + mu2*(1-p)
  x <- rep(mu, rep(10, length(mu))) + sigma.W*epsilon
  x <- matrix(x, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}
  
qcc.demo.jumps <- function()
{
  cat("# 6) Sudden jumps  \n")
  cat("#        x_ij = mu_i + sigma_W * epsilon_ij  \n")
  cat("# where mu_i is the mean of the process for state i (i=1,...,k) \n")
  mu <- rep(c(95,110,100,90), c(20,35,25,20))
  sigma.W <- 10
  epsilon <- rnorm(500)
  x <- rep(mu, rep(5, length(mu))) + sigma.W*epsilon
  x <- matrix(x, ncol=10, byrow=TRUE)
  pause()
  print(qcc(x, type="xbar"))
  pause()
  print(qcc(x, type="R"))
  pause()
  print(qcc(x, type="S"))
  pause()
  invisible(x)
}

# Run demo
qcc.demo.menu()
