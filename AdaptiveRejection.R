#This is the final file for the adaptive rejection sampling algorithm

###################### initG ######################
initG <- function(f, upperB = Inf, lowerB = -Inf) {
  ### initialize the bounding function g of the distribution of interest f
  ### picks -5 and 5 as starting values, adjusts if nessecary
  if(upperB == Inf || lowerB == -Inf) {
    #If either bound is infinite, use 5 or -5 as starting point for upper or lower bounds
    u <- (max(lowerB, -6)+1)
    v <- (min(upperB, 6)-1)    
  } else {
    #If neither bound is infinite, divide interval into thirds
    u <- (2 * lowerB + upperB) / 3
    v <- (2 * upperB + lowerB) / 3
  }
  # browser()
  #Use the fPrime module to evaluate slopes at u and v
  logf <- function(x) { log(f(x)) }
  uSlope <- grad(logf, x=u, method='simple')
  vSlope <- grad(logf, x=v, method='simple')
  
  #If f'(u) <= 0, slide u towards lowerB until f'(u) > 0
  while(uSlope <= 0) {
    u <- (2 * u + max(-abs(3*u), lowerB)) / 3
    uSlope <- grad(f, x=u, method='simple')
  }
  #If f'(v) >= 0, slide v towards upperB until f'(v) < 0
  while(vSlope >= 0) {
    v <- (2 * v + min(abs(3*v), upperB)) / 3
    vSlope <- grad(f, x=v, method='simple')
  }
  uVal <- log(f(u))
  vVal <- log(f(v))
  
  #Find the intersection point of the tangent lines
  w <- (log(f(v)) - log(f(u)) - v*vSlope + u*uSlope)/(uSlope - vSlope)
  
  #Create g
  g <- rbind(c(lowerB, w, u, uSlope, uVal), c(w, upperB, v, vSlope, vVal))
  colnames(g) <- c('start', 'end', 'intersect', 'm', 'b')
  gLower <- matrix(updateGLower(0, g, f), nrow = 1, ncol = 5)
  colnames(gLower) = c('start', 'end', 'intersect', 'm', 'b')
  return(list(Upper=g, Lower=gLower))
}
############################################

###################### findIntersect ######################
findIntersect <- function(m1, x1, b1, m2, x2, b2){
  ### return x value of the intersection of the lines:
  ### y = m1(x - x1) + b1 and y = m2(x - x2) + b2
  if (m1 != m2) {
    intersect <- (b2 - b1 - m2*x2 + m1*x1)/(m1 - m2)
    return(intersect)
  } else {
    return(NA)
  }
}

###################### updateG ######################
updateGUpper <- function(x, g, f, fval) {
  ### Return g with an added row with intersect = x
  ### Calculates which other rows of g need updating and updates them
  ### 
  ### x: The point to update at
  ### f: The distribution being drawn from
  ### g: A matrix where each row's values are start, end, intersect, m and b
  ###   where start and end define the range g applies to, intersect is the
  ###   point x at which g is tangent to log(f(x)), and b is the value log(f(x))
  
  # find index of the function whose range includes x:
  toUpdate <- which(g[ ,'start'] <= x & g[ ,'end'] > x)
  logfval <- log(fval)
  logf <- function(x) { log(f(x)) }
  fprime <- grad(logf, x=x)
  
  # check if x is to the left or right of the intersection toUpdate
  # update either the g to the right or left.
  if ((toUpdate == 1) & (x < g[1, 'intersect'])) {
    newRangeLeft <- g[1, 'start']
    newRangeRight <- findIntersect(fprime, x, logfval, 
                                   g[1, 'm'], g[1, 'intersect'], g[1, 'b'])
    g[1, 'start'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  } else if ((toUpdate == nrow(g)) & (x > g[nrow(g), 'intersect'])) {
    newRangeRight <- g[nrow(g), 'end']
    newRangeLeft <- findIntersect(fprime, x, logfval, 
                                  g[nrow(g), 'm'], g[nrow(g), 'intersect'], 
                                  g[nrow(g), 'b'])
    g[nrow(g), 'end'] <- newRangeLeft
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  } else if (x < g[toUpdate, 'intersect']) {
    left <- toUpdate - 1
    newRangeLeft <- findIntersect(fprime, x, logfval, 
                                  g[left, 'm'], g[left, 'intersect'], g[left, 'b'])
    newRangeRight <- findIntersect(fprime, x, logfval, 
                                   g[toUpdate, 'm'], g[toUpdate, 'intersect'], 
                                   g[toUpdate, 'b'])
    g[left, 'end'] <- newRangeLeft
    g[toUpdate, 'start'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  } else if (x > g[toUpdate, 'intersect']) {
    right <- toUpdate + 1
    newRangeRight <- findIntersect(fprime, x, logfval, 
                                   g[right, 'm'], g[right, 'intersect'], 
                                   g[right, 'b'])
    newRangeLeft <- findIntersect(fprime, x, logfval, 
                                  g[toUpdate, 'm'], g[toUpdate, 'intersect'], 
                                  g[toUpdate, 'b'])
    g[right, 'start'] <- newRangeRight
    g[toUpdate, 'end'] <- newRangeLeft
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  }
  g <- g[sort.list(g[ ,1], ), ]
  return(g)
}

updateGLower <- function(x, gu, f) {
  ### Return the lower bound glower, based on gu, the upper bound
  ### Returns a matrix where each row is contain start, end, intersect, m, and b
  ### start and end define the range that g lower is valid on
  ### intersect is not used
  ### m is the slope of the line
  ### b is the value of log(f(x)) at x = 'start', which can be used to
  ###   calculate the equation of the line.
  
  glower <- gu
  glower[ , 'start'] <- gu[ , 'intersect']
  glower[ , 'end'] <- c(gu[-1, 'intersect'], Inf)
  # calculate slope:
  glower[ , 'm'] <- c(diff(glower[ , 'b']), 0)/
                     (glower[ , 'end'] - glower[ , 'start'])
  glower <- glower[-nrow(glower), ] # remove last row - it's meaningless
  glower # note the b column is the value of log(f) at 'start'
}

checkConcav <- function(fval, g_upperVal, glist){
  ### check for log concavity of f by checking three cases
  ### a non-log concave function would fail
  if(fval > g_upperVal) {
    # make sure g is bounding f
    stop("f is not log-concave; this method will not work.")
  } else if (any(diff(glist$Upper[ , 'end']) < 0)) {
    # make sure endpoints always increasing
    stop("f is not log-concave; this method will not work.")
  } else if (any(diff(glist$Upper[ , 'm']) > 0)) {
    # make sure slope always decreasing
    stop("f is not log-concave; this method will not work.")
  }
}

updateG <- function(x, glist, f){
  ### Return a list with elements Upper (the upper bound of g in matrix form)
  ### and Lower (the lower bound of g in matrix form)
  
  index_Upper <- which(glist$Upper[ ,'start'] <= x & glist$Upper[ ,'end'] > x)
  # upperX is the value of the upper bound of f
  # calculated here to check that g bounds f
  upperX <- exp(glist$Upper[index_Upper, 'm'] * 
                (x - glist$Upper[index_Upper, 'intersect']) + 
                glist$Upper[index_Upper, 'b'])
  # f is evaluted here as part of checking concavity
  # this value is passed to updateG to avoid recalculation
  fval <- f(x)
  checkConcav(fval, upperX, glist)
  gu <- updateGUpper(x, glist$Upper, f, fval)
  gLower <- updateGLower(x, gu, f)
  return(list(Upper=gu, Lower=gLower, fx=fval))
}
############################################

###################### generatePoints ######################
generatePoints <- function(N, glist, f){
  ### Attempt to sample N points from f using rejection sampling based on
  ### glist$Upper and glist$Lower, which bound f. Update glist and exit 
  ### if evaluation of f is required.
  X <- sampleSX(glist$Upper, N)
  # sample N points from upper bound of g function
  U <- runif(N)
  # generate N points from uniform (0,1) distribution
  sampleX <- NULL
  # a vector to store X's that can be accepted
  # print(glist)
  for (i in 1:N){
    if(X[i] < glist$Lower[1,1] || X[i] > glist$Lower[nrow(glist$Lower),2]) {
      index_Upper <- which(glist$Upper[ ,'start'] <= X[i] & glist$Upper[ ,'end'] > X[i])
      upperX <- glist$Upper[index_Upper, 'm'] * (X[i] - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
      glist <- updateG(X[i], glist, f)
      if (U[i] < glist$fx / exp(upperX)){
        sampleX <- c(sampleX, X[i])
        break
      } else {
        break
      } 
    } else {
      index_Upper <- which(glist$Upper[ ,'start'] <= X[i] & glist$Upper[ ,'end'] > X[i])
      index_Lower <- which(glist$Lower[ ,'start'] <= X[i] & glist$Lower[ ,'end'] > X[i])
      # in fact index_Upper and index_Lower have some relation
      upperX <- glist$Upper[index_Upper, 'm'] * (X[i] - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
      # value of upper bound at X[i]
      lowerX <- glist$Lower[index_Lower, 'm'] * (X[i] - glist$Lower[index_Lower, 'start']) + glist$Upper[index_Lower, 'b']
      # value of lower bound at X[i]
      if (U[i] < exp(lowerX - upperX)){
        sampleX <- c(sampleX, X[i])
      }
      else{
        glist <- updateG(X[i], glist, f)
        if (U[i] < glist$fx / exp(upperX)){
          sampleX <- c(sampleX, X[i])
          break
        }
        else{
          break
        } 
      }    
    }
  }
  return(list("sample" = sampleX, "g" = glist))   
}
############################################

###################### sampleSX ######################
sampleSX <- function(g,n) {
  #this function draws n random samples from the exponentiated upper envelope
  #xk is a vector of the k abscissae x_i
  #gk is a vector of g(x)=log(f(x)) the k abscissae x_i
  #gpx is a vector of g'(x) of the k abscissae x_i
  #zk is the intersection points of the tangents
  xk <- g[,'intersect']
  k <- length(xk)
  gk <- g[,'b']
  gPrime <- g[,'m']
  zk <- c(g[,'start'], g[k,'end']) #note that there are (k+1) z values
  
  #to calculate the value of the upper envelope at the z points, interpolate using the tangent lines:
  zkMinusxk <- zk - c(0, xk)
  guz <- gPrime*(zkMinusxk[-1]) + gk
  #the first value (at zk[1]) is calculated manually:
  guz0 <- guz[1] - gPrime[1]*(zk[2] - zk[1])
  guz <- c(guz0, guz)
  
  #here it is, a vectorized version of computing the areas:
  t <- exp(c(guz,1)) - exp(c(1,guz))
  areas <- t[-c(1, k+2)]/gPrime
  
  #correcting for values of gPrime approximately equal to 0
  gPrimeZero <- which(sapply(gPrime, function(x){ all.equal(x,0) } ) == TRUE)
  areas[gPrimeZero] <- exp(gk[gPrimeZero])*(zk[gPrimeZero+1] - zk[gPrimeZero])
  
  #normalizing factor of the exponentiated upper hull:
  normFactor <- sum(areas)
  scum <- c(0, cumsum(areas)/normFactor)
  
  #generate uniforms to plug into inverse cdf
  u <- runif(n)  
  whichChunk <- sapply(u, function(x) { which(x < scum)[1]-1 })
  
  #now for the inverse cdf, broken up into a few pieces for readability
  piece1 <- (gPrime[whichChunk])*normFactor*(u - scum[whichChunk])
  piece2 <- log(exp(guz[whichChunk]) + piece1) - guz[whichChunk]
  sample <- zk[whichChunk] + piece2/gPrime[whichChunk]
  
  #correcting the NaN's for samples with the uniform in the first chunk
  firstChunk <- which(whichChunk == 1)
  piece11 <- gPrime[1]*normFactor*(scum[2] - u[firstChunk])
  piece12 <- log(exp(guz[2]) - piece11) - guz[2]
  sample[firstChunk] <- zk[2] + piece12/gPrime[1]
  
  #correcting for values of gPrime equal to 0
  for (b in gPrimeZero) {
    chunkgPrimeZero <- which(whichChunk==b)
    sample[chunkgPrimeZero] <- zk[b] + (u[chunkgPrimeZero] - scum[b]) * 
                               normFactor*(zk[b+1] - zk[b])
  }
  return(sample)
}
############################################

###################### plotg ######################
plotg <- function(f, glist){
  # plot log(f) and the upper and lower bounds of g
  x2 <- seq(-6, 6, by=.01)
  plot(x2, log(f(x2)), type='l')
  g <- glist$Upper
  for (i in 1:nrow(g)){
    x1 <- x2[x2 < g[i, 'end'] & x2 > g[i, 'start']]
    lines(x1, g[i , 'm']*(x1 - g[i , 'intersect']) + g[i , 'b'], col='red')
    abline(v=g[i, 'end'], col='blue')
  }
  glow <- glist$Lower
  for (i in 1:nrow(glow)) {
    x1 <- x2[x2 < glow[i, 'end'] & x2 > glow[i, 'start']]
    lines(x1, glow[i, 'm']*(x1 - glow[i, 'start']) + glow[i, 'b'], col='green')
  }
}
############################################

###################### main ######################
ars <- function(f, n, showg = FALSE, upperB=Inf, lowerB=-Inf) {
  ### Sample n points from distribution f using adaptive rejection sampling
  library(numDeriv)
  #Initialize g based on the given f and bounds
  glist <- initG(f, upperB, lowerB)
  acceptedX <- NULL
  while(length(acceptedX) < n) {
    updates <- generatePoints(min(n - length(acceptedX), 10*nrow(glist$Upper)), glist, f)
    #Append acceptedX with newly generated points
    acceptedX <- c(acceptedX, updates$sample)
    glist <- updates$g
  }
  if (showg) {
    plotg(f, glist)
  }
  return(acceptedX)
}
############################################

###################### testing ######################
library(numDeriv)

testUpdateG <- function(f, frand, upperB = Inf, lowerB = -Inf) {
  #### test UpdateG
  logf <- function(x) { log(f(x)) }
  glist <- initG(f, upperB, lowerB)
  gAdd <- c(glist$Upper[1,3], glist$Upper[2,3])
  names(gAdd) <- NULL
  gx <- frand(1e3)
  for(i in 1:length(gx)) {
    glist <- updateG(gx[i], glist, f)
  }
  gx <- c(gx, gAdd)
  gx <- sort(gx)
  gVal <- logf(gx)
  gSlope <- grad(logf, gx)
  if(all.equal(gx, glist$Upper[,3]) == T) {
    print("Intersect values are accurate")
  } else {
    print(all.equal(gx, glist$Upper[,3]))
  }
  if(all.equal(gSlope, glist$Upper[,4]) == T) {
    print("Slope values are accurate")
  } else {
    cat("Slope", (all.equal(gSlope, glist$Upper[,4])), "\n")
  }
  if(all.equal(gVal, glist$Upper[,5]) == T) {
    print("Intercept values are accurate")
  } else {
    print(all.equal(gVal, glist$Upper[,5]))    
  }
}

testSample <- function(testname, gu, trueMean, trueVar, 
                       n=1e5, ks=FALSE, ksdist=NULL) {
  ### test sampleSX
  print(paste("Testing sampleSX() with:", testname))
  y <- sampleSX(gu, n)
  pass <- TRUE
  print(paste("Sample mean:", mean(y)))
  print(paste("Sample var:", var(y)))
  if (abs(mean(y) - trueMean) > .01) {
    print(paste("Failed test", testname, 'with mean', mean(y)))
    print(paste("True mean is", trueMean))
    pass <- FALSE
  }
  if (abs(var(y) - trueVar) > .01) {
    print(paste("Failed test", testname, 'with variance', var(y)))
    print(paste("True variance is", trueVar))
    pass <- FALSE
  }
  if (ks) {
    ksResult <- ks.test(y, ksdist)
    print(ksResult)
    if (ksResult$p.value < .01) {
      print(paste("Failed test", testname, 
                  "with ks test result", ksResult$p.value))
      pass <- FALSE
    }
  }
  if (pass) {
    print(paste("Passed test:", testname))
  }
}

testdist <- function(y, f, testname, trueMean, trueVar){
  #### test ars
  pass <- TRUE
  if (length(y)<1e3){
    stop("sample size should be no less than 1e3.")
  }
  print(paste("sample mean =",mean(y)))
  print(paste("sample variance =",var(y)))
  h <- hist(y, breaks=50)
  xhist <- c(min(h$breaks), h$breaks)
  yhist <- c(0, h$density, 0)
  xfit <- seq(min(y), max(y), length=50)
  yfit <- f(xfit)
  plot(xhist, yhist, type='s', ylim=c(0, max(yhist,yfit)), main='f(x) and histogram of sample points')
  lines(xfit, yfit, col='red')
  if (abs((mean(y) - trueMean)/trueMean) > .05) {
    print(paste("Failed test", testname, 'with mean', mean(y)))
    print(paste("True mean is", trueMean))
    pass <- FALSE
  }
  if (abs((var(y) - trueVar)/trueVar) > .05) {
    print(paste("Failed test", testname, 'with variance', var(y)))
    print(paste("True var is", trueVar))
    pass <- FALSE
  }
  if (pass) {
    print(paste("Passed test:", testname))
  }
}


test <- function() {
  j <- function(x) {dnorm(x)}
  k <- function(x) {dchisq(x, df = 4)}
  l <- function(x) {dgamma(x, shape = 3)}
  m <- function(x) {dbeta(x, shape1 = 3, shape2 = 2)}
  
  jrand <- function(n) {rnorm(n)}
  krand <- function(n) {rchisq(n, df = 4)}
  lrand <- function(n) {rgamma(n, shape = 3)}
  mrand <- function(n) {rbeta(n, shape1 = 3, shape2 = 2)}
  cat("\n***Testing UpdateG: ***\n")
  testUpdateG(j, jrand)
  testUpdateG(k, krand, lowerB = 0)
  testUpdateG(m, mrand, lowerB = 0, upperB = 1)
  testUpdateG(l, lrand, lowerB = 0)

  
  cat("\n*** Testing SampleSX: ***\n")
  # a simple exponential, exp(-x), split into two parts:
  gu <- rbind(c(0, 1.5, .5, -1, -.5), c(1.5, Inf, 2, -1, -2))
  colnames(gu) <- c('start', 'end', 'intersect', 'm', 'b')
  testSample('exponential', gu, 1, 1, n=1e5, ks=TRUE, ksdist=pexp)
  
  # exp(-x) split into three parts:
  gu <- rbind(c(0, 1, .5, -1, -.5), c(1, 2, 1.5, -1, -1.5), 
              c(2, Inf, 2.5, -1, -2.5))
  colnames(gu) <- c('start', 'end', 'intersect', 'm', 'b')
  testSample('exponential in three parts', gu, 1, 1, 
             n=1e5, ks=TRUE, ksdist=pexp)
  
  # two exponentials exp(5x) for x < 0 and exp(-5x) for x > 0
  testSample('two-sided exponential', initG(dnorm)$Upper, 
             0, 2/25, n=1e5)
  
  # same as above,
  # two exponentials exp(5x) for x < 0 and exp(-5x) for x > 0
  # but shifted to the right by one
  gu <- rbind(c(-Inf, 1, -4, 5, -13.4), c(1, Inf, 6, -5, -13.4))
  colnames(gu) <- c('start', 'end', 'intersect', 'm', 'b')
  testSample('shifted two-sided exponenital', gu, 1, 2/25, n=1e5)
  
  # an exponential truncated at 2:
  gu <- rbind(c(0, 1, .5, -1, -.5), c(1, 2, 1.5, -1, -1.5))
  colnames(gu) <- c('start', 'end', 'intersect', 'm', 'b')
  testSample('truncated exponential', gu, .687, .276, n=1e5)
  
  ### 2 exponentials and a flat region
  ### simulates exp(x) for x < 0, unif for 0 < x < 1, exp(-(x-1)) for x > 1
  gu <- rbind(c(-Inf, 0, -1, 1, -1), c(0, 1, .5, 0, 0), 
              c(1, Inf, 2, -1, -1))
  colnames(gu) <- c('start', 'end', 'intersect', 'm', 'b')
  testSample('2 exponentials and a flat region', gu, 
             0.5, 2.195, n=1e6)
  
  cat("\n*** Testing ars() ***\n")
  
  # examples to test whether the sample points are from the given distribution.
  f1 <- function(x) dbeta(x, 2, 2)
  sampleF <- ars(f1, upperB=1, lowerB=0, 10000)
  testdist(sampleF, f1, 'beta', 1/2, 1/20)
  
  f2 <- function(x) dnorm(x, 1, 3)
  sampleF <- ars(f2, 10000)
  testdist(sampleF, f2, 'nomal', 1, 9)
}
