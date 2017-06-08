#function to test significant departure in observed tip branch lengths from the null expectation
addsigtest <- function (nll,obs,dt,n) {
	medX <- ceiling(quantile(obs,0.5)/dt)
	pmed <- sum(sapply(c(1:medX),function (i) sapply(c(0:(ceiling((2*n+2)/4)-1)),function (j) sapply(c(0:(n-ceiling((2*n+2)/4))), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(n-j-z)*factorial(n)/factorial(j)/factorial(z)/factorial(n-j-z)))))
	quadX <- ceiling(quantile(obs,0.25)/dt)
	pquad <- sum(sapply(c(1:quadX),function (i) sapply(c(0:(ceiling((1*n+1)/4)-1)),function (j) sapply(c(0:(n-ceiling((1*n+1)/4))), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(n-j-z)*factorial(n)/factorial(j)/factorial(z)/factorial(n-j-z)))))
	quad3X <- ceiling(quantile(obs,0.75)/dt)
	pquad3 <- sum(sapply(c(1:quad3X),function (i) sapply(c(0:(ceiling((3*n+3)/4)-1)),function (j) sapply(c(0:(n-ceiling((3*n+3)/4))), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(n-j-z)*factorial(n)/factorial(j)/factorial(z)/factorial(n-j-z)))))
	c(pmed,pquad,pquad3)
}