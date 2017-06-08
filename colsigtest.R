#function to test significant departure in the inferred colonization times from the null expectation
colsigtest <- function (nll,obs) {
	medX <- min(which(cumsum(obs)>=0.5))
	pmed <- sum(sapply(c(1:medX),function (i) sapply(c(0:9),function (j) sapply(c(0:9), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(19-j-z)*factorial(19)/factorial(j)/factorial(z)/factorial(19-j-z)))))
	quadX <- min(which(cumsum(obs)>=0.25))
	pquad <- sum(sapply(c(1:quadX),function (i) sapply(c(0:4),function (j) sapply(c(0:14), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(19-j-z)*factorial(19)/factorial(j)/factorial(z)/factorial(19-j-z)))))
	quad3X <- min(which(cumsum(obs)>=0.75))
	pquad3 <- sum(sapply(c(1:quad3X),function (i) sapply(c(0:14),function (j) sapply(c(0:4), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(19-j-z)*factorial(19)/factorial(j)/factorial(z)/factorial(19-j-z)))))
	c(pmed,pquad,pquad3)
}