probtwo <- function (RES) {
	if (dim(RES[[2]])[2]==3) {
			bb <- max(RES[[2]][,2])
			aa <- max(RES[[2]][,3])
			cc <- dim(RES[[2]])[1]
			dd <- (cc-(aa-1)*(bb-1))/(aa+bb-1)
			prob <- numeric(cc)
			prob[sapply(seq(0,(aa-1)*(bb-3),aa-1),function (i) i+1:(aa-2))] <- RES[[1]][sapply(seq(0,(aa-1)*(bb-3),aa-1),function (i) i+1:(aa-2))]-RES[[1]][sapply(seq(0,(aa-1)*(bb-3),aa-1),function (i) i+aa+1:(aa-2))]
			prob[sapply(seq(0,(aa-1)*(bb-3),aa-1),function (i) i+aa-1)] <- RES[[1]][sapply(seq(0,(aa-1)*(bb-3),aa-1),function (i) i+aa-1)]-RES[[1]][(aa-1)*(bb-1)+1+1:(bb-2)]
			prob[(aa-1)*(bb-2)+1:(aa-1)] <- RES[[1]][(aa-1)*(bb-2)+1:(aa-1)]-RES[[1]][(aa-1)*(bb-1)+bb+1:(aa-1)]
			prob[(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+1:(bb-2))] <- RES[[1]][(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+1:(bb-2))]-RES[[1]][(aa-1)*(bb-1)+sapply(seq((aa+bb-1),(aa+bb-1)*(dd-1),(aa+bb-1)),function (i) i+1+1:(bb-2))] 
			prob[(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+1:(bb-2))] <- RES[[1]][(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+1:(bb-2))]-RES[[1]][(aa-1)*(bb-1)+sapply(seq((aa+bb-1),(aa+bb-1)*(dd-1),(aa+bb-1)),function (i) i+1+1:(bb-2))]
			prob[(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+bb-1)] <- RES[[1]][(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+bb-1)]-RES[[1]][(aa-1)*(bb-1)+sapply(seq((aa+bb-1),(aa+bb-1)*(dd-1),(aa+bb-1)),function (i) i+1+aa+bb-2)]
			prob[(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+bb-1+1:(aa-1))] <- RES[[1]][(aa-1)*(bb-1)+sapply(seq(0,(aa+bb-1)*(dd-2),(aa+bb-1)),function (i) i+bb-1+1:(aa-1))]-RES[[1]][(aa-1)*(bb-1)+sapply(seq((aa+bb-1),(aa+bb-1)*(dd-1),(aa+bb-1)),function (i) i+bb-1+2:aa)]
			prob[prob==0]<-RES[[1]][prob==0]
			idx <- numeric()
			probnew <- prob[1:((aa-1)*(bb-1))]
			for (i in 1:(dd-1)) {
				if (RES[[2]][(aa-1)*(bb-1)+(i+1)*(aa+bb-1),1]==1) {
					idx <- c(idx,(aa-1)*(bb-1)+i*(aa+bb-1)+1:(aa+bb-1))
					j <- i+1
					tmp <- prob[(aa-1)*(bb-1)+(i-1)*(aa+bb-1)+1:(aa+bb-1)]+prob[(aa-1)*(bb-1)+i*(aa+bb-1)+1:(aa+bb-1)]
					while(RES[[2]][(aa-1)*(bb-1)+(j+1)*(aa+bb-1),1]==1) {
						idx <- c(idx,(aa-1)*(bb-1)+j*(aa+bb-1)+1:(aa+bb-1))
						tmp <- tmp+prob[(aa-1)*(bb-1)+j*(aa+bb-1)+1:(aa+bb-1)]
						j <- j+1
					}
					i <-j
					probnew <- c(probnew,tmp)
				} else {
					probnew <- c(probnew,prob[(aa-1)*(bb-1)+(i-1)*(aa+bb-1)+1:(aa+bb-1)])
				}
			}
			probnew <- c(probnew,prob[(aa-1)*(bb-1)+(dd-1)*(aa+bb-1)+1:(aa+bb-1)])
			if (length(idx)>0) {RES[[2]] <- RES[[2]][-idx,]}
			RES[[2]][((aa-1)*(bb-1)+1):length(RES[[2]][,3]),1] <-as.numeric(sapply(c(1:((length(RES[[2]][,3])-(aa-1)*(bb-1))/(aa+bb-1))),function (i) rep(i,(aa+bb-1))))
			prob <- probnew
			prob[prob<0] <- 0
			prob <- prob/sum(prob)
	} else {
			bb <- max(RES[[2]][,1])
			aa <- max(RES[[2]][,2])
			prob <- numeric(length(RES[[1]]))
			prob[sapply(seq(0,aa*(bb-2),aa),function (i) i+1:(aa-1))] <- RES[[1]][sapply(seq(0,aa*(bb-2),aa),function (i) i+1:(aa-1))]-RES[[1]][sapply(seq(0,aa*(bb-2),aa),function (i) i+aa+2:aa)]
			prob[prob<0] <- 0
			prob <- prob/sum(prob)
	}
	list (prob,RES[[2]])
}