#function to calculate the likelhood of each subtree of the New Zealand passerine assemblage phylogenies. Note that this function is specifically designed for the topology of the New Zealand passerine phylogenies.
likcal <- function (pars,edge,category,Ngenus,is.AB,is.AB.root,N) {
	if (category!=3) {
		n <- dim(edge)[1]
		prob <- matrix(NA,n,6)
	}
	lg <- 0
	func <- function (time,state,par) {
		with(as.list(c(state,par)), {
			sA <- par[1]
			sB <- par[2]
			sAB <- par[3]
			xA <- par[4]
			xB <- par[5]
			dA <- par[6]
			dB <- par[7]
			deA <- -(sA+dA+xA)*eA+xA+dA*eAB+sA*eA*eA
			deB <- -(sB+dB+xB)*eB+xB+dB*eAB+sB*eB*eB
			deAB <- -(sA+sB+sAB+xA+xB)*eAB+xA*eB+xB*eA+sA*eAB*eA+sB*eAB*eB+sAB*eA*eB
			dnA <- -(sA+dA+xA)*nA+dA*nAB+2*sA*nA*eA
			dnB <- -(sB+dB+xB)*nB+dB*nAB+2*sB*nB*eB
			dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+sA*(eA*nAB+eAB*nA)+sB*(eB*nAB+eAB*nB)+sAB*(eA*nB+eB*nA)
			return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
		})
	}
if (category==1) {
	sampling.f <- N/Ngenus
	if (is.AB) {
		RES <- tryCatch(ode(c(nA=0,nB=0,nAB=1,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[2,4]),func,pars),warning=function(w) w)
		tmp <- try(sum(RES[2,2:4]),silent=T)
		if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
		lg <- lg+log(tmp)
		RES[2,2:4] <- RES[2,2:4]/tmp
		res <- tryCatch(ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[2,4]),func,pars),warning=function(w) w)
		tmp <- try(sum(res[2,2:4]),silent=T)
		if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
		lg <- lg+log(tmp)
		res[2,2:4] <- res[2,2:4]/tmp
		prob <- as.numeric(c(res[2,2]*RES[2,2]*pars[1],res[2,3]*RES[2,3]*pars[2],(res[2,2]*RES[2,4]+RES[2,2]*res[2,4])*pars[1]/2+(res[2,4]*RES[2,3]+RES[2,4]*res[2,3])*pars[2]/2+(res[2,2]*RES[2,3]+RES[2,2]*res[2,3])*pars[3]/2,res[2,5:7]))
		res <- tryCatch(ode(c(nA=prob[1],nB=prob[2],nAB=prob[3],eA=prob[4],eB=prob[5],eAB=prob[6]),c(0,edge[1,which(edge[1,2:3]==edge[2,1])+3]),func,pars),warning=function(w) w)
		if (inherits(res,"simpleWarning")) {stop("negative probability")}
		lik <- res[2,3]
	}
	else {
 	for (i in n:1) {
		if (edge[i,6]<=1 && edge[i,7]<=1) {
			if (edge[i,6]==0 && edge[i,7]==0) {
				res <- tryCatch(ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)*2
				res[2,2:4] <- res[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res[2,2]^2*pars[1],res[2,3]^2*pars[2],res[2,4]*res[2,2]*pars[1]+res[2,4]*res[2,3]*pars[2]+res[2,2]*res[2,3]*pars[3],res[2,5:7]))
			}
			else if (edge[i,6]==1 && edge[i,7]==1) {
				res <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)*2
				res[2,2:4] <- res[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res[2,2]^2*pars[1],res[2,3]^2*pars[2],res[2,4]*res[2,2]*pars[1]+res[2,4]*res[2,3]*pars[2]+res[2,2]*res[2,3]*pars[3],res[2,5:7]))
			}
			else {
				res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*1,nB=(edge[i,6]==0)*sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=(edge[i,7]==1)*1,nB=(edge[i,7]==0)*sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
				}
		}
		else {
			if (edge[i,6]<=1) {
				res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*1,nB=(edge[i,6]==0)*sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
			else if (edge[i,7]<=1) {
				res1 <- tryCatch(ode(c(nA=(edge[i,7]==1)*1,nB=(edge[i,7]==0)*sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
			else {
			res1 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
			tmp <- try(sum(res1[2,2:4]),silent=T)
			if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
			lg <- lg+log(tmp)
			res1[2,2:4] <- res1[2,2:4]/tmp
			res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
			tmp <- try(sum(res2[2,2:4]),silent=T)
			if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
			lg <- lg+log(tmp)
			res2[2,2:4] <- res2[2,2:4]/tmp
			prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
		}
	}
	lik <- prob[1,2]
	}
}
if (category==2) {
	for (i in n:2) {
		if (edge[i,6]==1 && edge[i,7]==1) {
				res <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)*2
				res[2,2:4] <- res[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res[2,2]^2*pars[1],res[2,3]^2*pars[2],res[2,4]*res[2,2]*pars[1]+res[2,4]*res[2,3]*pars[2]+res[2,2]*res[2,3]*pars[3],res[2,5:7]))
			}
		else {
			if (edge[i,6]==1) {
				res1 <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
			else if (edge[i,7]==1) {
				res1 <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
			else {
			res1 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
			tmp <- try(sum(res1[2,2:4]),silent=T)
			if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
			lg <- lg+log(tmp)
			res1[2,2:4] <- res1[2,2:4]/tmp
			res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
			tmp <- try(sum(res2[2,2:4]),silent=T)
			if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
			lg <- lg+log(tmp)
			res2[2,2:4] <- res2[2,2:4]/tmp
			prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
			}
		}
	}
	res <- tryCatch(ode(c(nA=prob[2,1],nB=prob[2,2],nAB=prob[2,3],eA=prob[2,4],eB=prob[2,5],eAB=prob[2,6]),c(0,edge[i,c(4,5)[which(edge[1,6:7]==2)]]),func,pars),silent=T)
	if (inherits(res,"simpleWarning")) {stop("negative probability")}
	lik <- res[2,3]
	if (is.AB.root) {lik <- (res[2,2]*pars[3]/2+res[2,3]*pars[2]+res[2,4]*pars[2]/2)/(pars[3]/2+pars[2]+pars[2]/2)}
	}	
if (category==3) {
	res <- tryCatch(ode(c(nA=(!is.AB)*1,nB=0,nAB=(is.AB)*1,eA=0,eB=0,eAB=0),c(0,edge[(edge[2]>edge[3])+4]),func,pars),silent=T)
	if (inherits(res,"simpleWarning")) {stop("negative probability")}
	lik <- res[2,3]
	if (is.AB.root) {lik <- (res[2,2]*pars[3]/2+res[2,3]*pars[2]+res[2,4]*pars[2]/2)/(pars[3]/2+pars[2]+pars[2]/2)}
}
log(lik)+lg
}