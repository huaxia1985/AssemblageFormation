#function to infer colonisation times from New Zealand passerine assemblage phylogenies. Note that this function is specifically designed for the topology of the New Zealand passerines phylogenies.
colcal <- function(pars,edge,category,Ngenus,is.AB,dt,N) {
	func0 <- function (time,state1,par) {
		with(as.list(c(state1,par)), {
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
			dnB <- -(sB+dB+xB)*nB+0+2*sB*nB*eB
			dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+sA*(eA*nAB+eAB*nA)+sB*(eB*nAB+eAB*nB)+sAB*(eA*nB+eB*nA)
			return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
		})
	}
	func <- function (time,state1,par) {
		with(as.list(c(state1,par)), {
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
	dtcal <- function (state,pars,dt,t,t1) {
		if (t1>=t) {
			pp <- numeric(2)
			state1 <- matrix(NA,2,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func0,pars)
			pp[2] <- log(sum(tmp[2,2:4]))
			state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			t1 <- t1-t
		} else {
			Tt <- trunc((t-t1)/dt)
			pp <- numeric(Tt+3)
			state1 <- matrix(NA,Tt+3,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- pp[1]+log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func0,pars)
			pp[Tt+3] <- pp[Tt+3]+log(sum(tmp[2,2:4]))
			state1[Tt+3,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			if (t1>0) {
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t1),func0,pars)
				pp[2:(Tt+2)] <- pp[2:(Tt+2)]+log(sum(tmp[2,2:4]))
				state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(t1,t),func,pars)
				pp[2] <- pp[2]+log(sum(tmp[2,2:4]))
				state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			}
			if (Tt!=0) {
				for (i in 1:Tt) {
					tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,i*dt),func0,pars)
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
					state2 <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
					tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
					state1[i+2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				}
			}
				t1 <- dt-(t-Tt*dt-t1)
			if (is.na(state1[2,1])) {
				pp <- pp[-2]
				state1 <- state1[-2,]
			}
		}
		list(pp,state1,t1)
	}
		dtcal2 <- function (state,pars,dt,t,t1) {
		if (t1>=t) {
			pp <- numeric(2)
			state1 <- matrix(NA,1,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			t1 <- t1-t
		} else {
			Tt <- trunc((t-t1)/dt)
			pp <- numeric(Tt+2)
			state1 <- matrix(NA,Tt+2,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- pp[1]+log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			if (t1>0) {
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t1),func0,pars)
				pp[2:(Tt+2)] <- pp[2:(Tt+2)]+log(sum(tmp[2,2:4]))
				state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(t1,t),func,pars)
				pp[2] <- pp[2]+log(sum(tmp[2,2:4]))
				state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			}
			if (Tt!=0) {
				for (i in 1:Tt) {
					tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,i*dt),func0,pars)
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
					state2 <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
					tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
					state1[i+2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				}
			}
				t1 <- dt-(t-Tt*dt-t1)
			if (is.na(state1[2,1])) {
				pp <- pp[-2]
				state1 <- state1[-2,]
			}
		}
		list(pp,state1,t1)
	}
	odecal <- function (state,t,pars) {
		tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
		pp <- sum(tmp[2,2:4])
		state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
		list(pp,state)
	}
	if (category==1) {
		sampling.f <- N/Ngenus
		if (is.AB) {
			tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[2,4]),func,pars)
			tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
			res <- dtcal(c(0,0,1,0,1-sampling.f,0),pars,dt,edge[2,4],0)
			state <- res[[2]]
			state[,3] <- (state[,3]*tmp[1]+state[,1]*tmp[3])*pars[1]/2+(state[,2]*tmp[3]+state[,3]*tmp[2])*pars[2]/2+(state[,1]*tmp[2]+state[,2]*tmp[1])*pars[3]/2
			state[,1] <- state[,1]*tmp[1]*pars[1]
			state[,2] <- state[,2]*tmp[2]*pars[2]
			res[[1]] <- res[[1]]+log(rowSums(state[,1:3]))
			state[,1:3] <- state[,1:3]/rowSums(state[,1:3])
			ll <- dim(state)[1]
			res1 <- dtcal(state[ll,],pars,dt,edge[1,which(edge[1,2:3]==edge[2,1])+3],res[[3]])
			res2 <- apply(matrix(state[-ll,],ll-1,6),1,odecal,edge[1,which(edge[1,2:3]==edge[2,1])+3],pars)
			prob <- c(sapply(c(1:length(res2)),function (i) exp(res[[1]][i])*res2[[i]][[1]]*res2[[i]][[2]][2]/sum(res2[[i]][[2]][1:3])),exp(res[[1]][ll])*exp(res1[[1]])*res1[[2]][,2]/rowSums(res1[[2]][,1:3]))
			prob <- prob[-length(prob)]-prob[-1]
			if (edge[2,4]<dt) {
				prob[1] <- prob[1]+prob[2]
				prob <- prob[-2]
			} else {
				TT <- trunc(edge[2,4]/dt)
				prob[TT+1] <- prob[TT+1]+prob[TT+2]
				prob <- prob[-c(TT+2)]
			}
			prob <- prob/sum(prob)
			RES <- vector("list",4)
			RES[[4]] <-NA
			ll<- NA
		} else {
			n <- dim(edge)[1]
			state <- vector("list",n)
			res <- vector("list",n)
			for (i in n:1) {
				if (edge[i,6]==0 && edge[i,7]==0) {
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars)
					state[[i]] <- as.numeric(c(tmp[2,2]^2*pars[1],tmp[2,3]^2*pars[2],tmp[2,4]*tmp[2,2]*pars[1]+tmp[2,4]*tmp[2,3]*pars[2]+tmp[2,2]*tmp[2,3]*pars[3],tmp[2,5:7]))
					state[[i]][1:3] <- state[[i]][1:3]/sum(state[[i]][1:3])
					state[[i]] <- matrix(state[[i]],1,6)
				} else if (edge[i,6]<=1 && edge[i,7]<=1) {
					res[[i]] <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,edge[i,4],0)
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					Tt <- edge[i,4]
				} else if (edge[i,6]==1 && edge[i,7]>0) {
					res[[i]] <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,edge[i,4],0)
					tmp <- state[[edge[i,7]]]
					tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,5]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					Tt <- edge[i,4]
				} else if (edge[i,6]==0 && edge[i,7]>0) {
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					if (dim(state[[edge[i,7]]])[1]==1) {
						res1 <- state[[edge[i,7]]]
						res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,edge[i,5]),func,pars)
						state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					} else {
						ll <- dim(state[[edge[i,7]]])[1]
						res1 <- dtcal(state[[edge[i,7]]][ll,],pars,dt,edge[i,5],res[[edge[i,7]]][[3]])
						res2 <- apply(matrix(state[[edge[i,7]]][-ll,],(ll-1),6),1,odecal,edge[i,5],pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[edge[i,7]]][[1]][j]),(res[[edge[i,7]]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						Tt <- c(Tt,edge[i,5])
						state[[i]] <- res[[i]][[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					}
				} else if (edge[i,6]>0 && edge[i,7]==1) {
					res[[i]] <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,edge[i,5],0)
					tmp <- state[[edge[i,6]]]
					tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,4]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					Tt <- edge[i,5]
				} else if (edge[i,6]>0 && edge[i,7]==0) {
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,5]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					if (dim(state[[edge[i,6]]])[1]==1) {
						res1 <- state[[edge[i,6]]]
						res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,edge[i,4]),func,pars)
						state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					} else {
						ll <- dim(state[[edge[i,6]]])[1]
						res1 <- dtcal(state[[edge[i,6]]][ll,],pars,dt,edge[i,4],res[[edge[i,6]]][[3]])
						res2 <- apply(matrix(state[[edge[i,6]]][-ll,],(ll-1),6),1,odecal,edge[i,4],pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[edge[i,6]]][[1]][j]),(res[[edge[i,6]]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						Tt <- c(Tt,edge[i,4])
						state[[i]] <- res[[i]][[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					}
				} else {
					if (dim(state[[edge[i,6]]])[1]==1 && dim(state[[edge[i,7]]])[1]==1) {
						res1 <- state[[edge[i,6]]]
						res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,edge[i,4]),func,pars)
						tmp <- state[[edge[i,7]]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,5]),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					} else if (dim(state[[edge[i,6]]])[1]>1 && dim(state[[edge[i,7]]])[1]==1) {
						tmp <- state[[edge[i,7]]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,5]),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						ll <- dim(state[[edge[i,6]]])[1]
						res1 <- dtcal(state[[edge[i,6]]][ll,],pars,dt,edge[i,4],res[[edge[i,6]]][[3]])
						res2 <- apply(matrix(state[[edge[i,6]]][-ll,],(ll-1),6),1,odecal,edge[i,4],pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[edge[i,6]]][[1]][j]),(res[[edge[i,6]]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						Tt <- c(Tt,edge[i,4])
						state[[i]] <- res[[i]][[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					} else if (dim(state[[edge[i,6]]])[1]==1 && dim(state[[edge[i,7]]])[1]>1) {
						tmp <- state[[edge[i,6]]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,4]),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						ll <- dim(state[[edge[i,7]]])[1]
						res1 <- dtcal(state[[edge[i,7]]][ll,],pars,dt,edge[i,5],res[[edge[i,7]]][[3]])
						res2 <- apply(matrix(state[[edge[i,7]]][-ll,],(ll-1),6),1,odecal,edge[i,5],pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[edge[i,7]]][[1]][j]),(res[[edge[i,7]]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						Tt <- c(Tt,edge[i,5])
						state[[i]] <- res[[i]][[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					}
				}
			}
		Tt <- Tt[-length(Tt)]
		prob <- exp(res[[1]][[1]])*state[[1]][,2]
		prob <- prob[-length(prob)]-prob[-1]
		t1 <- 0
		prob1 <-numeric()
		for (i in Tt) {
			if (i<t1) {
				t1 <- t1-i
				prob1[length(prob1)] <- prob1[length(prob1)]+prob[1]
				prob <- prob[-1]
			} else if ((i-t1)<dt) {
				t1 <- dt-i+t1
				prob1 <- c(prob1,prob[1]+prob[2])
				prob <- prob[-c(1:2)]
			} else {
				TT <- trunc((i-t1)/dt)
				t1 <- dt-(i-t1-TT*dt)
				prob1 <- c(prob1,prob[1:(TT)],prob[TT+1]+prob[TT+2])
				prob <- prob[-c(1:(TT+2))]
			}
		}
		prob1 <- c(prob1,prob)
		prob1[prob1<0] <- 0
		prob <- prob1/sum(prob1)
		RES <- vector("list",4)
		RES[[4]] <-NA
		ll<- NA
		}
	}
	if (category==2) {
			n <- dim(edge)[1]
			res <- vector("list",n)
			for (i in n:2) {
				if (edge[i,6]==1 && edge[i,7]==1) {
					tmp <- dtcal(c(1,0,0,0,0,0),pars,dt,edge[i,4],0)
					ll <- length(tmp[[1]])
					idx <- sapply(c(1:ll),function(j) rep(j,ll))
					idx <- cbind(as.numeric(idx),as.numeric(t(idx)))
					state <- matrix(NA,dim(idx)[1],6)
					state[,1] <- tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
					state[,2] <- tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
					state[,3] <- (tmp[[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+tmp[[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
					state[,4] <- tmp[[2]][1,4]
					state[,5] <- tmp[[2]][1,6]
					state[,6] <- tmp[[2]][1,6]
					res[[i]] <- vector("list",7)
					res[[i]][[1]] <- tmp[[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(state[,1:3]))
					state[,1:3] <- state[,1:3]/rowSums(state[,1:3])
					res[[i]][[2]] <- state
					res[[i]][[3]] <- tmp[[3]]
					res[[i]][[4]] <- idx
					res[[i]][[5]] <- cbind(seq(ll,ll^2,ll),(ll*(ll-1)+1):(ll*ll))
					res[[i]][[6]] <- cbind(unlist(sapply(c((ll-1):1),function (i) seq(ll-i,ll*(ll-i),ll))),unlist(sapply(c(1:(ll-1)),function (i) (ll*(ll-1-i)+1):(ll*(ll-i)-i))))
					res[[i]][[7]] <- c(seq(ll,ll*(ll-1),ll),(ll*(ll-1)+1):(ll*ll))
				} else {
					a <- ifelse (edge[i,6]>1,edge[i,6],edge[i,7])
					tmp <- dtcal(c(1,0,0,0,0,0),pars,dt,ifelse(edge[i,6]>1,edge[i,5],edge[i,4]),0)
					res1 <- vector("list",max(res[[a]][[5]]))
					res1[res[[a]][[5]][,1]] <- apply(matrix(res[[a]][[2]][res[[a]][[5]][,1],],length(res[[a]][[5]][,1]),6),1,dtcal,pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
					res1[res[[a]][[5]][,2]] <- res1[res[[a]][[5]][,1]]
					res1[res[[a]][[6]][,1]] <- apply(matrix(res[[a]][[2]][res[[a]][[6]][,1],],length(res[[a]][[6]][,1]),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
					res1[res[[a]][[6]][,2]] <- res1[res[[a]][[6]][,1]]
					res2 <- res1[-res[[a]][[7]]]
					res1 <- res1[res[[a]][[7]]]
					res[[i]] <- vector("list",5)
					res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[7]]][j]),as.numeric(res[[a]][[1]][res[[a]][[7]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[7]])),function (j) res1[[j]][[1]][z]))))
					res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
					res[[i]][[3]] <- res1[[3]][[3]]
					res[[i]][[4]] <- rbind(cbind(0,matrix(res[[a]][[4]][-res[[a]][[7]],],dim(res[[a]][[4]])[1]-length(res[[a]][[7]]),dim(res[[a]][[4]])[2])),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[7]])))),apply(res[[a]][[4]][res[[a]][[7]],],2,rep,length(res1[[1]][[1]]))))
					res[[i]][[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[7]])+1:length(res[[a]][[7]])
					idx <- rbind(cbind(rep(c(1:(length(res[[i]][[1]])-length(res[[i]][[5]]))),length(tmp[[1]])),as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]]))))),cbind(rep((length(res[[i]][[1]])-length(res[[i]][[5]])+1):length(res[[i]][[1]]),length(tmp[[1]])),as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[5]]))))))
					state <- matrix(NA,dim(idx)[1],6)
					state[,1] <- res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
					state[,2] <- res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
					state[,3] <- (res[[i]][[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+res[[i]][[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
					state[,4] <- tmp[[2]][1,4]
					state[,5] <- tmp[[2]][1,6]
					state[,6] <- tmp[[2]][1,6]
					RES <- vector("list",7)
					RES[[1]] <- res[[i]][[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(state[,1:3]))
					state[,1:3] <- state[,1:3]/rowSums(state[,1:3])
					RES[[2]] <- state
					RES[[3]] <- tmp[[3]]
					RES[[4]] <- rbind(cbind(as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]])))),apply(res[[i]][[4]][-res[[i]][[5]],],2,rep,length(tmp[[1]]))),cbind(as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[5]])))),apply(res[[i]][[4]][res[[i]][[5]],],2,rep,length(tmp[[1]]))))
					if (i==(n-1)) {
						aa <- (ll-1)^2+(2*ll-1)*(length(res1[[1]][[1]])-1)
						RES[[5]] <- rbind(cbind(as.numeric((length(tmp[[1]])-1)*aa+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1)))),as.numeric((length(tmp[[1]])-1)*aa+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1)))),cbind(as.numeric((length(tmp[[1]])-1)*aa+as.numeric(sapply(c(1:(length(res1[[1]][[1]])-1)),function (i) (ll-1)^2+(i-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),as.numeric((length(tmp[[1]])-1)*aa+sapply(c(1:(length(res1[[1]][[1]])-1)),function (i) (ll-1)^2+(i-1)*(2*ll-1)+c(ll:(2*ll-1))))),cbind(as.numeric(sapply(c(1:length(tmp[[1]])), function (z) (z-1)*(2*ll-1)+length(tmp[[1]])*aa+c(1:(ll-1),2*ll-1))),as.numeric(sapply(c(1:length(tmp[[1]])), function (z) (z-1)*(2*ll-1)+length(tmp[[1]])*aa+c(ll:(2*ll-1))))))
						RES[[6]] <- rbind(cbind(as.numeric(sapply(c(1:(length(tmp[[1]])-1)), function (z) (z-1)*aa+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1))))),as.numeric(sapply(c(1:(length(tmp[[1]])-1)), function (z) (z-1)*aa+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1))))),cbind(as.numeric(sapply(c(1:(length(tmp[[1]])-1)), function (z) (z-1)*aa+as.numeric(sapply(c(1:(length(res1[[1]][[1]])-1)),function (i) (ll-1)^2+(i-1)*(2*ll-1)+c(1:(ll-1),2*ll-1))))),as.numeric(sapply(c(1:(length(tmp[[1]])-1)), function (z) (z-1)*aa+sapply(c(1:(length(res1[[1]][[1]])-1)),function (i) (ll-1)^2+(i-1)*(2*ll-1)+c(ll:(2*ll-1)))))))
					}
					if (i==(n-2)) {
						aa <- (ll-1)^2+(2*ll-1)*(max(RES[[4]][,4])-1)
						bb <- aa*max(RES[[4]][,3])
						cc <- bb+(2*ll-1)*max(RES[[4]][,3])
						dd <- aa*(max(RES[[4]][,3])-1)+((ll-1)^2+(2*ll-1)*(max(RES[[4]][,4])-1)+(2*ll-1)*max(RES[[4]][,3]))*(max(RES[[4]][,2])-1)
						ee <- dd*(max(RES[[4]][,1])-1)
						RES[[5]] <-rbind(cbind(ee+as.numeric(sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1))))),ee+as.numeric(sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1))))),cbind(ee+as.numeric(sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),ee+as.numeric(sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1)))))),cbind(ee+aa*(max(RES[[4]][,3])-1)+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1))))),ee+aa*(max(RES[[4]][,3])-1)+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1))))),cbind(ee+aa*(max(RES[[4]][,3])-1)+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),ee+aa*(max(RES[[4]][,3])-1)+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1)))))),cbind(ee+bb+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),ee+bb+as.numeric(sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1)))))),cbind(ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1))))),ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1))))),cbind(ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1)))))),cbind(ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+(2*ll-1)*(max(RES[[4]][,4])-1)+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1)))),ee+dd+as.numeric(sapply(c(0:(max(RES[[4]][,1])-1)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+(2*ll-1)*(max(RES[[4]][,4])-1)+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1)))))))
						RES[[6]] <- rbind(cbind(as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1)))))),as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1)))))),cbind(as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1))))),as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+sapply(c(1:(max(RES[[4]][,3])-1)), function (z) (z-1)*aa+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1))))))),cbind(as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+aa*(max(RES[[4]][,3])-1)+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(1:(ll-1)),function (i) c(0:(ll-i-1))+i+(ll-1)*(i-1)))))),as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+aa*(max(RES[[4]][,3])-1)+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+unlist(sapply(c(0:(ll-2)),function (i) seq((ll-1)*i,(ll-1)*(ll-2),(ll-1))+i+1)))))),cbind(as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+aa*(max(RES[[4]][,3])-1)+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1))))),as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+aa*(max(RES[[4]][,3])-1)+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+(ll-1)^2+sapply(c(1:(max(RES[[4]][,4])-1)), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1))))))),cbind(as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+bb+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(1:(ll-1),2*ll-1))))),as.numeric(sapply(c(1:(max(RES[[4]][,1])-1)), function (z) (z-1)*dd+bb+sapply(c(0:(max(RES[[4]][,2])-2)), function (z) z*(aa+(2*ll-1)*max(RES[[4]][,3]))+sapply(c(1:max(RES[[4]][,3])), function (z) (z-1)*(2*ll-1)+c(ll:(2*ll-1))))))))
						}
						RES[[7]] <- (length(res[[i]][[1]])-length(res[[i]][[5]]))*(length(tmp[[1]])-1)+1:(length(res[[i]][[1]])-length(res[[i]][[5]])+length(tmp[[1]])*length(res[[i]][[5]]))
						res[[i]] <- RES
						}
			}
					res1 <- vector("list",max(res[[2]][[7]]))
					res1[res[[2]][[5]][,1]] <- apply(matrix(res[[2]][[2]][res[[2]][[5]][,1],],length(res[[2]][[5]][,1]),6),1,dtcal2,pars,dt,edge[1,which(edge[1,6:7]==2)+3],res[[2]][[3]])
					res1[res[[2]][[5]][,2]] <- res1[res[[2]][[5]][,1]]
					res1[res[[2]][[6]][,1]] <- apply(matrix(res[[2]][[2]][res[[2]][[6]][,1],],length(res[[2]][[6]][,1]),6),1,odecal,edge[1,which(edge[1,6:7]==2)+3],pars)
					res1[res[[2]][[6]][,2]] <- res1[res[[2]][[6]][,1]]
					res2 <- res1[-res[[2]][[7]]]
					res1 <- res1[res[[2]][[7]]]
					RES <- vector("list",5)
					RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[2]][[1]][-res[[2]][[7]]][j]),as.numeric(sapply(c(1:length(res1)),function (j) res[[2]][[1]][res[[2]][[7]][j]]+res1[[j]][[1]])))
					RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1)),function (i) as.matrix(res1[[1]][[2]]))))
					RES[[1]] <- exp(RES[[1]])*RES[[2]][,2]
					RES[[4]] <- rbind(cbind(0,matrix(res[[2]][[4]][-res[[2]][[7]],],dim(res[[2]][[4]])[1]-length(res[[2]][[7]]),dim(res[[2]][[4]])[2])),cbind(rep(c(1:length(res1[[1]][[1]])),length(res[[2]][[7]])),apply(res[[2]][[4]][res[[2]][[7]],],2,function (i) rep(i,rep(length(res1[[1]][[1]]),length(res[[2]][[7]]))))))
				if (n==3) {prob <- probthree(RES)}
				if (n==4) {prob <- probfour(RES)}
				a <- dim(RES[[4]])
				ll <- sapply(c(1:a[2]),function (i) max(RES[[4]][,i]))
			}
	if (category==3) {
		RES <- dtcal2(c((!is.AB)*1,0,(is.AB)*1,0,0,0),pars,dt,max(edge[4:5]),0)
		prob <- exp(RES[[1]])*RES[[2]][,2]
		prob <- prob[-length(prob)]-prob[-1]
		prob[prob<0] <- 0
		prob <- prob/sum(prob)
		RES[[4]] <- NA
		ll <- NA
		}
	if (category==4) {
		sampling.f <- N/Ngenus
		n <- dim(edge)[1]
			state <- vector("list",n)
			res <- vector("list",n)
			Tt <- vector("list",2)
			for (i in n:1) {
				if (edge[i,6]==0 && edge[i,7]==0) {
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars)
					state[[i]] <- as.numeric(c(tmp[2,2]^2*pars[1],tmp[2,3]^2*pars[2],tmp[2,4]*tmp[2,2]*pars[1]+tmp[2,4]*tmp[2,3]*pars[2]+tmp[2,2]*tmp[2,3]*pars[3],tmp[2,5:7]))
					state[[i]][1:3] <- state[[i]][1:3]/sum(state[[i]][1:3])
					state[[i]] <- matrix(state[[i]],1,6)
				} else if (edge[i,6]==1 && edge[i,7]==1) {
					tmp <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,edge[i,4],0)
					ll <- length(tmp[[1]])
					idx <- sapply(c(1:ll),function(j) rep(j,ll))
					idx <- cbind(as.numeric(idx),as.numeric(t(idx)))
					state[[i]] <- matrix(NA,dim(idx)[1],6)
					state[[i]][,3] <- (tmp[[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+tmp[[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
					state[[i]][,1] <- tmp[[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
					state[[i]][,2] <- tmp[[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
					state[[i]][,4] <- tmp[[2]][1,4]
					state[[i]][,5] <- tmp[[2]][1,6]
					state[[i]][,6] <- tmp[[2]][1,6]
					res[[i]] <- vector("list",4)
					res[[i]][[1]] <- tmp[[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(state[[i]][,1:3]))
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					res[[i]][[2]] <- idx
					res[[i]][[3]] <- c(seq(ll,ll*(ll-1),ll),(ll*(ll-1)+1):(ll*ll))
					res[[i]][[4]] <- tmp[[3]]
					Tt <- list(c(i,edge[i,4]),c(i,edge[i,4]))
				} else if ((edge[i,6]==0 && edge[i,7]==1)||(edge[i,6]==1 && edge[i,7]==0)) {
					tmp1 <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,edge[i,4],0)
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,edge[i,4]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- tmp1[[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- tmp1[[1]]*rowSums(state[[i]][,1:3])
					res[[i]][[2]] <- c(1:length(tmp1[[1]]))
					res[[i]][[3]] <- length(tmp1[[1]])
					res[[i]][[4]] <- tmp1[[3]]
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					Tt[[ifelse(is.null(Tt[[1]]),1,2)]] <- c(i,edge[i,4])
				} else if ((edge[i,6]==1 && edge[i,7]>1)||(edge[i,6]>1 && edge[i,7]==1)) {
					a <- ifelse (edge[i,6]>1,edge[i,6],edge[i,7])
					tmp1 <- dtcal(c(1,0,0,0,1-sampling.f,0),pars,dt,ifelse(edge[i,6]==1,edge[i,4],edge[i,5]),0)
					if (dim(state[[a]])[1]==1) {
						tmp <- state[[a]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,ifelse(edge[i,6]==1,edge[i,5],edge[i,4])),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						state[[i]] <- tmp1[[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- tmp1[[1]]*rowSums(state[[i]][,1:3])
						res[[i]][[2]] <- c(1:length(tmp1[[1]]))
						res[[i]][[3]] <- length(tmp1[[1]])
						res[[i]][[4]] <- tmp1[[3]]
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					} else {
						res1 <- dtcal(state[[a]][res[[a]][[3]],],pars,dt,ifelse(edge[i,6]==1,edge[i,5],edge[i,4]),res[[a]][[4]])
						res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,ifelse(edge[i,6]==1,edge[i,5],edge[i,4]),pars)
						tmp <- vector("list",2)
						tmp[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),res[[a]][[1]][res[[a]][[3]]]+res1[[1]])
						tmp[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						idx <- cbind(as.numeric(sapply(c(1:length(tmp1[[1]])),function(j) rep(j,length(tmp[[1]])))),as.numeric(t(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(tmp1[[1]]))))))
						state[[i]] <- matrix(NA,dim(idx)[1],6)
						state[[i]][,3] <- (tmp1[[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+tmp1[[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
						state[[i]][,1] <- tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
						state[[i]][,2] <- tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
						state[[i]][,4] <- tmp1[[2]][1,4]
						state[[i]][,5] <- tmp1[[2]][1,6]
						state[[i]][,6] <- tmp1[[2]][1,6]
 						res[[i]] <- vector("list",4)
						res[[i]][[1]] <- tmp1[[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(state[[i]][,1:3]))
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
						res[[i]][[2]] <- idx
						res[[i]][[3]] <- c(seq(length(tmp[[1]]),length(tmp[[1]])*(length(tmp1[[1]])-1),length(tmp[[1]])),(length(tmp[[1]])*(length(tmp1[[1]])-1)+1):(length(tmp[[1]])*length(tmp1[[1]])))
						res[[i]][[4]] <- res1[[3]]
						Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]] <- cbind(Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]],c(i,ifelse(edge[i,6]==1,edge[i,5],edge[i,4])))
					}
					Tt[[ifelse(is.element(a,Tt[[1]]),2,1)]] <- c(i,ifelse(edge[i,6]==1,edge[i,4],edge[i,5]))
				} else if ((edge[i,6]==0 && edge[i,7]>1)||(edge[i,7]==0 && edge[i,6]>1)) {
					a <- ifelse(edge[i,6]==0,edge[i,7],edge[i,6])
					tmp <- ode(c(nA=0,nB=sampling.f,nAB=0,eA=0,eB=1-sampling.f,eAB=0),c(0,ifelse(edge[i,6]==0,edge[i,4],edge[i,5])),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					if (dim(state[[a]])[1]==1) {
						res1 <- state[[a]]
						res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,ifelse(edge[i,6]==0,edge[i,5],edge[i,4])),func,pars)
						state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					} else if (length(res[[a]][[2]])==length(res[[a]][[1]])) {
						res1 <- dtcal(state[[a]][res[[a]][[3]],],pars,dt,ifelse(edge[i,6]==0,edge[i,5],edge[i,4]),res[[a]][[4]])
						res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,ifelse(edge[i,6]==0,edge[i,5],edge[i,4]),pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),res[[a]][[1]][res[[a]][[3]]]+res1[[1]])
						res[[i]][[2]] <- c(res[[a]][[2]],max(res[[a]][[2]])+1:(length(res1[[1]])-1))
						res[[i]][[3]] <- max(res[[i]][[2]])
						res[[i]][[4]] <- res1[[3]]
 						state[[i]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
						Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]] <- cbind(Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]],c(i,ifelse(edge[i,6]==0,edge[i,5],edge[i,4])))
					} else {
						res1 <- apply(state[[a]][res[[a]][[3]],],1,dtcal,pars,dt,ifelse(edge[i,6]==0,edge[i,5],edge[i,4]),res[[a]][[4]])
						res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,ifelse(edge[i,6]==0,edge[i,5],edge[i,4]),pars)
						res[[i]] <- vector("list",4)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),as.numeric(res[[a]][[1]][res[[a]][[3]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[3]])),function (j) res1[[j]][[1]][z]))))
						if (dim(res[[a]][[2]])[2]==2) {
								res[[i]][[2]] <- rbind(cbind(0,res[[a]][[2]][-res[[a]][[3]],]),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[3]])))),apply(res[[a]][[2]][res[[a]][[3]],],2,rep,length(res1[[1]][[1]]))))
						} else {
								res[[i]][[2]] <- rbind(res[[a]][[2]][-res[[a]][[3]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[3]])))),apply(res[[a]][[2]][,c(2:3)][res[[a]][[3]],],2,rep,length(res1[[1]][[1]]))))
						}
						res[[i]][[3]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[3]])+1:length(res[[a]][[3]])
						state[[i]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
						res[[i]][[4]] <- res1[[1]][[3]]
 					 	state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						res[[i]][[1]] <- res[[i]][[1]]+log(rowSums(state[[i]][,1:3]))
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
					}
				} else {
					if (dim(state[[edge[i,6]]])[1]==1 && dim(state[[edge[i,7]]])[1]==1) {
						res1 <- state[[edge[i,6]]]
						res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,edge[i,4]),func,pars)
						tmp <- state[[edge[i,7]]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,5]),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					} else if ((dim(state[[edge[i,6]]])[1]>1 && dim(state[[edge[i,7]]])[1]==1)||(dim(state[[edge[i,6]]])[1]==1 && dim(state[[edge[i,7]]])[1]>1)) {
						a <- ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,7],edge[i,6])
						tmp <- state[[ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,6],edge[i,7])
]]
						tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,4],edge[i,5])),func,pars)
						tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
						if (length(res[[a]][[2]])==length(res[[a]][[1]])) {
							res1 <- dtcal(state[[a]][res[[a]][[3]],],pars,dt,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,5],edge[i,4]),res[[a]][[4]])
							res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,5],edge[i,4]),pars)
							res[[i]] <- vector("list",3)
							res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),res[[a]][[1]][res[[a]][[3]]]+res1[[1]])
							res[[i]][[2]] <- c(res[[a]][[2]],max(res[[a]][[2]])+1:(length(res1[[1]])-1))
							res[[i]][[3]] <- max(res[[i]][[2]])
							res[[i]][[4]] <- res1[[3]]
 							state[[i]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
							state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
							state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
							state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
							res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
							state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
							Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]] <- cbind(Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]],c(i,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,5],edge[i,4])))
						} else {
							res1 <- apply(matrix(state[[a]][res[[a]][[3]],],length(res[[a]][[3]]),6),1,dtcal,pars,dt,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,5],edge[i,4]),res[[a]][[4]])
							res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,ifelse(dim(state[[edge[i,6]]])[1]==1,edge[i,5],edge[i,4]),pars)
							res[[i]] <- vector("list",4)
							res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),as.numeric(res[[a]][[1]][res[[a]][[3]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[3]])),function (j) res1[[j]][[1]][z]))))
							if (dim(res[[a]][[2]])[2]==2) {
								res[[i]][[2]] <- rbind(cbind(0,res[[a]][[2]][-res[[a]][[3]],]),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[3]])))),apply(res[[a]][[2]][res[[a]][[3]],],2,rep,length(res1[[1]][[1]]))))
							} else {
								res[[i]][[2]] <- rbind(res[[a]][[2]][-res[[a]][[3]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[3]])))),apply(res[[a]][[2]][,c(2:3)][res[[a]][[3]],],2,rep,length(res1[[1]][[1]]))))
							}
							res[[i]][[3]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[3]])+1:length(res[[a]][[3]])
							state[[i]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
							res[[i]][[3]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[3]])+1:length(res[[a]][[3]])
							res[[i]][[4]] <- res1[[1]][[3]]
 						 	state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
							state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
							state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
							res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
							state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
						}
					} else {
						a <- edge[i,6]
						b <- edge[i,7]
						res1 <- dtcal(state[[a]][res[[a]][[3]],],pars,dt,edge[i,4],res[[a]][[4]])
						res2 <- apply(matrix(state[[a]][-res[[a]][[3]],],(max(res[[a]][[3]])-length(res[[a]][[3]])),6),1,odecal,edge[i,4],pars)
						tmp1 <- vector("list",2)
						tmp1[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[3]]][j]),res[[a]][[1]][res[[a]][[3]]]+res1[[1]])
 						tmp1[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
 						res1 <- dtcal(state[[b]][res[[b]][[3]],],pars,dt,edge[i,5],res[[b]][[4]])
						res2 <- apply(matrix(state[[b]][-res[[b]][[3]],],(max(res[[b]][[3]])-length(res[[b]][[3]])),6),1,odecal,edge[i,5],pars)
						tmp <- vector("list",2)
						tmp[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[b]][[1]][-res[[b]][[3]]][j]),res[[b]][[1]][res[[b]][[3]]]+res1[[1]])
						tmp[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						idx <- cbind(as.numeric(sapply(c(1:length(tmp1[[1]])),function(j) rep(j,length(tmp[[1]])))),as.numeric(t(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(tmp1[[1]]))))))
						state[[i]] <- matrix(NA,dim(idx)[1],6)
						state[[i]][,3] <- (tmp1[[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+tmp1[[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
						state[[i]][,1] <- tmp1[[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
						state[[i]][,2] <- tmp1[[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
						state[[i]][,4] <- tmp1[[2]][1,4]
						state[[i]][,5] <- tmp1[[2]][1,6]
						state[[i]][,6] <- tmp1[[2]][1,6]
 						res[[i]] <- vector("list",4)
						res[[i]][[1]] <- tmp1[[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(state[[i]][,1:3]))
						state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
						res[[i]][[2]] <- idx
						res[[i]][[3]] <- c(seq(length(tmp[[1]]),length(tmp[[1]])*(length(tmp1[[1]])-1),length(tmp[[1]])),(length(tmp[[1]])*(length(tmp1[[1]])-1)+1):(length(tmp[[1]])*length(tmp1[[1]])))
						res[[i]][[4]] <- res1[[3]]
						Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]] <- cbind(Tt[[ifelse(is.element(a,Tt[[2]]),2,1)]],c(i,edge[i,4]))
						Tt[[ifelse(is.element(b,Tt[[2]]),2,1)]] <- cbind(Tt[[ifelse(is.element(b,Tt[[2]]),2,1)]],c(i,edge[i,5]))
				}
				}
			}
			res[[1]][[1]] <- exp(res[[1]][[1]])*state[[1]][,2]
			RES <- res[[1]]
			prob <- probtwo(res[[1]])
			RES[[4]] <- prob[[2]]
			prob <- prob[[1]]
			ll <- Tt
		}
	list(prob,RES[[4]],ll)
	}