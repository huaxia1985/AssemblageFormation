#this file includes the executable codes for the analysis of temporal patterns in species addition events of the New Zealand passerine assemblage. See details of the analysis in paper "Tracking the formation of species assemblages over time: phylogenetic reconstruction of patterns of colonization and speciation" by Xia Hua, Robert Lanfear, and Lindell Bromham

library(ape)
library(deSolve)
#load workspace file, workspace file includes all the functions and a list of 300 trees for the New Zealand passerine assemblage phylogenies. The workspace file is available from the authors.
load("workspace")

result <- vector("list",300)
# for the ith tree
for (i in c(1:300)) {
	dt <- 3 # set time step
	rr <- 0.5 # set prior parameter
	w <- rep(1,7)
	# enter tip numbers of each subtree
	native<-c("1","5","7","8","30","32","35","42","48","49","57","60","61","62","66","70","71","73","76","82","90","97","98","100","101")
	acanthisittidae<-c("1","35","70","71","97","100")
	anthus<-c("6","7")
	corvus<-as.character(18:31)
	gerygone<-as.character(38:47)
	megalurus<-as.character(54:58)
	mohoua<-as.character(60:62)
	petroica<-as.character(73:80)
	rhipidura<-as.character(88:92)
	zoserops<-"101"
	hirundo<-"49"
	turnagra<-"98"
	anthornis_pro<-c("5","82","9","10","11","86","87")
	callaeatidae<-c("8","32","48")
	notiomystis<-c("66")
	callaeatidae_notiomystis<-c("8","32","48","66")
	phy <- drop.tip (tree[[i]],c("38","58"))
	tot <- 0
	parent <- min(which(phy$edge[,1]==101))
	while (phy$edge[parent,2]>99) {
		tot <- tot+phy$edge.length[parent]
		parent <- min(which(phy$edge[,1]==phy$edge[parent,2]))
	}
	tot <- tot+phy$edge.length[parent]
	N <- length(phy$tip.label)
	#extract subtrees
	native <- which(is.element(phy$tip.label,native))
	acanthisittidae <- extractcom(phy,which(is.element(phy$tip.label,acanthisittidae)),category=2,N,native)
	anthus <- extractcom(phy,which(is.element(phy$tip.label,anthus)),category=2,N,native)
	corvus <- extractcom(phy,which(is.element(phy$tip.label,corvus)),category=1,N,native)
	gerygone <- extractcom(phy,which(is.element(phy$tip.label,gerygone)),category=1,N,native)
	megalurus <- extractcom(phy,which(is.element(phy$tip.label,megalurus)),category=1,N,native)
	mohoua <- extractcom(phy,which(is.element(phy$tip.label,mohoua)),category=2,N,native)
	petroica <- extractcom(phy,which(is.element(phy$tip.label,petroica)),category=1,N,native)
	rhipidura <- extractcom(phy,which(is.element(phy$tip.label,rhipidura)),category=1,N,native)
	anthornis_pro <- extractcom(phy,which(is.element(phy$tip.label,anthornis_pro)),category=1,N,native)
	callaeatidae_notiomystis <- extractcom(phy,which(is.element(phy$tip.label,callaeatidae_notiomystis)),category=2,N,native)
	zoserops <- extractcom(phy,which(is.element(phy$tip.label,zoserops)),category=3,N,native)
	hirundo <- extractcom(phy,which(is.element(phy$tip.label,hirundo)),category=3,N,native)
	turnagra <- extractcom(phy,which(is.element(phy$tip.label,turnagra)),category=3,N,native)
	#add known age of colonization
	zoserops[4] <- 0.0002
	hirundo[5] <- 0.0001
	# get tip branch lengths of the tree
	addtime <- c(acanthisittidae[acanthisittidae[,6]==1,4],acanthisittidae[acanthisittidae[,7]==1,5],
				anthus[anthus[,6]==1,4],anthus[anthus[,7]==1,5],
				corvus[corvus[,6]==1,4],corvus[corvus[,7]==1,5],
				gerygone[gerygone[,6]==1,4],gerygone[gerygone[,7]==1,5],
				megalurus[megalurus[,6]==1,4],megalurus[megalurus[,7]==1,5],
				mohoua[mohoua[,6]==1,4],mohoua[mohoua[,7]==1,5],
				petroica[petroica[,6]==1,4],petroica[petroica[,7]==1,5],
				rhipidura[rhipidura[,6]==1,4],rhipidura[rhipidura[,7]==1,5],
				anthornis_pro[anthornis_pro[,6]==1,4],anthornis_pro[anthornis_pro[,7]==1,5],
				callaeatidae_notiomystis[callaeatidae_notiomystis[,6]==1,4],callaeatidae_notiomystis[callaeatidae_notiomystis[,7]==1,5],max(turnagra[4:5]))
	# function to calculate the likelihood as the product of the likelhood of each subtree
	posterior <- function(pars) {
		out <- try(likcal(pars,acanthisittidae,category=2,Ngenus=NA,is.AB=F,is.AB.root=T,NA)+
		likcal(pars,anthus,category=1,Ngenus=39,is.AB=T,is.AB.root=F,1)+
		likcal(pars,corvus,category=1,Ngenus=42,is.AB=F,is.AB.root=F,13)+
		likcal(pars,gerygone,category=1,Ngenus=18,is.AB=F,is.AB.root=F,8)+
		likcal(pars,megalurus,category=1,Ngenus=43,is.AB=F,is.AB.root=F,3)+
		likcal(pars,mohoua,category=2,Ngenus=NA,is.AB=F,is.AB.root=F)+
		likcal(pars,petroica,category=1,Ngenus=8,is.AB=F,is.AB.root=F,6)+
		likcal(pars,rhipidura,category=1,Ngenus=41,is.AB=F,is.AB.root=F,4)+
		likcal(pars,anthornis_pro,category=1,Ngenus=6,is.AB=F,is.AB.root=F,5)+
		likcal(pars,callaeatidae_notiomystis,category=2,Ngenus=NA,is.AB=F,is.AB.root=F,NA)+
		likcal(pars,zoserops,category=3,Ngenus=NA,is.AB=T,is.AB.root=F,NA)+
		likcal(pars,hirundo,category=3,Ngenus=NA,is.AB=T,is.AB.root=F,NA)+
		likcal(pars,turnagra,category=3,Ngenus=NA,is.AB=F,is.AB.root=T,NA)
		if ((inherits(out, "try-error") || (!is.finite(out))||(is.na(out))))
     	 100000
   		else -out
	}
	# calculate the maximum likelihood estimate as the initial value of MCMC
	x.init <- nlm(posterior,p=rep(0.1,7),hessian=F)$estimate
	# function to calculate the posterior probability of the tree
	posterior <- function(pars) {
		out <- try(likcal(pars,acanthisittidae,category=2,Ngenus=NA,is.AB=F,is.AB.root=T,NA)+
		likcal(pars,anthus,category=1,Ngenus=39,is.AB=T,is.AB.root=F,1)+
		likcal(pars,corvus,category=1,Ngenus=42,is.AB=F,is.AB.root=F,13)+
		likcal(pars,gerygone,category=1,Ngenus=18,is.AB=F,is.AB.root=F,8)+
		likcal(pars,megalurus,category=1,Ngenus=43,is.AB=F,is.AB.root=F,3)+
		likcal(pars,mohoua,category=2,Ngenus=NA,is.AB=F,is.AB.root=F)+
		likcal(pars,petroica,category=1,Ngenus=8,is.AB=F,is.AB.root=F,6)+
		likcal(pars,rhipidura,category=1,Ngenus=41,is.AB=F,is.AB.root=F,4)+
		likcal(pars,anthornis_pro,category=1,Ngenus=6,is.AB=F,is.AB.root=F,5)+
		likcal(pars,callaeatidae_notiomystis,category=2,Ngenus=NA,is.AB=F,is.AB.root=F,NA)+
		likcal(pars,zoserops,category=3,Ngenus=NA,is.AB=T,is.AB.root=F,NA)+
		likcal(pars,hirundo,category=3,Ngenus=NA,is.AB=T,is.AB.root=F,NA)+
		likcal(pars,turnagra,category=3,Ngenus=NA,is.AB=F,is.AB.root=T,NA)+
		sum(dexp(pars, rr, log=TRUE)),silent=T)
		if ((inherits(out, "try-error") || (!is.finite(out))||(is.na(out))))
     	 -100000
   	else out
	}
	# start MCMC
	y.init <- posterior(x.init)
	if ( y.init == -100000 )
    stop("Starting point must have finite probability")
    output <- matrix(NA,100,length(x.init)+1)
	for ( ii in seq_len(100) ) {
      for (j in seq_along(x.init)) {
      	z <- y.init - rexp(1)
      	L <- x.init
      	R <- x.init
      	while (L[j]<0)
      	   {u <- runif(1) * w[j]
  			L[j] <- L[j] - u
  			R[j] <- R[j] + (w[j]-u)}
  		while ( posterior(L) > z )
    		{L[j] <- L[j] - w[j]}
  		while ( posterior(R) > z )
    		{R[j] <- R[j] + w[j]}
    	L<-max(0,L[j])
    	R<-R[j]
    	xs <- x.init
    	repeat {
    		xs[j] <- runif(1, L, R)
    		ys <- posterior(xs)
    		if ( ys > z )
      		break
    		if ( xs[j] < x.init[j] )
      		L <- xs[j]
    		else
      		R <- xs[j]
      	}
      	x.init <- xs
      	y.init <- ys
      }
      output[ii,] <- c(x.init,y.init)
  	}
  	for (j in seq_along(x.init)) {
	w[j] <- diff(quantile(output[,j],c(0.025, 0.975)))
	}
	output <- matrix(NA,2000,length(x.init)+1)
	for ( ii in seq_len(2000) ) {
      for (j in seq_along(x.init)) {
      	z <- y.init - rexp(1)
      	L <- x.init
      	R <- x.init
      	while (L[j]<0)
      	   {u <- runif(1) * w[j]
  			L[j] <- L[j] - u
  			R[j] <- R[j] + (w[j]-u)}
  		while ( posterior(L) > z )
    		{L[j] <- L[j] - w[j]}
  		while ( posterior(R) > z )
    		{R[j] <- R[j] + w[j]}
    	L<-max(0,L[j])
    	R<-R[j]
    	xs <- x.init
    	repeat {
    		xs[j] <- runif(1, L, R)
    		ys <- posterior(xs)
    		if ( ys > z )
      		break
    		if ( xs[j] < x.init[j] )
      		L <- xs[j]
    		else
      		R <- xs[j]
      	}
      	x.init <- xs
      	y.init <- ys
      }
      output[ii,] <- c(x.init,y.init)
    }
	  coltest <- matrix(NA,1000,3)
	  addtest <- matrix(NA,1000,3)
	  colprob <- matrix(NA,1000,floor(tot/dt))
	  colnllprob <- matrix(NA,1000,floor(tot/dt))
	  addnllprob <- matrix(NA,1000,floor(tot/dt))
	  # for the jth set of parameter values
	  for (j in 1:1000) {
	  	x.init<-as.numeric(output[j+1000,1:7])
	  	# calculate null expectations for tip branch length
      	addnll <- addcal(states=c(22/25,0,3/25,1-1/22,1,1-1/3),pars=x.init,tot=tot,dt=dt)
      	# calculate null expectations for colonization times
      	colnll <- colnllcal(states=c(16/19,0,3/19,1-1/16,1,1-1/3),pars=x.init,tot=tot,dt=dt)
      	# infer probability distribution of colonization time for the tree
      	coltime <- vector("list",10)
      	coltime[[1]] <- colcal(x.init,anthus,category=1,Ngenus=39,is.AB=T,dt,1)
        coltime[[2]] <- colcal(x.init,corvus,category=1,Ngenus=42,is.AB=F,dt,13)
		coltime[[3]] <- colcal(x.init,gerygone,category=1,Ngenus=18,is.AB=F,dt,8)
		coltime[[4]] <- colcal(x.init,megalurus,category=1,Ngenus=43,is.AB=F,dt,3)
		coltime[[5]] <- colcal(x.init,mohoua,category=2,Ngenus=NA,is.AB=F,dt,NA)
		coltime[[6]] <- colcal(x.init,petroica,category=4,Ngenus=8,is.AB=F,dt,6)
		coltime[[7]] <- colcal(x.init,rhipidura,category=1,Ngenus=41,is.AB=F,dt,4)
		coltime[[8]] <- colcal(x.init,anthornis_pro,category=4,Ngenus=6,is.AB=F,dt,5)
		coltime[[9]] <- colcal(x.init,callaeatidae_notiomystis,category=2,Ngenus=NA,is.AB=F,dt,NA)
		coltime[[10]] <- colcal(x.init,turnagra,category=3,Ngenus=NA,is.AB=F,dt,NA)
		for (ii in 1:10) {
			coltime[[ii]][[1]][coltime[[ii]][[1]]<0] <- 0
			coltime[[ii]][[1]] <- coltime[[ii]][[1]]/sum(coltime[[ii]][[1]])
		}
		coltime[[5]][[2]] <- apply(coltime[[5]][[2]],1,setprobthree,coltime[[5]][[3]])
		coltime[[6]][[2]] <- apply(coltime[[6]][[2]],1,setprobtwo,coltime[[6]][[3]],dt,coltime[[6]][[2]])
		coltime[[8]][[2]] <- apply(coltime[[8]][[2]],1,setprobtwo,coltime[[8]][[3]],dt,coltime[[8]][[2]])
		coltime[[9]][[2]] <- apply(coltime[[9]][[2]],1,setprobfour,coltime[[9]][[3]])
		prob <- numeric(floor(tot/dt))
		for (ii in c(1,2,3,4,7,10)) {
			prob[1:length(coltime[[ii]][[1]])] <- prob[1:length(coltime[[ii]][[1]])]+coltime[[ii]][[1]]
		}
		for (ii in c(5,6,8,9)) {
			prob1 <- numeric(length(colnll))
			tmp <- as.numeric(t(matrix(rep(coltime[[ii]][[1]][1:dim(coltime[[ii]][[2]])[2]],dim(coltime[[ii]][[2]])[1]),dim(coltime[[ii]][[2]])[2],dim(coltime[[ii]][[2]])[1])))
			prob1 <- sapply(c(1:length(colnll)),function (z) sum(tmp[coltime[[ii]][[2]]==z]))
			prob <- prob+prob1/sum(prob1)
		}
		prob[1] <- prob[1]+2
		prob <- prob/12
		colprob[j,]<- prob
		colnllprob[j,] <- colnll
		addnllprob[j,] <- addnll
		# test if the observed tip branch lengths are random draws from the null expectations
		addtest[j,] <- addsigtest(addnll,addtime,dt,n=23)
		# test if the inferred colonization times are random draws from the null expectations
 	  	coltest[j,] <- colsigtest(colnll,prob)
 	  }
 	  # results include: parameter values, clade size test, colonisation time test, inferred colonisation times, null expectations for colonisation times, tip length test, observed tip branch length, null expectation for tip branch length
	result[[i]] <- list(ParameterList=output,ColonisationResult=coltest,InferredColonisationTime=colprob,NullColonisationTime=colnllprob,TipLengthResult=addtest,TipLength=addtime,NullTipLength=addnllprob)
}
