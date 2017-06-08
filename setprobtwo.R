setprobtwo <- function (sample,ll,dt,RES) {
	if (length(sample)==2) {
		if (is.matrix(ll[[1]])) {
			b <- 0
			t1 <- 0
			for (i in 1:dim(ll[[1]])[2]) {
				if (ll[[1]][2,i]<t1) {
					idx <- b+1
 				} else {
					a <- floor(ll[[1]][2,i]/dt)
					idx <- c(b+1:(a+1))
					b <- b+a+1
					t1 <- dt-ll[[1]][2,i]+a*dt
				}
				if (is.element(sample[1],idx)) {
					if (i!=1) {
					sample[1] <- floor(sum(ll[[1]][2,1:(i-1)])/dt)+which(idx==sample[1])
					}
					break
				}
			}
		} 
		if (is.matrix(ll[[2]])) {
			b <- 0
			t1 <- 0
			for (i in 1:dim(ll[[2]])[2]) {
				if (ll[[2]][2,i]<t1) {
					idx <- b+1
 				} else {
					a <- floor(ll[[2]][2,i]/dt)
					idx <- c(b+1:(a+1))
					b <- b+a+1
					t1 <- dt-ll[[2]][2,i]+a*dt
				}
				if (is.element(sample[2],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(ll[[2]][2,1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
		}
		out<-sample
	} else {
		if (sample[1]==0) {
			if (is.matrix(ll[[1]])) {
			b <- 0
			t1 <- 0
			for (i in 1:dim(ll[[1]])[2]) {
				if (ll[[1]][2,i]<t1) {
					idx <- b+1
 				} else {
					a <- floor(ll[[1]][2,i]/dt)
					idx <- c(b+1:(a+1))
					b <- b+a+1
					t1 <- dt-ll[[1]][2,i]+a*dt
				}
				if (is.element(sample[2],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(ll[[1]][2,1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
			} 
			if (is.matrix(ll[[2]])) {
			b <- 0
			t1 <- 0
			for (i in 1:dim(ll[[2]])[2]) {
				if (ll[[2]][2,i]<t1) {
					idx <- b+1
 				} else {
					a <- floor(ll[[2]][2,i]/dt)
					idx <- c(b+1:(a+1))
					b <- b+a+1
					t1 <- dt-ll[[2]][2,i]+a*dt
				}
				if (is.element(sample[3],idx)) {
					if (i!=1) {
					sample[3] <- floor(sum(ll[[2]][2,1:(i-1)])/dt)+which(idx==sample[3])
					}
					break
				}
			}
			}
		out<-sample[2:3]
		} else {
			sample[1] <- sample[1]+floor(ifelse(is.matrix(ll[[2]]),sum(ll[[2]][2,]),ll[[2]][2])/dt)
			if (sample[2]<max(RES[,2])) {
				if (is.matrix(ll[[1]])) {
				b <- 0
				t1 <- 0
				for (i in 1:dim(ll[[1]])[2]) {
					if (ll[[1]][2,i]<t1) {
						idx <- b+1
 					} else {
						a <- floor(ll[[1]][2,i]/dt)
						idx <- c(b+1:(a+1))
						b <- b+a+1
						t1 <- dt-ll[[1]][2,i]+a*dt
					}
					if (is.element(sample[2],idx)) {
						if (i!=1) {
						sample[2] <- floor(sum(ll[[1]][2,1:(i-1)])/dt)+which(idx==sample[2])
						}
						break
					}
				}
				}
				out<-sample[1:2] 
			}
			if (sample[3]<max(RES[,3])) {
				if (is.matrix(ll[[2]])) {
				b <- 0
				t1 <- 0
				for (i in 1:dim(ll[[2]])[2]) {
					if (ll[[2]][2,i]<t1) {
						idx <- b+1
 					} else {
						a <- floor(ll[[2]][2,i]/dt)
						idx <- c(b+1:(a+1))
						b <- b+a+1
						t1 <- dt-ll[[2]][2,i]+a*dt
					}
					if (is.element(sample[3],idx)) {
						if (i!=1) {
						sample[3] <- floor(sum(ll[[2]][2,1:(i-1)])/dt)+which(idx==sample[3])
						}
						break
					}
				}
				}
				out<-sample[c(1,3)]
			}
			if (sample[2]==max(RES[,2])&&sample[3]==max(RES[,3])) {
				out<-c(sample[1],sample[1])
			}
		}
}
out
}