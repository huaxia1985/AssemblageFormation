setprobthree <- function (sample,ll) {
	if (sample[3]==0) {
		t1 <- sample[4]
		t2 <- sample[5]
		if (sample[2]<ll[2]) {
			t3 <- sample[2]
		} else {
			t3 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
		}
	} else if (sample[4]<ll[4]) {
		t1 <- sample[4]
		if (sample[1]==0) {
			t3 <- sample[2]
			t2 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
		} else {
			if (sample[2]<ll[2]) {
				t3 <- sample[2]
			} else {
				t3 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[3]<ll[3]) {
				t2 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
			} else {
				t2 <- sample[1]+ifelse(ll[5]==3,0,ll[5]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
		}
	} else if (sample[5]<ll[5]) {
		t2 <- sample[5]
		if (sample[1]==0) {
			t3 <- sample[2]
			t1 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
		} else {
			if (sample[2]<ll[2]) {
				t3 <- sample[2]
			} else {
				t3 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[3]<ll[3]) {
				t1 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
			} else {
				t1 <- sample[1]+ifelse(ll[5]==3,0,ll[5]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
		}	
	} else {
		if (sample[1]==0) {
			t3 <- sample[2]
			t1 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
		} else {
			if (sample[2]<ll[2]) {
				t3 <- sample[2]
			} else {
				t3 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[3]<ll[3]) {
				t1 <- sample[3]+ifelse(ll[5]==3,0,ll[5]-2)		
			} else {
				t1 <- sample[1]+ifelse(ll[5]==3,0,ll[5]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
		}	
		t2 <- t1
	}
    c(t1,t2,t3)
}