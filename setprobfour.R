setprobfour <- function (sample,ll) {
	if (sample[5]==0) {
		t1 <- sample[6]
		t2 <- sample[7]
		if (sample[3]==0) {
			t3 <- sample[4]
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else if (sample[3]<ll[3]) {
			t3 <- sample[3]+ifelse(ll[4]==3,0,ll[4]-2)
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else {
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[1]+ifelse(ll[4]==3,0,ll[4]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}

		}
	} else if (sample[6]<ll[7]) {
		t1 <- sample[6]
		if (sample[3]==0) {
			t3 <- sample[4]
			t2 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else if (sample[3]<ll[3]) {
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[3]+ifelse(ll[4]==3,0,ll[4]-2)
			}
			if (sample[5]<ll[5]) {
				t2 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t2 <- sample[3]+ifelse(ll[7]==3,0,ll[7]-2)+ifelse(ll[5]==2,0,ll[5]-2)
			}
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else {
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[1]+ifelse(ll[4]==3,0,ll[4]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
			if (sample[5]<ll[5]) {
				t2 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t2 <- sample[1]+ifelse(ll[3]==2,0,ll[3]-2)+ifelse(ll[5]==2,0,ll[5]-2)+ifelse(ll[7]==3,0,ll[7]-2)
			}
		} 	
	} else if (sample[7]<ll[7]) {
		t2 <- sample[7]
		if (sample[3]==0) {
			t3 <- sample[4]
			t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}	
		} else if (sample[3]<ll[3]) {
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[3]+ifelse(ll[4]==3,0,ll[4]-2)
			}
			if (sample[5]<ll[5]) {
				t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t1 <- sample[3]+ifelse(ll[7]==3,0,ll[7]-2)+ifelse(ll[5]==2,0,ll[5]-2)
			}
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else {
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[1]+ifelse(ll[4]==3,0,ll[4]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
			if (sample[5]<ll[5]) {
				t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t1 <- sample[1]+ifelse(ll[3]==2,0,ll[3]-2)+ifelse(ll[5]==2,0,ll[5]-2)+ifelse(ll[7]==3,0,ll[7]-2)
			}
		}	
	} else {
		if (sample[3]==0) {
			t3 <- sample[4]
			t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else if (sample[3]<ll[3]) {
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[3]+ifelse(ll[4]==2,0,ll[4]-2)
			}
			if (sample[5]<ll[5]) {
				t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t1 <- sample[3]+ifelse(ll[7]==3,0,ll[7]-2)+ifelse(ll[5]==2,0,ll[5]-2)
			}
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
		} else {
			if (sample[2]<ll[2]) {
				t4 <- sample[2]
			} else {
				t4 <- sample[1]+ifelse(ll[2]==3,0,ll[2]-2)
			}
			if (sample[4]<ll[4]) {
				t3 <- sample[4]
			} else {
				t3 <- sample[1]+ifelse(ll[4]==3,0,ll[4]-2)+ifelse(ll[3]==2,0,ll[3]-2)
			}
			if (sample[5]<ll[5]) {
				t1 <- sample[5]+ifelse(ll[7]==3,0,ll[7]-2)		
			} else {
				t1 <- sample[1]+ifelse(ll[3]==2,0,ll[3]-2)+ifelse(ll[5]==2,0,ll[5]-2)+ifelse(ll[7]==3,0,ll[7]-2)
			}
		}
		t2 <- t1
	}
	c(t1,t2,t3,t4)
}