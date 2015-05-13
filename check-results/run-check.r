
library(gleeR)

load("A.RData")
for (j in 1:ncol(A)) {
	ind = which(is.na(A[,j]))
	A[ind,j] = 0
}

load("Acc.RData")

# comp = "prep-wt"; Data = cbind(Acc, A[,c(2,6,10,1,5,9)])
# comp = "oop-wt"; Data = cbind(Acc, A[,c(3,7,11,1,5,9)])
# comp = "triple-wt"; Data = cbind(Acc, A[,c(4,8,12,1,5,9)])
# comp = "prep-triple"; Data = cbind(Acc, A[,c(2,6,10,4,8,12)])
comp = "oop-triple"; Data = cbind(Acc, A[,c(3,7,11,4,8,12)])

# number of replicates
nA = 3
nB = 3

# options
fit_type = "cubic"
num_iter = 10000
num_digits = 4

if (!data_ok(Data,nA,nB)) {
	msg = paste("check input data! spreadsheet must have",
			"(1) the right number of columns",
			"(2) positive finite values",
			"\n",sep="\n")
	stop(msg)
}

Prot = as.character(Data[,1])
A = as.matrix(Data[,1+(1:nA)])
B = as.matrix(Data[,1+nA+(1:nB)])

m = fit_model(A, B, fit_type)
model_fit_plots(m, outfile=paste0(comp,".fitplots.png"))
stn_pval = calc_stn_pval(A, B, m, num_iter)
stn_pval_plots(stn_pval, outfile=paste0(comp,".stn-pval.png"))
T = diff_exp_table(stn_pval, Prot, num_digits)
write.csv(T, file=paste0(comp,".glee-diffexp.csv"), row.names=FALSE)
