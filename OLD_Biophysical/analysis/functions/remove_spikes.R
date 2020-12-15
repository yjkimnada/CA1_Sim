
remove.spikes <- function(X, spt, dt){
	i.ths <- spt$st - spt$delay/dt
	n.cut <- spt$dur/dt
	## i.th: indices of the first points to cut
	## n.cut: number of points to cut
	
	if (n.cut<1) stop("n.cut must be greater than 0")
	if (is.vector(X)) X <- matrix(X, 1, length(X))
	nx <- ncol(X)
	for (i.th in i.ths){
		# spt=10; delay=-2; dur=5; [10-2]-[10-2+5]; cut 9-10-11-12
		X[,(max(i.th,1)):(min(i.th + n.cut -1, nx))] <- X[,max(i.th,1)]
	}
	X
}
