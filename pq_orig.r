# Simulation of p and q

require(PBSmodelling)

logit <- function(p) {log(p/(1-p))}		 #T2.1

invlogit <- function(x) {exp(x)/(1+exp(x))}

p2q <- function(p) {
  getWinVal(scope="L");
      invlogit(a - b*logit(p))}

q2p <- function(qq) {
  getWinVal(scope="L");
  invlogit((a - logit(qq))/b)}

pqSim <- function() {
  getWinVal(scope="L");
  age = 1:n;
  Si = exp(-Z*(age-1)); Si[n] = Si[n]/(1 - exp(-Z))
  betai = rep(1,n); 
  ii = 1:(A-1);
  betai[ii] = 1 - (1-beta1)*((A-ii)/(A-1))^alpha;
  Ri = rep(1,n);
  pvec = Si*betai*Ri; pvec = pvec / sum(pvec);
  a = logit(qtil) + b*logit(1/n); # pbar = 1/n	   T2.3

  
  setWinVal(list(a=a))
  qvec=p2q(pvec);  #T2.6
   cat("pvec\n")
   print(pvec)
  
  cat("qvec\n")
   print(qvec)

  age <<- age; Si <<- Si; betai <<- betai; 
  pvec <<- pvec; qvec <<- qvec;
  return(cbind(age,pvec,qvec,Si,betai)); }

pqPairs <- function() {
  pairs(cbind(age,pvec,qvec,Si,betai)); };

pqPlot <- function() {
  getWinVal(scope="L"); # for qtil
  plot(pvec,qvec); lines(pvec,qvec,col="blue",lwd=2);
  abline(h=qtil,col="red"); abline(v=1/n,col="red"); };

pqPlotFull <- function() {
  getWinVal(scope="L"); # for qtil
  minp = min(.001,pvec);
  maxp = max(.999,pvec);
  pp = seq(from=minp,to=maxp,length=1000);
  qq = p2q(pp);
  plot(pp,qq,type="l"); points(pvec,qvec,col="green");
  abline(h=qtil,col="red"); abline(v=1/n,col="red"); };

pqPlotLogit <- function() {
  getWinVal(scope="L"); # for qtil
  minp = min(.001,pvec);
  maxp = max(.999,pvec);
  pp = seq(from=minp,to=maxp,length=1000);
  qq = p2q(pp);
  lp = logit(pp); lq =logit(qq);
  plot(lp,lq,type="l");
  points(logit(pvec),logit(qvec),col="green");
  abline(h=logit(qtil),col="red"); abline(v=logit(1/n),col="red");
  abline(a=0,b=-1,col="blue") };

# Simulation of data y (with ns samples)

# Composition operator
# Converts vector v to proportions
conv2p <- function(v, rm.zero=TRUE) {
	if (rm.zero) v <- v[v>0 & !is.na(v)]
	y <- abs(v)/sum(abs(v)); 
	y; };

# Dirichlet random sample
#   Input:  ns = number of samples
#           N = effective sample size
#           pvec = probability vector (sums to 1)
#   Output: matrix (ns x length(pvec))
#           each row is a Dirichlet sample that sums to 1
rdirich <- function(ns,N,pvec) {
	pvec <- conv2p(pvec); # forces sum to 1
	np <- length(pvec);
	yvec <- rgamma(ns*np,shape=N*pvec,scale=N);
	y1 <- matrix(yvec,nrow=ns,ncol=np,byrow=TRUE);
	ys <- apply(y1,1,"sum");
	y2 <- sweep(y1,1,ys,"/");
	if (ns==1) as.vector(y2) else y2 
}

# Simulate x, y' once
xySim <- function() {
  getWinVal(scope="L"); 
  pqSim();
  uvec = runif(n);
  zvec = (uvec>=qvec) # zvec is F when xi=0
  xvec = as.numeric(zvec) # vector x 
  pp = conv2p(pvec*xvec,rm.zero=TRUE)  # zeros not included
  yy = rdirich(1,N,pp)
  yvec = rep(0,n); 
  yvec[zvec] = yy;
  xvec <<- xvec; 
  yvec <<- yvec; 
  return(rbind(xvec,yvec)); };

# Simulate ns sample vectors y'
xyGen <- function() {
  getWinVal(scope="L"); 
  pqSim();
  ymat <- matrix(0,nrow=n,ncol=ns);
  for (i in (1:ns)) {xySim(); ymat[,i] = yvec; };
  ymat <<- ymat; return(ymat) }

yBub <- function() {
  getWinVal(scope="L");
  yB <- cbind(pvec,qvec,ymat);
  plotBubbles(yB,clrs=c("black","red","red"),size=sz,
    powr=powr,prettyaxis=TRUE); };

xBub <- function() {
  getWinVal(scope="L");
  xmat = max(qvec)*matrix(as.numeric(ymat==0),nrow=n,ncol=ns,byrow=F);
  xB <- cbind(pvec,qvec,xmat);
  plotBubbles(xB,clrs=c("black","red","red"),size=sz,
    powr=powr,prettyaxis=TRUE); };

restart <- function() {source("pq.r")}  
  
createWin("pqWin.txt");
  