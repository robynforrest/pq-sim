#******************************************************
#	Programmers: Robyn Forrest, Jon Schnute and Rowan Haigh
#	Authors: Jon Schnute, Robyn Forrest, Rowan Haigh: Pacific Biological Station
#	Project Name: Dealing with those pesky zeros: Simulation of Bernouilli-Dirichlet model for age-composition data
#	Date:	April 6, 2014
#	Version:  1.0
#	Comments:   Script file for batching up pq simulation in R and estimation in ADMB R 
#	
#******************************************************/
rm(list=ls(all=T))
graphics.off()
require(PBSmodelling)
source("read.admb.R")	  #Martell code for reading rep file

#_____________________________________________
#GUI FUNCTIONS
#_____________________________________________
restart <- function() {source("pqBatch.r")}  
createWin("pqWin.txt");

#_____________________________________________
#TRANSFORMATION FUNCTIONS
#_____________________________________________
#Logit
logit <- function(p) {log(p/(1-p))}		 #T2.1

#Inverse Logit
invlogit <- function(x) {exp(x)/(1+exp(x))}

#Get q from p 
p2q <- function(p) {
  getWinVal(scope="L");
      invlogit(a - b*logit(p))}

q2p <- function(qq) {
  getWinVal(scope="L");
  invlogit((a - logit(qq))/b)}

# Composition operator
# Converts vector v to proportions
conv2p <- function(v, rm.zero=TRUE) {
	if (rm.zero) v <- v[v>0 & !is.na(v)]
	y <- abs(v)/sum(abs(v)); 
	y; };

#_____________________________________________
#FUNCTIONS TO SIMULATE DATA
#_____________________________________________
pqSim <- function() {
	  getWinVal(scope="L");
	  age = 1:n;
	  Zi <- vector(length=n) #Total mortality at age
	  survi <- vector(length=n) #Eqm fished survival rate at age
	  Si <- vector(length=n)	#Eqm fished survivorship at age
	  
	  #Selectivity
	  betai <- plogis(age,ah,gh)
	  
	  #Total mortality
	  Zi <- M + betai*Fish
		  
	  #Survival rate -- at eqm with constant mortality
	 survi <- exp(-Zi)
	  
	  #Eqm survivorship
	  Si[1]<-1
	  for(j in 2:n) Si[j] <- Si[j-1]*survi[j-1]
            Si[n] = Si[n]/(1 - survi[n])
	  
	  Ri = rep(1,n);
	  pvec = Si*betai*Ri;    
	  pvec = pvec / sum(pvec);	 
	  a = logit(qtil) + b*logit(1/n); # pbar = 1/n	   T2.3

	  setWinVal(list(a=a))
	  qvec=p2q(pvec);  #T2.6
	  age <<- age; Si <<- Si; betai <<- betai; 
	  pvec <<- pvec; qvec <<- qvec;

	  #print(pvec)
	  #print(qvec)
	  return(cbind(age,pvec,qvec,Si,betai)); 
  }#End Function pqSim

# Simulation of data y (with ns samples)
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
  return(yvec); };

# Simulate ns sample vectors y'
xyGen <- function() {
  getWinVal(scope="L"); 
  pqSim();
  ymat <- matrix(0,nrow=n,ncol=ns);
  for (i in (1:ns)) {xySim(); ymat[,i] = yvec; };
  ymat <<- ymat; return(ymat) }

#_____________________________________________
#PLOTTING FUNCTIONS
#_____________________________________________
#SIMULATION DATA
pqPairs <- function() {
  pairs(cbind(age,pvec,qvec,Si,betai)); };

pqPlot <- function() {
  getWinVal(scope="L"); # for qtil
  z = order(pvec)
  plot(pvec[z],qvec[z]); lines(pvec[z],qvec[z],col="blue",lwd=2);
  abline(h=qtil,col="red"); abline(v=1/n,col="red"); };

pqPlotFull <- function() {
  getWinVal(scope="L"); # for qtil
  minp = min(.001,pvec);
  maxp = max(.999,pvec);
  pp = seq(from=minp,to=maxp,length=1000);
  qq = p2q(pp);
  plot(pp,qq,type="l",col="green4",lwd=2); points(pvec,qvec,pch=21,col="green4",bg="green");
  abline(h=qtil,col="red"); abline(v=1/n,col="red"); };

pqPlotLogit <- function() {
  getWinVal(scope="L"); # for qtil
  minp = min(.001,pvec);
  maxp = max(.999,pvec);
  pp = seq(from=minp,to=maxp,length=1000);
  qq = p2q(pp);
  lp = logit(pp); lq =logit(qq);
  plot(lp,lq,type="l",col="green4",lwd=2);
  points(logit(pvec),logit(qvec),pch=21,col="green4",bg="green");
  abline(h=logit(qtil),col="red"); abline(v=logit(1/n),col="red");
  abline(a=0,b=-1,col="blue") };

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

plotFit <- function(age,sim,est,ii,typ,nam,onePg=FALSE, Xlab="", Ylab=""){
	# barplot(sim, names.arg=paste(age), las=1,col="gray", xlab="Age", ylab="Proportion", ylim=c(0,0.15)) 
	plot(age,sim, type=typ, las=1,col=1, xlab=Xlab, ylab=Ylab, lwd=2, xaxt=ifelse(onePg,"n","s"), yaxt=ifelse(onePg,"n","s"), ylim=c(0, 1.1*max(sim,est)))#, ylim=c(0,0.15)
	if (onePg) {
		axis(1,labels=FALSE,tcl=0.5)
		axis(2,labels=FALSE,tcl=0.5)
		if (par()$mfg[1]==par()$mfg[3]) axis(1,tick=FALSE,mgp=c(1.75,0.5,0),las=1)
		if (par()$mfg[2]==1) axis(2,tick=FALSE,mgp=c(1.75,0.5,0),las=1)
	}
	lines(age,est,type="o",pch=20,cex=2,col=2)
	legend("topright", legend=c("Sim","Est"),col=c(1,2), pch=c(15,19), bty="n", cex=1.1,inset=0.05)
	mtext(paste("Sample",ii), side=3, outer=F, line=-1.25, adj=ifelse(onePg&&typ!="h",0.05,0.5))
	mtext("Age", side=1, outer=T, line=ifelse(onePg,2,-0.05), cex=1.3)
	mtext("Proportion", side=2, outer=T, line=ifelse(onePg,2.5,0.5), cex=1.3)
	mtext(nam, side=3, outer=T, line=ifelse(onePg,0.5,-1), cex=1.3)
}

panel.hist <- function(x, ...){
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="tan", ...)
  }
 
 plotFitPairs <- function(Pars) {
	 pairs(Pars,pch=20,upper.panel=panel.smooth,diag.panel=panel.hist, lower.panel=panel.smooth)
  }

#_____________________________________________
#ADMB FUNCTIONS
#_____________________________________________
#Function to rewrite the data and pin file for pq.exe
write_dat_pin=function(yobs)
{
	getWinVal(scope="L");
	dfile="pq.dat"
	write("#Parameters and simulated data from p-q simulation",dfile,1)
	write("#n Maximum age class",dfile,1,append=T)
	write(n,dfile,1,append=T)
	write("#rbar Average recruitment ",dfile,1,append=T)
	write(rbar,dfile,1,append=T)
	write("#Fish Constant fishing mortality -- will be estimated eventually ",dfile,1,append=T)
	write(Fish,dfile,1,append=T)
	write("#yObs Simulated age proportions from Bernoulli Dirichlet distribution",dfile,1,append=T)
	write(yobs,dfile,1,append=T)
	write(" #debug	Switches on cout statements",dfile,1,append=T)
	write(0,dfile,1,append=T)
	write("#eof",dfile,1,append=T)
	write(999,dfile,1,append=T)

	dfile="pq.pin"
	write("#Initial parameter values for logM, beta1, alpha, qtil, b and N",dfile,1)
	write(log(M),dfile,1,append=T)
	write(ah,dfile,1,append=T)
	write(gh,dfile,1,append=T)
	write(qtil,dfile,1,append=T)
	write(b,dfile,1,append=T)
	write(N,dfile,1,append=T)
}


#_____________________________________________
#This function is called from the GUI
#MPD
callADMB <- function() {
	getWinVal(scope="L");
	ymatEst <- matrix(0,nrow=n, ncol=ns+1)
	selEst  <- matrix(0,nrow=n, ncol=ns+1)
	parEst <-  matrix(0,nrow=ns, ncol=6)
	colnames(parEst) <- c("logM", "ah","gh","qtil","b","N")
	xyGen(); #generate ns samples 
	
	#Loop over ns samples
       	for(i in 1:ns) { 
		yprime <- ymat[,i]
		write_dat_pin(yprime);
        	system("pq.exe -maxfn 2000 -nox",show.output.on.console=F,invisible=F)
		
		#Read the report file from ADMB and put proportions at age and selectivity from this run into a matrix
		out<-readList("pq.rep")
		if(i==1) {
			ymatEst[,1] <- out$age
			selEst[,1]    <- out$age}

		ymatEst[,(i+1)] <- out$pvec
		selEst[,(i+1)] <- out$Betai
		parEst[i,] <-c(out$log_M, out$ah, out$gh, out$qtil, out$b, out$N)
	}
       parEst <<- parEst
       selEst <<- selEst
       ymatEst<<-ymatEst
  }

#_____________________________________________
#This function is called from the GUI
#MCMC
callADMBmc <- function() {
	getWinVal(scope="L");
		
	yprime <- xySim();
	write_dat_pin(yprime);
        admbcall <- paste("pq.exe -maxfn 2000 -mcmc", mcsamples, "-mcsave", mcthin, "-mcscale", mcsamples)
        system(admbcall,show.output.on.console=F,invisible=F)            
        system("pq.exe -mceval",show.output.on.console=F,invisible=F );
        
        parEstmc <<- read.table("pqpars.csv", header=T, sep=",")
        selEstmc  <<- read.table("pqselectivity.csv", header=F, sep=",")
	ymatEstmc <<- read.table("pqproportions.csv", header=F, sep=",")
  }


#Call the function to plot proportions at age
fitProp<-function() {
	oldpar=par(no.readonly=TRUE); on.exit(par(oldpar))
	getWinVal(scope="L");
	graphcount<-0
	if (onePage) {
		rc = .findSquare(ns)
		par(mfrow=rc, mar=c(0,0,0,0), oma=c(4,4,2,1)) # All graphs
	} else {
		par(mfrow=c(2,2), oma=c(2,2,1,1), mai=c(.35,.35,.3,.3)) #4 graphs
	}
	for(i in 1:ns) {
		graphcount<-graphcount+1
		plotFit(ymatEst[,1],ymat[,i], ymatEst[,(i+1)],i, "h", "Proportions-at-age",onePg=onePage, Xlab="", Ylab="")
		if(!onePage && graphcount==4) {
			if(ns>4){
				windows()
				par(mfcol=c(2,2), oma=c(2,3,1,1), mai=c(.35,.35,.3,.3))
				graphcount <-0} # end if
		}# end if
	} #end for
} #end function

#Call the function to plot selectivity at age
fitSel<-function() {
	oldpar=par(no.readonly=TRUE); on.exit(par(oldpar))
	getWinVal(scope="L");
	graphcount<-0
	if (onePage) {
		rc = .findSquare(ns)
		par(mfrow=rc, mar=c(0,0,0,0), oma=c(4,4,2,1)) # All graphs
	} else {
		par(mfrow=c(2,2), oma=c(2,2,1,1), mai=c(.35,.35,.3,.3)) #4 graphs
	}
	for(i in 1:ns) {
		plotFit(selEst[,1],betai, selEst[,(i+1)],i, "l", "Selectivity",onePg=onePage, Xlab="", Ylab="")
		if(!onePage && graphcount==4) {
			if(ns>4){
				windows()
				par(mfcol=c(2,2), oma=c(2,3,1,1), mai=c(.35,.35,.3,.3))
				graphcount <-0} # end if
		}# end if
	} #end for
} #end function

 fitPairs <- function() {
       plotFitPairs(parEst)
 }
 
 fitPropmc<-function() {
 	getWinVal(scope="L");
 	medY <- vector(length=n)
 	for(i in 1:n) medY[i] <- median(ymatEstmc[,i])
 	plotFit(1:n,yvec, medY,1, "h", "Proportions-at-age",onePg=TRUE, Xlab="Age", Ylab="Proportions")
 	mtext("MCMC Median fitted proportions", side=3, line=1, cex=1.5)
 	
 } #end function
 
 #Call the function to plot selectivity at age
 fitSelmc<-function() {
 	getWinVal(scope="L");
 	medSel <- vector(length=n)
 	for(i in 1:n) medSel[i] <- median(selEstmc[,i])
 	
 	matplot(1:n,t(selEstmc), type="l", xlab="Age", ylab="Proportion", col="darkgray", las=1)
 	lines(1:n,medSel, col=2,lwd=2)
 	lines(1:n,betai, col=1,lwd=2)
 	legend("topleft", legend=c("MCMC Samples","MCMC Median", "True"), col=c("darkgray",2:1), lty=1, lwd=2, bty="n")
 		
 } #end function
 
  fitPairsmc <- function() {
        plotFitPairs(parEstmc)
 }

