//******************************************************
//	Programmer: Robyn Forrest
//	Authors: Jon Schnute, Robyn Forrest, Rowan Haigh: Pacific Biological Station
//	Project Name: Dealing with those pesky zeros: Simulation of Bernouilli-Dirichlet model for age-composition data
//	Date:	April 6, 2014
//	Version:  1.0
//	Comments:   Called from R script pqBatch.r
//	
//******************************************************/

DATA_SECTION
	  int nf; //function evaluations

	   //Read in parameters from data file
	   init_int n;	//Maximum age class				
	   init_number rbar; //Average recruitment
	   init_number Fish; //constant fishing mortality -- estimate eventually
	   init_vector yObs(1,n);
	   init_number ubM; //upper bound of log_M
	   init_int debug; //debugging flag
	   init_int eof;	  //end of file marker

	//Trap for incorrect dimensions in data file
	LOC_CALCS
		 cout<<"eof = "<<eof<<endl;
	   	  if(eof==999){
	   		cout<<"\n -- SUCCESSFULLY READ DATA FILE -- \n"<<endl;
	   	  }else{
	   		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	   	  }
	END_CALCS

	vector age(1,n);	
	vector xi(1,n);
	
	LOC_CALCS
	       age.fill_seqadd(1,1);
	       	if(debug) COUT(age);
	      
	      //Get xi :: Eq T4.1
	      for(int i=1; i<=n; i++) {
			if(yObs(i)) {
				xi(i) =1;
			}else xi(i) =0;
	       }// end for loop
	       	if(debug) COUT(xi);

	       nf=0;
	END_CALCS

PARAMETER_SECTION
	init_bounded_number log_M(-2.99,ubM,1);              //Natural mortality	      
	init_bounded_number ah(1,n,1);      //Age at 50% selectivity
	init_bounded_number gh(0.0001,2.,1);             //SD in selectivity
	init_bounded_number qtil(0.0001,0.999,1);	   //Parameter relating p to q
	init_bounded_number b(0.0001,3.,1);		 //Parameter relating p to q	 
	init_bounded_number N(5,500.,1);	         //Effective sample size for Dirichlet distribution -- doesn't work very well if estimating this
	
	number a;
	number M;
 	vector Si(1,n);
        vector Betai(1,n);
 	vector Ri(1,n);		
 	vector pvec(1,n);		
	vector qvec(1,n);
	vector pPrime(1,n);
	vector log_resid(1,n);

	objective_function_value f;

PROCEDURE_SECTION
	getQi(); 	//Calls getPi 
	getpPrime();
	calcObjectiveFunction();
	if (mceval_phase()) mcOutput();
		
FUNCTION calcObjectiveFunction
 {
       //Calculate negative log likelihood from T4.4 
       //PSEUDOCODE FOR LIKELIHOOD
       //like = sum(gammaln(N*pPrime) - (N*pPrime - 1)*log(yObs)) - sum((1 - xi)*log(qvec) + xi*log(1-qvec) - gammaln(N));
       
       dvariable like1;	   //First part of objective function - sum over pPrime >0
       dvariable like;
       
       //Build components of objective function
       dvar_vector Npp = N*pPrime;			//N*pPrime
       dvar_vector Nppminus = N*pPrime - 1.;	//N*pPrime - 1
       dvar_vector minusxi = -xi;
       dvar_vector minusqi = -qvec;
       minusxi +=1.; 							// 1-xi
       minusqi +=1.; 							// 1-qvec 
      
      //Calculate objective function
      like1=0.;
      //Calculate like1, the first part of objective function :: only for pPrime>0
      for(int i=1; i<=n;i++){
	   if(pPrime(i)>0) like1 += gammln(Npp(i)) - (Nppminus(i))*log(yObs(i));
      }
      //Calculate the objective function
      like = like1 - sum((minusxi)*log(qvec) + elem_prod(xi,log(minusqi)) - gammln(N));
      
      f = like;
      nf++;
  }// end function

FUNCTION void getpPrime()
  {
	  pPrime=elem_prod(xi,pvec);
	  pPrime/=sum(pPrime); //EqT4.2
  }

FUNCTION void getPi()
  {
		//Function to estimate proportions at age and sampled proportions at age. 		
		//This is an equilibrium model with the population in equilibrium with constant fishing mortality
                //This function is called by getQi()
		
		//Declare and initialize variables
		int j;
		dvar_vector surv(1,n); 
		
		dvariable aminus;
		Si.initialize(); Betai.initialize(); Ri.initialize(); pvec.initialize(); qvec.initialize();
               	
               	//1. Logistic selectivity
               	Betai = plogis<dvar_vector>(age, ah,gh);

		//2. Get fished survival rate at age 
		M=mfexp(log_M);
		for(j=1; j<=n; j++) surv(j)=mfexp(-M - Betai(j)*Fish); 
		
		//3. Get survivorship -- this is eqm vulnerable numbers so no need for selectivity when calculating pi
		Si(1)=1.;
		for(j=2; j<=n; j++)  Si(j)=Si(j-1)*surv(j-1);  
		Si(n) /=(1. - surv(n));
		
		if(debug) {
				COUT(Fish);
				COUT(M);
			        COUT(surv);
				COUT(Si);
				COUT(Betai);
				COUT(ah);
				COUT(gh);
			}
				
		//Recruits to each age class -- placeholder -- add time and anomalies later 
		Ri = rbar;	  //Eq T6.3
		
		//Get pvec -- proportions at age in the population
		pvec = elem_prod(elem_prod(Si,Betai),Ri); 
		pvec /= sum(pvec);   //EqT6.4	 -- normalise
			if(debug) COUT(pvec);
 }

FUNCTION void getQi()
  {
	dvariable invn;
	getPi();
	 //Get qvec
	 invn=1./n;
	 a = logit(qtil) + b*logit(invn);

	qvec=p2q(pvec, a, b);      //EqT2.6
	        if(debug) COUT(a);
       		if(debug) COUT(qvec);
  }

FUNCTION mcOutput
  // Output parameters sampled; printed to 'cout' during the '-eval' phase.
  // Pipe to a file, such as vbmc.dat. Each line has 6 values.
    cout << log_M << " " << ah << " " << gh << " " << qtil << " "<< b << " "<< N << " "<< f << endl;
    
//______________________________
//Functions to link p and q -- overloaded
//______________________________

FUNCTION dvar_vector p2q(const dvar_vector pp,const dvariable aa, const dvariable bb)
	 {
		   return invlogit(aa - bb*logit(pp));
	 }

FUNCTION dvector p2q(const dvector pp,const double aa, const double bb)
	 {
		   return invlogit(aa - bb*logit(pp));
	 }

FUNCTION dvar_vector q2p(const dvar_vector qq,const double aa, const double bb)
	 {
		   return invlogit((aa - logit(qq))/bb);
	 }
	 
FUNCTION dvector q2p(const dvector qq,const double aa, const double bb)
	 {
		   return invlogit((aa - logit(qq))/bb);
	 }

//______________________________
//	Statistical Functions
//______________________________
//Logit -- overloaded
FUNCTION double logit(const double p)
 {
	    return log(p/(1.-p));
 }// end function

FUNCTION dvariable logit(const dvariable p)
 {
	    return log(p/(1.-p));
 }// end function

FUNCTION dvector logit(const dvector p)
  {
 	    return log(elem_div(p,(1.-p)));
  }// end function

FUNCTION dvar_vector logit(const dvar_vector p)
  {
 	    return log(elem_div(p,(1.-p)));
  }// end function

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Invlogit  -- overloaded
FUNCTION double invlogit(const double x)
  {
	  return mfexp(x)/(1. + mfexp(x));
  }// end function

FUNCTION dvariable invlogit(const dvariable x)
  {
	  return mfexp(x)/(1. + mfexp(x));
  }// end function

FUNCTION dvector invlogit(const dvector x)
   {
 	 return elem_div(mfexp(x),(1. + mfexp(x)));
   }// end function

FUNCTION dvar_vector invlogit(const dvar_vector x)
   {
 	 return elem_div(mfexp(x),(1. + mfexp(x)));
   }// end function


REPORT_SECTION
	report<<"#Objective function value"<<endl;
	REPORT(f);
	report<<"#Data"<<endl;
	REPORT(Fish);
	REPORT(yObs);
	REPORT(xi);
	REPORT(age);
	report<<"#Estimated Parameters"<<endl;
	REPORT(log_M);
	REPORT(ah);
	REPORT(gh);
	REPORT(qtil);
	REPORT(b);
	REPORT(a);
	REPORT(N);
	report<<"#Model Dimensions"<<endl;
	REPORT(n);
	report<<"#Other variables"<<endl;
	REPORT(M);
	REPORT(Si);
	REPORT(Betai);
	REPORT(rbar);
	REPORT(Ri);
	REPORT(pvec);
	REPORT(qvec);
	report << "$mcnames" << endl;
	report << "logM ah gh qtil b N fval" << endl;
 	report << "$mcest" << endl;
	report<< log_M << " " << ah << " " << gh << " " << qtil << " "<< b << " "<< N << " "<< f << endl;
	  	
GLOBALS_SECTION
	 #undef REPORT
	#define REPORT(object) report <<#object <<"\n" << object << endl;

	#undef COUT
	#define COUT(object) cout << #object "\n" << object <<endl;

	#if defined(_WIN32) && !defined(__linux__)
		const char* PLATFORM = "Windows";
	#else
		const char* PLATFORM = "Linux";
	#endif
	
	#include <admodel.h>
	#include <string.h>
	#include <statsLib.h>
	
TOP_OF_MAIN_SECTION
  arrmblsize = 1000000;  										
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000000); 	
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e5);  		
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2100000000); 			
  gradient_structure::set_MAX_NVAR_OFFSET(1000); 	

