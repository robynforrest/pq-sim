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
	init_bounded_number Phi1(-2.99,1.,1);          	//log_M Natural mortality	      
	init_bounded_number Phi2(1,n,1);      			//ah Age at 50% selectivity
	init_bounded_number Phi3(0.0001,3.,1);             	//gh SD in selectivity
	init_bounded_number Phi4(0.0001,0.999,1);	 	//qtil Parameter relating p to q
	init_bounded_number Phi5(0.0001,3.,1);			//b Parameter relating p to q	 
	init_bounded_number Phi6(5.,500.,1);	         		//N Effective sample size for Dirichlet distribution -- doesn't work very well if estimating this
	
	//Estimated parameters -- pick them up in initModel() -- some ADMB functions aren't overloaded for init_numbers
	number log_M;  
	number ah;
	number gh;
	number qtil;
	number b;
	number N;
	
	//Derived parameters
	number a;
	number M;
 	vector Si(1,n);
        vector Betai(1,n);
 	vector Ri(1,n);		
 	vector pvec(1,n);		
	vector qvec(1,n);
	vector pPrime(1,n);
	//vector log_resid(1,n);

	objective_function_value f;
	sdreport_number sdreport;

PROCEDURE_SECTION
	initModel();
	getQi(); 	//Calls getPi 
	getpPrime();
	calcObjectiveFunction();
	sdreport=f;
	if(mceval_phase()) mcOutput();
		
FUNCTION initModel
	 log_M = Phi1;              //Natural mortality	      
	 ah = Phi2;      //Age at 50% selectivity
	 gh = Phi3;             //SD in selectivity
	 qtil = Phi4;	   //Parameter relating p to q
	 b = Phi5;		 //Parameter relating p to q	 
	 N = Phi6;	 

FUNCTION calcObjectiveFunction
 {
       //Calculate negative log likelihood from T4.4 
       //PSEUDOCODE FOR LIKELIHOOD
       //like = sum(gammaln(N*pPrime) - (N*pPrime - 1)*log(yObs)) - sum((1 - xi)*log(qvec) + xi*log(1-qvec)) - gammaln(N);
       
       dvariable like1;	 //First part of objective function - sum over pPrime >0
       dvariable like2;	//Second part of objective function - sum over pPrime >0
      dvariable like;
       
       //Build components of objective function
       dvar_vector Npp = N*pPrime;			//N*pPrime
       dvar_vector Nppminus = N*pPrime - 1.;	//N*pPrime - 1
       //dvar_vector minusxi = -xi;
      // dvar_vector minusqi = -qvec;
       //minusxi +=1.; 							// 1-xi
      // minusqi +=1.; 							// 1-qvec 
      
      //Calculate objective function
      like1=0.;
      like2=0.;
      //Calculate like1 and like2, the first parts of objective function :: only for pPrime>0
      for(int i=1; i<=n;i++){
	   if(pPrime(i)>0) {
	   	like1 += gammln(Npp(i)) - (Nppminus(i))*log(yObs(i)); //
		like2 += (1.-xi(i))*log(qvec(i)) + xi(i)*log(1.-qvec(i));
		} //end if
	} //end for i
      //Calculate the objective function
     like = like1 - like2 - gammln(N);
    
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
		dvar_vector Z(1,n); 
		dvar_vector surv(1,n); 
				
		dvariable aminus;
		Si.initialize(); Betai.initialize(); Ri.initialize(); pvec.initialize(); qvec.initialize(); surv.initialize(); Z.initialize(); 
               	
               	//1. Logistic selectivity
               	Betai = plogis<dvar_vector>(age,ah,gh);
       		
       		//2. Get fished survival rate at age 
		M=mfexp(log_M);
		
		//Get Zi = M + Betai*Fi , where fishing mortality at age is modified by selectivity at age (0 < Betai <=1)
		Z=M;
		Z+=Betai*Fish;
		surv=mfexp(-Z); 
				
		//3. Get survivorship -- a function of natural mortality, fishing mortality and selectivity
		Si(1)=1.;
		for(j=2; j<=n; j++)  Si(j)=Si(j-1)*surv(j-1);  
		Si(n) /=(1. - surv(n));
		
		if(debug) {
				COUT(Fish);
				COUT(M);
				COUT(Z)
			        COUT(surv);
				COUT(Si);
				COUT(Betai);
				COUT(ah);
				COUT(gh);
			}
				
		//Recruits to each age class -- placeholder -- add time and anomalies later 
		Ri = rbar;	  //Eq T6.3
		
		//Get pvec -- vulnerable proportions at age in the population
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
 {
	  // Output parameters sampled; printed to 'cout' during the '-eval' phase.
	  // Pipe to a file, such as vbmc.dat. Each line has 6 values.
	   if(nf==1){
			cout<<"\n -- RUNNING MCEVAL -- \n"<<endl;
			ofstream ofs("pqpars.csv");
			ofs<< "log_M, " << "ah, " << "gh, " << "qtil,"<<" b, "<< " N,"<< "f"<< endl;
			ofstream ofprop("pqproportions.csv");
			ofstream ofsel("pqselectivity.csv");
	   }

	  ofstream ofs("pqpars.csv",ios::app); 
	  ofs << log_M << ", " << ah << ", " << gh << ", " << qtil << " ,"<< b << ", "<< N << " ,"<< f << endl;
	  
	 //Do this with a loop -- ADMB putting an NA in first column otherwise
	 ofstream ofsel("pqselectivity.csv",ios::app);
	  for(int jj=1; jj<=n;jj++) {
	  	if(jj<n) ofsel<<Betai(jj)<<","; 
	  	if(jj==n) ofsel<<Betai(jj)<<endl;
          }
          ofstream ofprop("pqproportions.csv",ios::app);
	  for(int jj=1; jj<=n;jj++) {
		if(jj<n) ofprop<<pvec(jj)<<","; 
		if(jj==n) ofprop<<pvec(jj)<<endl;
          }
  }
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
	report<<"$f"<<"\n"<<f<<endl;
	report<<"#Data"<<endl;
	report << "$Fish" << "\n" << Fish <<endl;
	report << "$yObs" << "\n" << yObs<<endl;
	report << "$xi" << "\n" << xi<<endl;
	report << "$age" << "\n" << age<<endl;
	report<<"#Estimated Parameters"<<endl;
	report << "$log_M" << "\n" <<log_M <<endl;
	report << "$ah" << "\n" <<ah <<endl;
	report << "$gh" << "\n" <<gh <<endl;
	report << "$qtil" << "\n" <<qtil <<endl;
	report << "$b" << "\n" <<b <<endl;
	report << "$a" << "\n" <<a <<endl;
	report << "$N" << "\n" <<N <<endl;
	report<<"#Model Dimensions"<<endl;
	report << "$n" << "\n" <<n<<endl;
	report<<"#Other variables"<<endl;
	report << "$M" << "\n" <<M <<endl;
	report << "$Si" << "\n" <<Si <<endl;
	report << "$Betai" << "\n" <<Betai <<endl;
	report << "$rbar" << "\n" <<rbar <<endl;
	report << "$Ri" << "\n" <<Ri <<endl;
	report << "$pvec" << "\n" <<pvec <<endl;
	report << "$qvec" << "\n" <<qvec <<endl;
	report << "$mcnames" << endl;
	report << "logM ah gh qtil b N fval" << endl;
 	report << "$mcest" << endl;
	report<< log_M << " " << ah << " " << gh << " " << qtil << " "<< b << " "<< N << " "<< f << endl;
	  	
GLOBALS_SECTION
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
  gradient_structure::set_MAX_NVAR_OFFSET(30000000); 	

