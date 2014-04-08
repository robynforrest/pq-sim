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
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pq.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  A.allocate("A");
  rbar.allocate("rbar");
  yObs.allocate(1,n,"yObs");
  debug.allocate("debug");
  eof.allocate("eof");
		 cout<<"eof = "<<eof<<endl;
	   	  if(eof==999){
	   		cout<<"\n -- SUCCESSFULLY READ DATA FILE -- \n"<<endl;
	   	  }else{
	   		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	   	  }
  age.allocate(1,n);
  xi.allocate(1,n);
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
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_Z.allocate(-3.5,-0.2,1,"log_Z");
  beta1.allocate(0.0001,1,1,"beta1");
  alpha.allocate(0.0001,10.,1,"alpha");
  qtil.allocate(0.000001,1,1,"qtil");
  b.allocate(0.000001,4.,1,"b");
  N.allocate(0.000001,500.,-1,"N");
  a.allocate("a");
  #ifndef NO_AD_INITIALIZE
  a.initialize();
  #endif
  Z.allocate("Z");
  #ifndef NO_AD_INITIALIZE
  Z.initialize();
  #endif
  Si.allocate(1,n,"Si");
  #ifndef NO_AD_INITIALIZE
    Si.initialize();
  #endif
  Betai.allocate(1,n,"Betai");
  #ifndef NO_AD_INITIALIZE
    Betai.initialize();
  #endif
  Ri.allocate(1,n,"Ri");
  #ifndef NO_AD_INITIALIZE
    Ri.initialize();
  #endif
  pvec.allocate(1,n,"pvec");
  #ifndef NO_AD_INITIALIZE
    pvec.initialize();
  #endif
  qvec.allocate(1,n,"qvec");
  #ifndef NO_AD_INITIALIZE
    qvec.initialize();
  #endif
  pPrime.allocate(1,n,"pPrime");
  #ifndef NO_AD_INITIALIZE
    pPrime.initialize();
  #endif
  log_resid.allocate(1,n,"log_resid");
  #ifndef NO_AD_INITIALIZE
    log_resid.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
	getQi(); 	//Calls getPi 
	getpPrime();
	calcObjectiveFunction();
}

void model_parameters::calcObjectiveFunction(void)
{
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
}

void model_parameters::getpPrime()
{
  {
	  pPrime=elem_prod(xi,pvec);
	  pPrime/=sum(pPrime); //EqT4.2
  }
}

void model_parameters::getPi()
{
  {
		//Function to estimate proportions at age and sampled proportions at age
		//This function is called by getQi()
		//Declare and initialize variables
		int j;
		dvariable surv; 
		dvariable aminus;
		Si.initialize(); Betai.initialize(); Ri.initialize(); pvec.initialize(); qvec.initialize();
		//Proportions at age -- RF thinks the F component of Z needs to be modified by the selectivity
		//Survivorship at age	 Eq. T6.1
		Z=mfexp(log_Z);
		surv=mfexp(-Z); 
		Si(1)=1.;
		for(j=1+1; j<=n; j++) Si(j)=surv*Si(j-1);
		Si(n) /=(1. - surv);
			if(debug) {
				COUT(Z);
			        COUT(surv);
				COUT(Si);
			}
		//Selectivity	Eq T6.2
		aminus=A-1;
		Betai=1.; //selectivity - fill with 1: ages < A overwritten in next line 
		for(j=1; j<=aminus; j++) Betai(j) = 1. - (1.-beta1) * pow(((A-j)/aminus),alpha);
			     if(debug) {
			     	COUT(Betai);
			   	COUT(beta1);
				COUT(alpha);
			   }
		//Recruits to each age class -- placeholder -- add time and anomalies later 
		Ri = rbar;	  //Eq T6.3
		//Get pvec -- proportions at age in the population
		pvec = elem_prod(elem_prod(Si,Betai),Ri); 
		pvec /= sum(pvec);   //EqT6.4	 -- normalise
			if(debug) COUT(pvec);
 }
}

void model_parameters::getQi()
{
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
}

dvar_vector model_parameters::p2q(const dvar_vector pp,const dvariable aa, const dvariable bb)
{
	 {
		   return invlogit(aa - bb*logit(pp));
	 }
}

dvector model_parameters::p2q(const dvector pp,const double aa, const double bb)
{
	 {
		   return invlogit(aa - bb*logit(pp));
	 }
}

dvar_vector model_parameters::q2p(const dvar_vector qq,const double aa, const double bb)
{
	 {
		   return invlogit((aa - logit(qq))/bb);
	 }
}

dvector model_parameters::q2p(const dvector qq,const double aa, const double bb)
{
	 {
		   return invlogit((aa - logit(qq))/bb);
	 }
}

double model_parameters::logit(const double p)
{
 {
	    return log(p/(1.-p));
 }// end function
}

dvariable model_parameters::logit(const dvariable p)
{
 {
	    return log(p/(1.-p));
 }// end function
}

dvector model_parameters::logit(const dvector p)
{
  {
 	    return log(elem_div(p,(1.-p)));
  }// end function
}

dvar_vector model_parameters::logit(const dvar_vector p)
{
  {
 	    return log(elem_div(p,(1.-p)));
  }// end function
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //Invlogit  -- overloaded
}

double model_parameters::invlogit(const double x)
{
  {
	  return mfexp(x)/(1. + mfexp(x));
  }// end function
}

dvariable model_parameters::invlogit(const dvariable x)
{
  {
	  return mfexp(x)/(1. + mfexp(x));
  }// end function
}

dvector model_parameters::invlogit(const dvector x)
{
   {
 	 return elem_div(mfexp(x),(1. + mfexp(x)));
   }// end function
}

dvar_vector model_parameters::invlogit(const dvar_vector x)
{
   {
 	 return elem_div(mfexp(x),(1. + mfexp(x)));
   }// end function
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	report<<"#Objective function value"<<endl;
	REPORT(f);
	report<<"#Data"<<endl;
	REPORT(yObs);
	REPORT(xi);
	REPORT(age);
	report<<"#Estimated Parameters"<<endl;
	REPORT(log_Z);
	REPORT(beta1);
	REPORT(alpha);
	REPORT(qtil);
	REPORT(b);
	REPORT(a);
	report<<"#Model Dimensions"<<endl;
	REPORT(n);
	REPORT(A);
	REPORT(N);
	report<<"#Other variables"<<endl;
	REPORT(Z);
	REPORT(Si);
	REPORT(Betai);
	REPORT(rbar);
	REPORT(Ri);
	REPORT(pvec);
	REPORT(qvec);
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 1000000;  										
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000000); 	
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e5);  		
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2100000000); 			
  gradient_structure::set_MAX_NVAR_OFFSET(1000); 	
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
