#include <cmath>
#include <iostream>
#include <string>

#include "params.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

//it would not be efficient to select "distrib for beam x" every time: must find a way to select it and tell it to sub functions

//

//int gsl_integration_qagi (gsl_function * f, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)

//params of nucleon distribution
struct nuc_params {

	double rho0;
	double Rhalf;
	double a;
};

//params for the geometry of the collision
struct geom_params {
	
	//info already contained in gen but useful
	double A1,p1,n1;
	double A2,p2,n2;

	//abserr != 0 is irrelevant here (I think)
	double abserr;
	double relerr;
	int maxiter;

	//
	std::string nucdatafile;
	
	//'h' for Hulthen, 'w' for Woods-Saxon, ...
	char distrib_1,distrib_2;
	
	//Params if Woods-Saxon
	nuc_params fpar_1, fpar_2;
	
	//experimental input
	double sigma_NN;

	//centrality class
	double lcc;
	double ucc;

};

//complete params
struct full_params {
	geom_params * geom;
	params * gen;
};

//Compute \int_bx^by{b*TAB(b)}/\int_bx^by{b*(1-exp(-sigma_NN*TAB(b)))} (2*pi factors cancel)
//which is the average over a x-y centrality bin of the nuclear overlap function from Glauber's model
double geom_factor(params * genp) {

	// Declaration of functions //
	double find_bx(double, double, full_params *);
	double glauber_inel_cs(double, double, full_params *, bool);

	// Declaration of variables //
	//params local_params;
	//fermi_params tmp_fpars;

	double num, denum;
	double res;
	double bx, by;
	
	full_params fullpars;
	geom_params geompars;

	// Initialization of variables //
	//local_params = * (params *) pars;
	
	// Core //

	//initialize geom params
	
	geompars.abserr=0.1;
	geompars.relerr=0.0;
	geompars.maxiter=1000;
	
	geompars.p1=(genp->beam_1)[0];
	geompars.n1=(genp->beam_1)[1];
	geompars.A1=geompars.p1+geompars.n1;
	
	geompars.p2=(genp->beam_2)[0];
	geompars.n2=(genp->beam_2)[1];
	geompars.A2=geompars.p2+geompars.n2;

	geompars.sigma_NN=7.7; //fm^2

	geompars.lcc=0.0;
	geompars.ucc=1.0;

	//
	fullpars.geom=&geompars;
	fullpars.gen=genp;

	//choose distrib for beam 1
	//choose_distrib(1,&local_params,&gpars);
	//equivalent to FindNucleusData?
	
	//choose distrib for beam 2
	//choose_distrib(2,&local_params,&gpars);
	//equivalent to FindNucleusData?
	
	//tabulate TAB

	//find bx and by of centrality class to an accuracy dx (second parameter)
	bx=find_bx(geompars.lcc, .01, &fullpars);
	by=find_bx(geompars.ucc, .01, &fullpars);

	std::cout << "bx=" << bx << "& by=" << by << "\n";

	//
	num=glauber_inel_cs(bx, by, &fullpars, 1);
	denum=glauber_inel_cs(bx, by, &fullpars, 0);
	
	//integrate geom_factor_sub_x over x from -inf to inf
	
	//for(int i=1;i < 10;i++) {
	//	std::cout << "bx at " << i*10 << "% centrality=" << find_bx(.1*i, .01, &fullpars) << "\n";
		//std::cout << "bx at " << .99999999 << "% centrality=" << find_bx(0.0, .01, &fullpars) << "\n";
	//}

	//std::cout << "208=?" << 2.0*M_PI*num << "\n";

	std::cout << "res=" << num/denum << "\n";

	res=1.0;
	
	return res;

}





void choose_distrib(int beam_id, void * pars, void * gpars) {

	// Declaration of functions //
	/*fermi_params set_fermi_params(int []);*/

	// Declaration of variables //
	/*char deuteron_dist = 'h';
	params local_params;
	geom_params export_gpars;*/
	
	/*int beam[2];
	char * distrib;
	nuc_params * fpars;*/
		
	// Initialization of variables //
	/*local_params = * (params *) pars;
	export_gpars = * (geom_params *) gpars;*/
	
	// Core //
	
	/*
	if (1 == beam_id) {
		beam[0]=local_params.beam_1[0];
		beam[1]=local_params.beam_1[1];
		distrib=&export_gpars.distrib_1;
		fpars=&export_gpars.fpar_1;
	}
	else {
		if (2 == beam_id) {
			beam[0]=local_params.beam_2[0];
			beam[1]=local_params.beam_2[1];
			distrib=&export_gpars.distrib_2;
			fpars=&export_gpars.fpar_2;
		}
	}
	
	
	//if deuteron
	if ((1 == beam[0])&&(1 == beam[1])) {
		
		//if W-S distrib
		if ('w' == deuteron_dist) {
			
			//set params
			*fpars=set_fermi_params(beam);
	
			//tell it's W-S
			*distrib='w';
			
		}
		else {
			//if Hulthen
			if ('h' == deuteron_dist) {
				
				//tell it's Hulthen
				*distrib='h';
				
			}
		}
	}
	//if anything else
	//be careful! Woods-Saxon distribution may not be appropriate for every nuclei.
	//also, only a handful of parameters for specific nuclei are specified in set_fermi_params(beam).
	//while a general formula is implemented in that function, the values it returns are known to be
	//inadequate for small (A<20) nuclei. anyway, do not rely on this formula if you can avoid it.
	else {
		
		//set params
		*fpars=set_fermi_params(beam);
	
		//tell it's W-S
		*distrib='w';
		
	}
	
	*/
}

/*
//Derived from Hulthen wave function
//Params from
double deuteron_density(double r) {
		
	// Declaration of functions //

	// Declaration of variables //
	double res;
	double a,b,rho0;

	// Initialization of variables //
	a=0.457;
	b=2.35;
	rho0=0.0788;
	
	// Core //

	res=rho0*pow((exp(-a*r)+exp(-b*r))/r,2);

	return res;
		
}
*/
/*
//Fermi distribution for nuclear density
double nucleon_density(double r, fermi_params pars) {

	// Declaration of functions //

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	res=(pars.rho0)/(1.0+exp((r-(pars.Rhalf))/(pars.a)));

	return res;

}
*/
/*
//a and Rhalf are in fm
//rho0 in fm^-3

fermi_params set_fermi_params(int beam[2]) {
	
	//, fermi_params * fpars

	// Declaration of functions //

	// Declaration of variables //
	double rho0, Rhalf, a;
	int n, p, tot;
	fermi_params fpars;

	// Initialization of variables //
	p=beam[0];
	n=beam[1];
	tot=n+p;

	// Core //

	switch (tot) {
	
		//Deuteron (and diproton and dineutron?)
		//a and Rhalf from projects.hepforge.org/tglaubermc/TGlauberMC_v1.1.pdf
		//rho0 from normalization
		case 2:
			a=0.5882;
			Rhalf=0.01;
			rho0=0.427;
			break;
		
		//Gold (from De Vries et al, 1987)
		case 197:
			a=0.535;
			Rhalf=6.38;
			rho0=0.17;
			break;	
		
		//Bertulani, p.100
		//According to Bertulani, works fairly well if nb>20
		default:
			if (tot > 20) {
				a=0.54;
				rho0=0.17;
				Rhalf=1.128*pow(tot,1.0/3.0)-0.89*pow(tot,-1.0/3.0);
			}
			else {
				std::cout << "No known parameters of beam in set_fermi_params(). Aborting...";
				exit(1);
			}
			break;
		
	}
	
	fpars.a=a;
	fpars.Rhalf=Rhalf;
	fpars.rho0=rho0;
	
	return fpars;

}
*/

//Another ^*(&$ struct b/c need to pass void pointer
//There _must_ be something I'm not doing right b/c it seems very stranger that I have to define a new struct each time I use GSL
//f(x)=a
struct solver_params {
	full_params * params;
	double a;
};


//Find bx such that \int_0^bx {db b (1-exp(-sigma_NN*TAB(b)))} = (X+/-dX) * \int_0^inf{db b (1-exp(-sigma_NN*TAB(b)))}
double find_bx(double x, double dx, full_params * fullp) {
	
	// Declaration of functions //
	double glauber_inel_cs(double, double, full_params *, bool);
	double find_bx_sub(double bx, void * params);
	
	// Declaration of variables //
	double res;
	double tot_cs;
	double bmax;
	int iter, maxiter, status, check_iter;
	
	solver_params solverp;
	
	gsl_function rf_funct;
	
	// Initialization of variables //
	bmax=2.0; //fm
	
	iter=0;
	maxiter=100;
	
	solverp.params=fullp;

	// Core //

	//preliminary verifications
	if ((x > 1)||(x < 0)) {
		std::cout << "x outside range ([0,1]) in find_bx(). Aborting...\n";
		exit(0);
	}
	else {

	//compute tot_cs=\int_0^inf{db b (1-exp(-sigma_NN*TAB(b)))} (actually, total cross-section modulo 2*PI)
	tot_cs=glauber_inel_cs(0, GSL_POSINF, fullp, 0);
	
	//std::cout << "totcs=" << tot_cs << "\n";
	
	//if x = 1, res=infinity (theoretically), although a small finite value should do due to the fast decreases of the nucleon distrib
	if (1.0 == x) {
		res=GSL_POSINF;
	}
	else {
	
		solverp.a=x*tot_cs;
		
		//find a b such that f(x) > x TAB(b)
		//required (AFAIK) for the type of root solving method used below (bracketing)
		do {
			bmax*=bmax;
		}
		while (find_bx_sub(bmax,&solverp) < 0);
		
		//std::cout << "bmax=" << bmax << "\n";
		
		//use GSL to find the root of abs(\int_0^bx {db b (1-exp(-sigma_NN*TAB(b)))} - X Z) such that abs(...) <= dx*tot_cs
			
		//
		rf_funct.function=&find_bx_sub;
		rf_funct.params=&solverp;
		
		//type of root solver
		//const gsl_root_fsolver_type * rf_solver_type = gsl_root_fsolver_bisection;
		const gsl_root_fsolver_type * rf_solver_type = gsl_root_fsolver_falsepos;
		//const gsl_root_fsolver_type * rf_solver_type = gsl_root_fsolver_brent;
		
		//This function returns a pointer to a newly allocated instance of a solver of type T
		gsl_root_fsolver * rf_solver = gsl_root_fsolver_alloc(rf_solver_type);
		
		//This function initializes, or reinitializes, an existing solver s to use the function f and the initial search interval [x_lower, x_upper]
		gsl_root_fsolver_set(rf_solver, &rf_funct, 0, bmax);
		
		do {
			//These functions perform a single iteration of the solver s.
			check_iter=gsl_root_fsolver_iterate (rf_solver);
		
			if (check_iter == GSL_EBADFUNC) {
				std::cout << "The iteration encountered a singular point where the function or its derivative evaluated to Inf or NaN in find_bx(). Aborting...\n";
				exit(0);
			}
		
			//These functions return the current estimate of the root for the solver s. 
			res=gsl_root_fsolver_root (rf_solver);
		
			//This function tests the residual value f against the absolute error bound epsabs. The test returns GSL_SUCCESS if the following condition is achieved, |f| < epsabs
			status=gsl_root_test_residual (find_bx_sub(res,&solverp), dx*tot_cs);
		
			//std::cout << iter << ": res=" << res << "\n";
		
			//
			if (iter >= maxiter) {
				std::cout << "Can't find bx within the specified maximum number of iterations in find_bx(). Aborting...\n";
				exit(0);
			}
			
			iter++;
		
		} while (status != GSL_SUCCESS);
		
		//These functions free all the memory associated with the solver s
		gsl_root_fsolver_free(rf_solver);
	
	}

	}

	return res;

}

//
double find_bx_sub(double by, void * params) {
	
	// Declaration of functions //
	double glauber_inel_cs(double, double, full_params *, bool);

	// Declaration of variables //
	double res;
	solver_params solverp;

	// Initialization of variables //
	solverp = * (solver_params *) params;

	// Core //
	res=glauber_inel_cs(0, by, (solverp.params), 0)-(solverp.a);
	
	return res;
	
}

//Compute the inelastic cross-section (modulo 2*PI) using Glauber's model: \int_0^bx {db b (1-exp(-sigma_NN*TAB(b)))}
//if small_sigma_approx = 1, integrate TAB(b) instead of (1-exp(-sigma_NN*TAB(b)))
//sigma_NN*TAB(b) is actually a first order approximation of (1-exp(-sigma_NN*TAB(b))) if sigma_NN is small (as in hard NN collisions)
//use qag for the integration for bx < infinity, and qagi for bx=infinity
double glauber_inel_cs(double bx, double by, full_params * fullp, bool small_sigma_approx) {

	// Declaration of functions //
	double glauber_inel_cs_sub(double, void *);
	double glauber_inel_cs_sub_approx(double, void *);

	// Declaration of variables //
	double res;
	size_t maxiter;
	gsl_integration_workspace * qag_ws;
	gsl_function qagf;
	
	double abserr, relerr, err;

	// Initialization of variables //
	maxiter=fullp->geom->maxiter;
	abserr=fullp->geom->abserr;
	relerr=fullp->geom->relerr;

	// Core //

	//function to integrate
	if (small_sigma_approx) {
		qagf.function = &glauber_inel_cs_sub_approx;
	}
	else {
		qagf.function = &glauber_inel_cs_sub;
	}
	//params to pass to the function
	qagf.params = fullp;
		
	//This function allocates a workspace sufficient to hold n double precision intervals, their integration results and error estimates.
	//gsl_integration_workspace * gsl_integration_workspace_alloc (size_t n);
 	qag_ws=gsl_integration_workspace_alloc(maxiter);
			
	
	//compute \int_0^bx{db b (1-exp(-sigma_NN*TAB(b)))} 
	if ((GSL_POSINF == by)&&(bx>=0.0)) {
		//...using QAGI
		//This function computes the integral of the function f over the semi-infinite interval (a,+\infty).
		//int gsl_integration_qagiu (gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)
		gsl_integration_qagiu(&qagf, bx, abserr, relerr, maxiter, qag_ws, &res, &err);

		std::cout << "sigma_geo?=" << res << "\n";

	}
	else {

		if ((by >= bx)&&(bx >= 0)) {
			//...using QAG
			//This function applies an integration rule adaptively until an estimate of the integral of f over (a,b) is achieved within the desired absolute and relative error limits, epsabs and epsrel.
			//int gsl_integration_qag (const gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr)
			gsl_integration_qag(&qagf, bx, by, abserr, relerr, maxiter, 1, qag_ws, &res, &err);
		}
		else {
			std::cout << "by <= bx or bx negative in glauber_inel_cs(). Aborting...\n";
			exit(0);
		}
	}
	
	//This function frees the memory associated with the workspace w.
	//void gsl_integration_workspace_free (gsl_integration_workspace * w)
	gsl_integration_workspace_free(qag_ws);

	//need to multiply the result by sigma_NN if small_sigma_approx is used
	return res;

}

//integrand of glauber_inel_cs() (if not using small_sigma_approx): b*(1-exp(-sigma_NN*TAB(b)))
//sigma_NN is the inelastic NN (usually pp) cross-section
//glauber_inel_cs_sub() will be in fm^2 (if all the parameter have the right units)
double glauber_inel_cs_sub(double b, void * params) {

	// Declaration of functions //
	double TAB_Wong(double, void *);

	// Declaration of variables //
	double res;
	double sigma_NN;
	full_params * fullp;

	// Initialization of variables //
	fullp = (full_params *) params;

	//inelastic NN (usually pp) cross-section
	sigma_NN=fullp->geom->sigma_NN; //fm^2
	
	//std::cout << "A=" << A << "&B=" << B <<"\n";
	
	// Core //

	//b (1-exp(-sigma_NN*TAB(b)))
	res=b*(1.0-exp(-1.0*sigma_NN*TAB_Wong(b,fullp)));

	return res;

}

//integrand of glauber_inel_cs() (if using small_sigma_approx): b*TAB(b)
//glauber_inel_cs_sub() will be in fm^2 (if all the parameter have the right units)
double glauber_inel_cs_sub_approx(double b, void * params) {

	// Declaration of functions //
	double TAB_Wong(double, void *);

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	//b (1-exp(-sigma_NN*TAB(b)))
	res=b*TAB_Wong(b,params);

	return res;

}



//From Wong - Introduction to High-Energy Heavy_ion Collisions, chap 13
//Careful about the normalization, Wong's density function are normalized to 1, not the number of nucleon
//s in fm
//Return TAB in fm^-2
double TAB_Wong(double s, void * fullp) {

	// Declaration of functions //

	// Declaration of variables //
	double res;
	double beta_sq, p1,p2,p3, r0;
	double A1,A2;
	
	full_params fullpars;

	// Initialization of variables //
	fullpars = * (full_params *) fullp;

	A1=fullpars.geom->A1;
	A2=fullpars.geom->A2;

	// Core //
	
	r0=1.05; //fm
	p1=r0*pow(A1,1.0/3.0)/sqrt(3.0);
	p2=r0*pow(A2,1.0/3.0)/sqrt(3.0);
	p3=0.68; //fm	
	
	beta_sq=p1*p1+p2*p2+p3*p3;
	
	//integrate geom_factor_sub_y over y from -inf to inf
	res=A1*A2*exp(-s*s/(2.0*beta_sq))/(2.0*M_PI*beta_sq);

	return res;
	
}

//TAB tabulated
double TAB(double s, void * fullp) {

	// Declaration of functions //

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	res=1.0;

	return res;
}

//
void tabulate_TAB() {
	
	// Declaration of functions //

	// Declaration of variables //

	// Initialization of variables //

	// Core //

	
	

}