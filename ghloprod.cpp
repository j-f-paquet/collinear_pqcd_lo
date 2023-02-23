#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_integration.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//#include "LHAPDF/LHAPDF.h"

#include "params.h"

using namespace std;

int main(int argc, char *argv[]) {

	// Declaration of functions //
	double partprod_pp(void *);
	double kfrag(double);
	double kdirect(double);
	double unit_fct(double);
	void save_info(ofstream *, params *);
	//double geom_factor(params *);

	// Declaration of variables //
	params pars;
  	ofstream outfile;
	char filename[200], origfn[200];
	int id;
	bool selname;
	double nnprod,res;

	// Initialization of variables //
	strcpy(filename,"./example.dat");
	strcpy(origfn,filename);
	id=1;

	// Core //

			
	//autoselect filename
	selname=0;
	do {
		outfile.open (filename,ios::in);
		if (outfile.is_open()) {
			outfile.close();
			sprintf(filename,"%s%d",origfn,id);
			id++;
			selname=1;
		}
		else {
			selname=0;
		}
	
	} while (selname);
	
	outfile.close();
	
	//
	outfile.open (filename,ios::out);
		
	//beams
	//nb of protons
	pars.beam_1[0]=1;
	//nb of neutrons
	pars.beam_1[1]=0;

	//beam 2
	//gold is 79 protons, 118 neutrons
	//Pb is 82 protons,125 neutrons
	pars.beam_2[0]=1;
	pars.beam_2[1]=0;
	
	//general params
	
		//number of flavours
		pars.nf=3;

		//lambda QCD
		pars.lambda=.2;

		//
		pars.sqrts=200.0;

		//
		pars.rapidity=0.0;
		
		//unit conversion factor
		pars.units_conv=0.389379;
		
		
	//what do you want?
	//for direct photon, part_type=10 and source=0
	//for frag photon, part_type=10 and source=1
	//for frag pions, part_type=5 and source=1
	//I suspect that other combinations give meaningless results
		
	//direct prod (source=0)
	
		//particle type
		pars.part_type=10;

		//direct or fragmentation?
		pars.source=0;
		
		//factorisation and fragmentation scales
		pars.Qfac='j';
		pars.Qren='j';
		pars.Qfac_norm=1.0/sqrt(2.0);
		pars.Qren_norm=1.0/sqrt(2.0);

		//K factor
		pars.Kfunct=unit_fct;
		
		//
		pars.maxiter=5000;
	
		//
		pars.qag_method=1;
	
		//
		pars.abserr=0.0;
		pars.relerr=1e-3;
		
		//
		pars.cteq_table=1;
		strcpy(pars.pdf_name,"cteq5l.LHgrid");
		
		//save params to file
		save_info(&outfile,&pars);
		outfile << "#\n#pT\tcs\n";

		//
		//LHAPDF::initPDFByName(pars.pdf_name,1);
		
		//
		for(int i=10;i<=40;i++) {

			pars.pJt=i*.5;

			nnprod=partprod_pp(&pars);

			//I don't think geom_factor has to be recalculated every times
			//res=nnprod*geom_factor(&pars);
			res=nnprod; //*geom_factor(&pars);

			//cout << partprod_pp(&pars) << "\n";
		
			outfile << pars.pJt << "\t" << res << "\n";

		}
	
	//frag prod (source=1)
	
	/*
		//particle type
		pars.part_type=10;

		//direct or fragmentation?
		pars.source=1;

		//factorisation and fragmentation scales
		pars.Qfac='j';
		pars.Qren='j';
		pars.Qfac_norm=1.0/sqrt(2.0);
		pars.Qren_norm=1.0/sqrt(2.0);

		//K factor
		pars.Kfunct=kfrag;
		
		//
		pars.frag_int_method='q';
		pars.maxiter=5000;
	
		//
		pars.qag_method=2;
	
		//
		pars.abserr=0.0;
		pars.relerr=1e-3;
		
		//
		//pars.cteq_table=3;
		strcpy(pars.pdf_name,"cteq5l.LHgrid");
		pars.kkp_order=0;
		pars.photon_set=1;
		
		//
		save_info(&outfile,&pars);
		outfile << "#\n#pT\tcs\n";
		
		LHAPDF::initPDFByName(pars.pdf_name,1);
		
		//
		for(int i=2;i<=2;i++) {
			
			pars.pHt=i*5.0;

			outfile << pars.pHt << "\t" << partprod_pp(&pars) << "\n";

		}	
	
	*/
	outfile.close();
		
}

//E_h d^3 sigma/d p_h^3
/*
source	direct	frag
id	0	1

part_type :
id	0	1		2		3		4		5	6		7
frag	-	(pi^+ + pi^-)/2	(K^+ + K^-)/2	(K^0+K^0_bar)/2	(p+p_bar)/2	(pi^0)	(n+n_bar)/2	(h^+ + h^-)
direct	gluon	u		d		s		c		b	t

id	10	-1	-2	-3	-4	-5	-6
frag	photon	-	-	-	-	-	-
direct	photon	ub	db	sb	cb	bb	tb

*/
double partprod_pp(void * pars) {

	// Declaration of functions //
	int two_two_type(int, int, int, int);
	double direct_prod(double, void *);
	double frag_prod_monte(double *, size_t, void *);
	double frag_prod_qag (double, void *);
	double cteq5_pdf(int, int, double, double);

	// Declaration of variables //
	params curr_set;
	int type;
	size_t maxiter;
	double lowbound, upbound;
	double abserr, relerr;
	int method;
	double result, error, total;
	double average_error;
	int nb;
	
	// Initialization of variables //
	 curr_set = * (params *) pars;
	
	 //
	 total=0.0;
	 nb=0;
	 average_error=0.0;
	 
	// Core //

	//


	 //direct
	 if (0 == curr_set.source) {

		 //
		 maxiter=curr_set.maxiter;

		int c;
		 
		gsl_function F;

	 	//initialisation for the numerical integration
		 gsl_integration_workspace * w = gsl_integration_workspace_alloc(maxiter);
		 
		 F.function = &direct_prod;
		 F.params = &curr_set;
		 
		lowbound=(2.0*curr_set.pJt/curr_set.sqrts)*exp(curr_set.rapidity)/(2.0-(2.0*curr_set.pJt/curr_set.sqrts)*exp(-1*curr_set.rapidity));
		 upbound=1.0;
		 
		 abserr=curr_set.abserr;
		 relerr=curr_set.relerr;
		 
		 method=curr_set.qag_method;

		for(int a = -1*curr_set.nf;a <= curr_set.nf;a++) {
		for(int b = -1*curr_set.nf;b <= curr_set.nf;b++) {
		for(int d = -1*curr_set.nf;d <= curr_set.nf;d++) {
	
			//c = part_type
			c = curr_set.part_type;
	
			//Determine the type
			type=two_two_type(a, b, c, d);

			//
			curr_set.parton_a=a;
			curr_set.parton_b=b;
			curr_set.parton_c=c;
			curr_set.parton_d=d;
			curr_set.type=type;

			result=0.0;

			//If type != 0
			if (type != 0) {
				
				//Integrate direct_prod over xa
				//cout << a << " + " << b << " + " << " -> " << c << " + " << d << " : " << type << "\n";

				gsl_integration_qag (&F, lowbound, upbound, abserr, relerr, maxiter, method, w, &result, &error); 

				nb++;
				average_error=((nb-1)*average_error+error/result)/double(nb);

			}
			
			total+=result;
	
		}
	
		}
	
		}

		//This function frees the memory associated with the workspace w. 
		gsl_integration_workspace_free(w);
		
		//
		total*=curr_set.units_conv*(*curr_set.Kfunct)(curr_set.pJt)*2.0/M_PI;
		
		//cout << "res=" << total << " and relerr~" << average_error << "\n";
		 
	}
	 
	 //frag
	if (1 == curr_set.source) {
		
		if ('m' == curr_set.frag_int_method) {
		
			//
			relerr=curr_set.relerr;
			
			double xl[2] = { (2.0*curr_set.pHt/curr_set.sqrts)*exp(curr_set.rapidity)/(2.0-(2.0*curr_set.pHt/curr_set.sqrts)*exp(-1.0*curr_set.rapidity)), 0.0}; 
			double xu[2] = { 1.0, 1.0};
				
			//Initialisation for the random number generator
			const gsl_rng_type *T;
			gsl_rng *r;
		
			//Initialisation for monte carlo
			gsl_monte_function G = { &frag_prod_monte, 2, &curr_set };
			
			gsl_rng_env_setup ();
			T = gsl_rng_default;
			
			r = gsl_rng_alloc (T);
			
			gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
			
			//Loop over a,b,c,d
			for(int a = -1*curr_set.nf;a <= curr_set.nf;a++) {
			for(int b = -1*curr_set.nf;b <= curr_set.nf;b++) {
			for(int c = -1*curr_set.nf;c <= curr_set.nf;c++) {	
			for(int d = -1*curr_set.nf;d <= curr_set.nf;d++) {
		
				//Determine the type
				type=two_two_type(a, b, c, d);
		
				//
				curr_set.parton_a=a;
				curr_set.parton_b=b;
				curr_set.parton_c=c;
				curr_set.parton_d=d;
				curr_set.type=type;
		
				result=0.0;
		
				//If type != 0
				if (type != 0) {
					
					//Integrate direct_prod over xa
					//cout << a << " + " << b << " + " << " -> " << c << " + " << d << " : " << type << "\n";
		
					//Set some max number of iteration
					maxiter=curr_set.maxiter;
					
					//Try to find the number of iterations required to achieve a precision of "relerr"
					//Don't care about the khi^2 right now
					do {
						gsl_monte_vegas_integrate (&G, xl, xu, 2, maxiter, r, s, &result, &error);
						
						if (error/result > relerr)  {
							maxiter*=1.5;
							//cout << "Max number of iter:" << maxiter << "\n";
						}
						
					} while (error/result > relerr);
					
					//cout << "Max number of iter:" << maxiter << "\n";
					
					//Try to get a khi-squared close to 1 
					//(error estimates may not be reliable for a khi-squared significantly different
					// of 1, although they generally are)
					do {
					
						do {
							gsl_monte_vegas_integrate (&G, xl, xu, 2, maxiter, r, s, &result, &error);
							
							//cout << "Chi-squared:" << s->chisq << "\n";
						
						} while (fabs(s->chisq - 1.0) > 0.25);
						
						if (error/result > relerr)  {
							maxiter*=1.5;
						}
					
					} while (error/result > relerr);
		
				}
				
				total+=result;
						
			}
		
			}
		
			}
		
			}
			
			gsl_monte_vegas_free (s);
			
			gsl_rng_free (r);
			
			total*=curr_set.units_conv*(*curr_set.Kfunct)(curr_set.pHt)*1.0/M_PI;
			
		}
		
		else {
		
		if ('q' == curr_set.frag_int_method) {
			
			
			//
			maxiter=curr_set.maxiter;
		 
			//
			gsl_function F;

	 		//initialisation for the numerical integration
			gsl_integration_workspace * w = gsl_integration_workspace_alloc(maxiter);
		 
			F.function = &frag_prod_qag;
			F.params = &curr_set;
		 
			lowbound=(2.0*curr_set.pHt/curr_set.sqrts)*(exp(curr_set.rapidity)+exp(-1.0*curr_set.rapidity))/2.0;
			upbound=1.0;
		 
			//
			abserr=curr_set.abserr;
			relerr=curr_set.relerr;
		 
			method=curr_set.qag_method;

			for(int a = -1*curr_set.nf;a <= curr_set.nf;a++) {
			for(int b = -1*curr_set.nf;b <= curr_set.nf;b++) {
			for(int c = -1*curr_set.nf;c <= curr_set.nf;c++) {
			for(int d = -1*curr_set.nf;d <= curr_set.nf;d++) {

				//Determine the type
				type=two_two_type(a, b, c, d);

				//
				curr_set.parton_a=a;
				curr_set.parton_b=b;
				curr_set.parton_c=c;
				curr_set.parton_d=d;
				curr_set.type=type;

				result=0.0;

				//If type != 0
				if (type != 0) {
		
					//Integrate direct_prod over xa
					gsl_integration_qag (&F, lowbound, upbound, abserr, relerr, maxiter, method, w, &result, &error); 

				}
	
				total+=result;

			}

			}
			
			}
	
			}

			//This function frees the memory associated with the workspace w. 
			gsl_integration_workspace_free(w);
		
			//
			total*=curr_set.units_conv*(*curr_set.Kfunct)(curr_set.pHt);
			
		}
		
		else {
			cout << "What method should be used for the numerical integration?\n";
			cout << "'m' for Monte Carlo (VEGAS), 'q' for adaptive quadrature\n";
			exit(1);
		}
		
		}
	
	}
	
	return total;

}


//Direct parton/photon production
double direct_prod (double xa, void * pars) {

	
	// Declaration of functions //
	double cteq5_pdf(int, int, double, double);
	double twotwocs(void *, double, double, double, double);
	void set_scales(double *, double *, void *);
	//double lhapdf(int,double,double);
	//double mod_pdf(int [], int, double , double);

	// Declaration of variables //
	params curr_set;
	double result;
	double Qren, Qfac;
	double xJt,xb,y,sqrts,pJt;
	double pdfa, pdfb, partoncs;
	int cteq_table;
	double sh, th, uh;

	// Initialization of variables //
	curr_set = * (params *) pars;
	
	//
	cteq_table=curr_set.cteq_table;
	
	//
	sqrts=curr_set.sqrts;
	y=curr_set.rapidity;
	pJt=curr_set.pJt;

	// Core //

	//
	set_scales(&Qfac, &Qren, &curr_set);
	
	//xJt
	xJt=2.0*pJt/curr_set.sqrts;
	
	//
	xb=xa*xJt*exp(-y)/(2.0*xa-xJt*exp(y));
	
	//Mandelstam variables for _massless_ partons
	sh=xa*xb*sqrts*sqrts;
	th=-xa*pJt*sqrts*exp(-y);
	uh=-xb*pJt*sqrts*exp(y);

	//pdfa
	pdfa=cteq5_pdf(curr_set.parton_a, cteq_table, xa, Qfac);
	//pdfa=lhapdf(curr_set.parton_a, xa, Qfac);
	//pdfa=mod_pdf(curr_set.beam_1, curr_set.parton_a, xa, Qfac);
	
	//pdf b
	pdfb=cteq5_pdf(curr_set.parton_b, cteq_table, xb, Qfac);
	//pdfb=lhapdf(curr_set.parton_b, xb, Qfac);
	//pdfb=mod_pdf(curr_set.beam_2, curr_set.parton_b, xb, Qfac);

	//if photon

		// ...

	//if hadron, compute initial hard jet spectra

		//...

	//parton-parton cs
	partoncs=twotwocs(&curr_set, sh, th, uh, Qren);

	//res
	result=pdfa*pdfb*xa*xb/(2.0*xa-xJt*exp(y))*partoncs;
	
	return result;

}

//Fragmentation hadron/photon production
double frag_prod_monte (double * x, size_t dim, void * pars)  {

	// Declaration of functions //
	double cteq5_pdf(int, int, double, double);
	double twotwocs(void *, double, double, double, double);
	double photons_ff(int, int, double, double);
	double kkp_ff(int, int, int, double, double);
	void set_scales(double *, double *, void *);
	//double lhapdf(int , double, double);
	//double mod_pdf(int [], int, double , double);

	// Declaration of variables //
	double result;
	double Qren, Qfac, Qfrag;
	double xHt,y,sqrts,pHt, pJt, zc;
	double sh, th, uh;
	double pdfa, pdfb, partoncs;
	double ff;
	params curr_set;
	int cteq_table, order, set;
	double xa, xb;
	
	// Initialization of variables //
	curr_set = * (params *) pars;
	
	//
	cteq_table=curr_set.cteq_table;
	order=curr_set.kkp_order;
	set=curr_set.photon_set;
	
	//
	sqrts=curr_set.sqrts;
	y=curr_set.rapidity;
	pHt=curr_set.pHt;

	//
	xa=x[0];
	xb=x[1];

	// Core //

	//xHt
	xHt=2.0*pHt/sqrts;
	
	if ((xb > xa*xHt*exp(-y)/(2.0*xa-xHt*exp(y)))&&(xa > xHt*exp(y)/(2.0-xHt*exp(-1.0*y)))) {

		//
		zc=xHt/2.0*(exp(-y)/xb+exp(y)/xa);
		
		//
		pJt=pHt/zc;
		
		curr_set.pJt=pJt;
	
		//
		set_scales(&Qfac, &Qren, &curr_set);
		Qfrag=pHt;
		
		//Mandelstam variables for _massless_ partons
		sh=xa*xb*sqrts*sqrts;
		th=-xa*pJt*sqrts*exp(-y);
		uh=-xb*pJt*sqrts*exp(y);
	
		//pdf a
		pdfa=cteq5_pdf(curr_set.parton_a, cteq_table, xa, Qfac);
		//pdfa=lhapdf(curr_set.parton_a, xa, Qfac);
		//pdfa=mod_pdf(curr_set.beam_1, curr_set.parton_a, xa, Qfac);
	
		//pdf b
		pdfb=cteq5_pdf(curr_set.parton_b, cteq_table, xb, Qfac);
		//pdfb=lhapdf(curr_set.parton_b, xb, Qfac);
		//pdfb=mod_pdf(curr_set.beam_2, curr_set.parton_b, xb, Qfac);
		
		//if photon
		if (10 == curr_set.part_type) {
			//ff = ...
			//kkp_ff(int hadron, int parton, int order, double x, double scale)
			//photons_ff(int parton, int set, double x, double scale) 
			ff=photons_ff(curr_set.parton_c, set, zc, Qfrag);
		}
		else {
			//if hadron
			if ((curr_set.part_type >= 1)&&(curr_set.part_type <= 7)) {
				//ff = ...
				//kkp_ff(int hadron, int parton, int order, double x, double scale)
				ff=kkp_ff(curr_set.part_type, curr_set.parton_c, order, zc, Qfrag);
			}
			else {
				cout << "Unknown particle type in frag_prod_monte(). Aborting...";
				exit(1);
			}
		}
		//parton-parton cs
		partoncs=twotwocs(&curr_set, sh, th, uh, Qren);
	
		//res
		result=pdfa*pdfb*ff*partoncs*1.0/zc;
		
	}
	else {
		
		result=0;
	}
	
	return result;

}

//Fragmentation hadron/photon production
double frag_prod_qag (double zc, void * pars)  {

	// Declaration of functions //
	double photons_ff(int, int, double, double);
	double kkp_ff(int, int, int, double, double);
	double frag_prod_qag_sub (double, void *);

	// Declaration of variables //
	double result;
	double Qfrag;
	double xHt,y,sqrts,pHt,pJt,xJt;
	double ff, direct_res, direct_err;
	params curr_set, tmpset;
	
	size_t maxiter;
	double abserr, relerr;
	int method;
	double lowbound, upbound;
	
	// Initialization of variables //
	curr_set = * (params *) pars;
	tmpset=curr_set;
	
	//
	sqrts=curr_set.sqrts;
	y=curr_set.rapidity;
	pHt=curr_set.pHt;
	
	//
	maxiter=curr_set.maxiter;
	abserr=curr_set.abserr;
	relerr=curr_set.relerr;
	method=curr_set.qag_method;
	
	// Core //

	//xHt
	xHt=2.0*pHt/sqrts;
	
	if (zc > (xHt*(exp(y)+exp(-y))/2.0)) {
		
		gsl_function F;

		//initialisation for the numerical integration
		gsl_integration_workspace * w = gsl_integration_workspace_alloc(maxiter);
	
		F.function = &frag_prod_qag_sub;
		F.params = &tmpset;

		//
		pJt=pHt/zc;
		xJt=2.0*pJt/sqrts;
		
		tmpset.pJt=pJt;
		tmpset.pHt=pHt;
		
		//
		lowbound=(xJt*exp(y)/(2.0-xJt*exp(-y)));
		upbound=1.0;
	
		//
		Qfrag=pHt;
		
		//if photon
		if (10 == curr_set.part_type) {
			//ff = ...
			//kkp_ff(int hadron, int parton, int order, double x, double scale)
			//photons_ff(int parton, int set, double x, double scale) 
			ff=photons_ff(curr_set.parton_c, curr_set.photon_set, zc, Qfrag);
		}
		else {
			//if hadron
			if ((curr_set.part_type >= 1)&&(curr_set.part_type <= 7)) {
				//ff = ...
				//kkp_ff(int hadron, int parton, int order, double x, double scale)
				ff=kkp_ff(curr_set.part_type, curr_set.parton_c, curr_set.kkp_order, zc, Qfrag);
			}
			else {
				cout << "Unknown particle type in frag_prod_qag(). Aborting...";
				exit(1);
			}
		}
		gsl_integration_qag (&F, lowbound, upbound, abserr, relerr, maxiter, method, w, &direct_res, &direct_err); 

		//res
		result=ff*direct_res/zc/zc;

		//This function frees the memory associated with the workspace w. 
		gsl_integration_workspace_free(w);

		
	}
	else {
		
		result=0;
	}

	return result;

}

//Fragmentation hadron/photon production
double frag_prod_qag_sub (double xa, void * pars)  {

	// Declaration of functions //
	double cteq5_pdf(int, int, double, double);
	double twotwocs(void *, double, double, double, double);
	void set_scales(double *, double *, void *);
	//double lhapdf(int , double, double);
	//double mod_pdf(int [], int, double , double);

	// Declaration of variables //
	double result;
	double xJt,y,sqrts, pJt;
	double sh, th, uh;
	double pdfa, pdfb, partoncs;
	params curr_set;
	int cteq_table;
	double xb;
	double qren, qfac;
	
	// Initialization of variables //
	curr_set = * (params *) pars;
	
	//
	cteq_table=curr_set.cteq_table;
	
	//
	sqrts=curr_set.sqrts;
	y=curr_set.rapidity;
	pJt=curr_set.pJt;

	// Core //

	//
	set_scales(&qfac, &qren, &curr_set);
	
	//xHt
	xJt=2.0*pJt/sqrts;
	
	if (xa > (xJt*exp(y)/(2.0-xJt*exp(-y)))) {

		//
		xb=xa*xJt*exp(-y)/(2.0*xa-xJt*exp(y));
		
		//Mandelstam variables for _massless_ partons
		sh=xa*xb*sqrts*sqrts;
		th=-xa*pJt*sqrts*exp(-y);
		uh=-xb*pJt*sqrts*exp(y);
	
		//pdf a
		pdfa=cteq5_pdf(curr_set.parton_a, cteq_table, xa, qfac);
		//pdfa=lhapdf(curr_set.parton_a, xa, qfac);
		//pdfa=mod_pdf(curr_set.beam_1, curr_set.parton_a, xa, qfac);
	
		//pdf b
		pdfb=cteq5_pdf(curr_set.parton_b, cteq_table, xb, qfac);
		//pdfb=lhapdf(curr_set.parton_b, xb, qfac);
		//pdfb=mod_pdf(curr_set.beam_2, curr_set.parton_b, xb, qfac);
		
		//parton-parton cs
		partoncs=twotwocs(&curr_set, sh, th, uh, qren);
	
		//res
		result=2.0/M_PI*pdfa*pdfb*partoncs*xa*xb/(2.0*xa-xJt*exp(y));
		
	}
	else {
		
		result=0;
	}
	
	return result;

}
