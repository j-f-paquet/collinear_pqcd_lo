#include <iostream>
#include <fstream>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "params.h"

using namespace std;

void set_scales(double * qfac, double * qren, void * pars) {

	// Declaration of functions //

	// Declaration of variables //
	params curr_set;
	
	// Initialization of variables //
	curr_set = * (params *) pars;
	
	// Core //

	if ('h' == curr_set.Qfac) {
		*qfac=curr_set.pHt*curr_set.Qfac_norm;
	}
	else {
		if ('j' == curr_set.Qfac) {
			*qfac=curr_set.pJt*curr_set.Qfac_norm;
		}
		else {
			cout << "Unknown factorisation scale. Aborting...\n";
			exit(1);
		}
	}
	
	if ('h' == curr_set.Qren) {
		*qren=curr_set.pHt*curr_set.Qren_norm;
	}
	else {
		if ('j' == curr_set.Qren) {
			*qren=curr_set.pJt*curr_set.Qren_norm;
		}
		else {
			cout << "Unknown renormalisation scale. Aborting...\n";
			exit(1);
		}
	}

}

void save_info(ofstream * file, params * pars) {

	// Declaration of functions //

	// Declaration of variables //
	//params curr_set;
	//ofstream outfile;
	
	// Initialization of variables //

	// Core //
	
	*file << "### General params ### \n";
	*file << "#lambda: " << pars->lambda << "\n";
	*file << "#sqrt(s): " << pars->sqrts << "\n";
	*file << "#rapidity: " << pars->rapidity << "\n";
	*file << "#unit conv factor: " << pars->units_conv << "\n";
	*file << "#\n";
	
	//
	*file << "### Particles ### \n";
	*file << "#Beam 1: " << (pars->beam_1)[0] << " protons and " << (pars->beam_1)[1] << " neutrons\n";
	*file << "#Beam 2: " << (pars->beam_2)[0] << " protons and " << (pars->beam_2)[1] << " neutrons\n";
	*file << "#Particle type: " << pars->part_type << "\n";
	*file << "#Particle source: " << pars->source << "\n";
	*file << "#\n";
	
	//
	*file << "### Params ### \n";
	*file << "#Qfac: " << pars->Qfac_norm << "*";
	if ('j' == pars->Qfac) {
		*file << "pJt\n";
	}
	else {
		if ('h' == pars->Qfac) {
			*file << "pHt\n";
		}
		else {
			*file << "[error]\n";
		}
	}
	*file << "#Qren: " << pars->Qren_norm << "*";
	if ('j' == pars->Qren) {
		*file << "pJt\n";
	}
	else {
		if ('h' == pars->Qren) {
			*file << "pHt\n";
		}
		else {
			*file << "[error]\n";
		}
	}
	if (0 == pars->source) {
		*file << "#QAG method: " << pars->qag_method << "\n";
	}
	else {
		if (1 == pars->source) {
			
			*file << "#Integration method: " << (*pars).frag_int_method << "\n";
			
			if ('q' == pars->frag_int_method) {
				*file << "#QAG method: " << pars->qag_method << "\n";
			}
			
			
			
		}
	}
		
	*file << "#Max num of iter: " << pars->maxiter << "\n";
	*file << "#Abs err: " << pars->abserr << "\n";
	*file << "#Rel err: " << pars->relerr << "\n";
	
	//
	//*file << "### Exp functions ### \n";
	//*file << "#PDF name: " << pars->pdf_name << "\n";

	if (1 == pars->source) {
		
		if ((pars->part_type <= 7)&&(pars->part_type >= 1)) {
			
			*file << "#KKP order: " << pars->kkp_order << "\n";
			
		}
		
		
		
		if (10 == pars->part_type) {
			
			*file << "#Photon FF set: " << pars->photon_set << "\n";
			
		}
	}

}

double unit_fct(double q) {
        return 1.0;
}

double kdirect(double q) {

	//Declaration of functions//
	
	//Declarations of variables//
	gsl_interp * intobj;
	const gsl_interp_type * inttype = gsl_interp_cspline;
	size_t nb = 190;
	int ident;
	gsl_interp_accel * intaccel;
	double res;

	double x[] = {1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42,42.5,43,43.5,44,44.5,45,45.5,46,46.5,47,47.5,48,48.5,49,49.5,50,50.5,51,51.5,52,52.5,53,53.5,54,54.5,55,55.5,56,56.5,57,57.5,58,58.5,59,59.5,60,60.5,61,61.5,62,62.5,63,63.5,64,64.5,65,65.5,66,66.5,67,67.5,68,68.5,69,69.5,70,70.5,71,71.5,72,72.5,73,73.5,74,74.5,75,75.5,76,76.5,77,77.5,78,78.5,79,79.5,80,80.5,81,81.5,82,82.5,83,83.5,84,84.5,85,85.5,86,86.5,87,87.5,88,88.5,89,89.5,90,90.5,91,91.5,92,92.5,93,93.5,94,94.5,95,95.5};
	double y[] = {1.5451830,1.8500339,1.9858160,1.8377326,1.7494822,1.6924863,1.6513450,1.6208120,1.5965786,1.5765117,1.5615237,1.5485356,1.5372588,1.5284614,1.5193268,1.5132944,1.5067669,1.5024014,1.4977407,1.4942148,1.4906463,1.4878970,1.4857789,1.4828951,1.4820492,1.4804452,1.4788713,1.4786246,1.4773642,1.4770257,1.4770175,1.4762853,1.4766084,1.4768701,1.4771131,1.4772723,1.4781454,1.4786891,1.4790142,1.4804800,1.4813783,1.4818605,1.4833702,1.4848481,1.4856152,1.4870465,1.4887775,1.4899971,1.4912763,1.4932795,1.4951041,1.4965948,1.4980753,1.5004132,1.5022117,1.5040446,1.5057727,1.5082706,1.5103142,1.5122161,1.5143650,1.5168418,1.5192011,1.5212513,1.5236088,1.5262320,1.5288522,1.5311524,1.5334733,1.5363680,1.5391536,1.5416636,1.5440105,1.5472074,1.5500706,1.5527842,1.5553302,1.5586582,1.5616964,1.5645515,1.5673423,1.5707992,1.5741752,1.5771062,1.5798337,1.5834996,1.5871341,1.5905715,1.5935158,1.5966657,1.6007658,1.6044671,1.6076685,1.6108328,1.6151439,1.6189865,1.6227675,1.6259776,1.6301226,1.6345927,1.6384600,1.6421139,1.6457922,1.6507530,1.6550202,1.6589134,1.6627394,1.6676943,1.6725191,1.6770470,1.6809577,1.6853555,1.6910151,1.6960567,1.7001801,1.7041689,1.7101886,1.7159787,1.7206474,1.7251203,1.7301603,1.7369451,1.7426433,1.7473982,1.7516842,1.7587400,1.7655703,1.7714500,1.7760683,1.7807620,1.7893640,1.7966810,1.8025639,1.8076285,1.8126774,1.8228048,1.8316361,1.8389739,1.8446655,1.8490775,1.8562999,1.8664868,1.8753288,1.8833209,1.8903235,1.8981563,1.9056605,1.9144396,1.9232488,1.9322082,1.9410062,1.9507982,1.9594735,1.9678717,1.9784420,1.9883042,1.9974059,2.0080764,2.0171137,2.0282766,2.0377333,2.0490153,2.0597273,2.0710639,2.0824778,2.0935676,2.1032202,2.1148486,2.1256090,2.1451098,2.1492553,2.1633024,2.1784391,2.1890581,2.2016011,2.2181809,2.2297393,2.2426884,2.2603271,2.2779565,2.2966523,2.3148218,2.3304423,2.3590673,2.3810545,2.4166506,2.4540165,2.4944825,2.5394321,2.6222995};
	
	//Initialisation of variables//
	
	//Core//
	
	//
	intobj=gsl_interp_alloc (inttype, nb);
	ident=gsl_interp_init (intobj, x, y, nb);
	intaccel=gsl_interp_accel_alloc();

	res=gsl_interp_eval (intobj, x, y, q, intaccel);

	//
	gsl_interp_free (intobj);
	gsl_interp_accel_free (intaccel);

	return res;

}

double kfrag(double) {

	// Declaration of functions //

	// Declaration of variables //
	
	// Initialization of variables //
	
	// Core //

	return 2.8;

}
