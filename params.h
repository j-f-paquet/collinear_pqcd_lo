#include <stdlib.h>

struct params {

	//beams
	//beam_1/2[# of proton,# of neutron]
	//beam_1/2[1,0] for a proton, beam_1/2[1,0] for a neutron, beam_1/2[1,1] for deuteron, ...
	int beam_1[2], beam_2[2];

	//particle type
	int part_type;

	//source
	int source;

	//pT
	double pJt, pHt;

	//rapidity
	double rapidity;

	//sqrts
	double sqrts;

	//parton a, b, c, d
	int parton_a, parton_b, parton_c, parton_d;

	//type
	int type;

	//lambda, nf
	double lambda;
	int nf;

	//K-factors
	double Kfactor;
	double (*Kfunct)(double);

	//unit conversion factor
	double units_conv;

	//scales
	char Qfac, Qren; //h for hadron, j for jet
	double Qfac_norm, Qren_norm;

	//integration params
	size_t maxiter;
	double abserr, relerr;
	int qag_method;
	char frag_int_method;

	//exp functions
	int kkp_order, cteq_table, photon_set;
	char pdf_name[200];
	//bool pdf_init;

};