#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "gsl/gsl_math.h"

#include "params.h"

//alpha
double alpha() {

	// Declaration of functions //

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	res=7.297352e-3;

	return res;

}



//alpha_s
double alpha_s(double q, double lambda, int nf) {

	// Declaration of functions //

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	//From arXiv:0802.4364, eq 3.2
	res=12.0*M_PI/( (33.0-2.0*nf) * log( pow(q,2) / pow(lambda,2) ) );
	
	return res;


}


//parton-parton cross-section
//could probably be optimised to reduce the number of if/then/else required, but I'm pretty sure that the code would then be unreadable
int two_two_type(int parton_a, int parton_b, int parton_c, int parton_d) {

	// Declaration of functions //

	// Declaration of variables //
	int type;

	// Initialisation of variables //

	type = 0;

	// Core //

	//Processes with no gluons (1,2,3,4 and maybe 11,12)
	if ((parton_a != 0)&&(parton_b != 0)&&(parton_c != 0)&&(parton_d != 0)) {

		//type 1
		if ((parton_a == parton_c)&&(parton_b == parton_d)&&(abs(parton_a) != abs(parton_b))) {
			type=1;
		}
		else {

		//type -1
		if ((parton_a == parton_d)&&(parton_b == parton_c)&&(abs(parton_a) != abs(parton_b))) {
			type=-1;
		}
		else {
			
		//type 2
		if ((parton_a == parton_b)&&(parton_b == parton_c)&&(parton_c == parton_d)) {
			type=2;
		}
		else {

		//type 3
		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(abs(parton_a) != abs(parton_c))) {
			type=3;
		}
		else {

		//type 4
		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(parton_a == parton_c)) {
			type=4;
		}
		else {

		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(parton_a == parton_d)) {
			type=-4;
		}

		}

		}

		}

		}
		
		}

	}
	else {

		//type 5
		if (((parton_a == 0)&&(parton_c == 0)&&(parton_b == parton_d)&&(parton_b != 0))||((parton_b == 0)&&(parton_d == 0)&&(parton_a == parton_c)&&(parton_a != 0))) {
			type=5;
		}
		else {
	
		//type -5
		if (((parton_a == 0)&&(parton_d == 0)&&(parton_b == parton_c)&&(parton_b != 0))||((parton_b == 0)&&(parton_c == 0)&&(parton_a == parton_d)&&(parton_a != 0))) {
			type=-5;
		}
		else { 
	
		//type 6
		if ((parton_a == -1*parton_b)&&(parton_c == 0)&&(parton_d == 0)&&(parton_a != 0)) {
			type=6;
		}
		else {
	
		//type 7
		if ((parton_c == -1*parton_d)&&(parton_a == 0)&&(parton_b == 0)&&(parton_c != 0)) {
			type=7;
		}
		else {
	
		//type 8
		if ((parton_a == 0)&&(parton_b == 0)&&(parton_c == 0)&&(parton_d == 0)) {
			type=8;
		}
		else {
	
		//type 9
		if (((parton_a == 0)&&(parton_b == parton_d)&&(parton_b != 0)&&(parton_c == 10))||((parton_b == 0)&&(parton_a == parton_c)&&(parton_a != 0)&&(parton_d == 10))) {
			type=9;
		}	
		else {
	
		//type -9
		if (((parton_a == 0)&&(parton_b == parton_c)&&(parton_b != 0)&&(parton_d == 10))||((parton_b == 0)&&(parton_a == parton_d)&&(parton_a != 0)&&(parton_c == 10))) {
			type=-9;
		}
		else {

		//type 10
		if (((parton_a == -1*parton_b)&&(parton_a != 0))&&(((parton_c == 10)&&(parton_d == 0))||((parton_c == 10)&&(parton_d == 0)))) {
			type=10;
		}
		else {
			type=0;
		}

		}
	
		}
	
		}
	
		}
	
		}
	
		}
	
		}
	
	}

	return type;

}



//type
double twotwocs(void * pars, double s, double t, double u, double scale) {

	/* Declaration of functions */
	double charge(int, int, int, int, int);

	/* Declaration of variables */
	double res;
	double sh, th, uh;
	double shs, ths, uhs;
	double valpha, valphas, vcharge;
	params plist;
	int type;


	/* Initialisation of variables */
	plist= *(params *) pars;

	sh=s;
	if (plist.type < 0) {
		uh=t;
		th=u;
		type=-1*plist.type;
	}
	else {
		th=t;
		uh=u;
		type=plist.type;
	}

	shs=sh*sh;
	ths=th*th;
	uhs=uh*uh;

	valpha=alpha();
	valphas=alpha_s(scale, plist.lambda, plist.nf);

	/* Core */

	//charge
	if ((9 == type)||(10 == type)) {
		vcharge=charge(plist.type, plist.parton_a, plist.parton_b, plist.parton_c, plist.parton_d);
	}	

	//From Owens, table I
	switch (type) {

		case 1:
			res = valphas*valphas*(4.0/9.0* ( shs+uhs ) /ths);
			//group of statements 1;
			break;
		case 2:
			//group of statements 2;
			res = valphas*valphas*(4.0/9.0* ( ( shs+uhs ) /ths + ( shs+ths ) /uhs )-8.0/27.0*shs/ ( th*uh ));
			break;

		case 3:
			res = valphas*valphas*(4.0/9.0* ( ths+uhs ) /shs);
			break;

		case 4:
			res = valphas*valphas*4.0/9.0*( ( ( shs+uhs ) /ths+ ( uhs+ths ) /shs )-8.0/27.0*uhs/ ( th*sh ) );
			break;

		case 5:
			res = valphas*valphas*(-4.0/9.0* ( sh/uh+uh/sh ) + ( shs+uhs ) /ths);
			break;

		case 6:
			res = valphas*valphas*(32.0/27.0* ( th/uh+uh/th )-8.0/3.0* ( ths+uhs ) /shs);
			break;

		case 7:
			res = valphas*valphas*(1.0/6.0* ( th/uh+uh/th )-3.0/8.0* ( ths+uhs ) /shs);
			break;

		case 8:
			res = valphas*valphas*(9.0/2.0* ( 3.0-th*uh/shs-sh*uh/ths-sh*th/uhs ));
			break;
		
		case 9:
			res = valpha*valphas*(-1.0/3.0*(vcharge*vcharge)*(uh/sh+sh/uh));
			break;

		case 10:
			res = valpha*valphas*(8.0/9.0*(vcharge*vcharge)*(uh/th+th/uh));
			break;

		default:
			std::cout << "Unknown type in twotwocs(): " << plist.type << " Aborting.../n";
			exit(1);
			break;

	}



	res=M_PI*res/shs;

	return res;

}


//charge
double charge(int type, int parton_a, int parton_b, int parton_c, int parton_d) {

	// Declaration of functions //

	// Declaration of variables //
	double res;
        int quark;

	// Initialization of variables //
	quark=0;
	
	// Core //
	
	//only need the charge for processes 9 and 10
	if (9 == abs(type)) {
		if (parton_a == 0) {
			quark=parton_b;
		}
		else {
			quark = parton_a;
		}
	}
	else {
		if (10 == type) {
			if (parton_a > 0) {
				quark=parton_a;
			}
			else {
				quark = parton_b;
			}
		}
	}

	//charge of each quark flavour
	switch(abs(quark)) {
		
		case 1:
			res = 2.0/3.0;
			break;

		case 2:
			res = -1.0/3.0;
			break;
			
		case 3:
			res = -1.0/3.0;
			break;

		case 4:
			res = 2.0/3.0;
			break;
			
		case 5:
			res = -1.0/3.0;
			break;

		case 6:
			res = 2.0/3.0;
			break;

		default:
			res = 0;
			break;
		
	}

	//if antiquark, -1*charge
	if (type < 0) {
		res*=-1;
	}
	
	return res;

}
