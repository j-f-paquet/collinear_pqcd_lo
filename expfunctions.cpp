#include <cmath>
#include <iostream>
#include <stdlib.h>

//#include "LHAPDF/LHAPDF.h"

//p.d.f.

/*
   
The function Ctq5Pdf (Iparton, X, Q)
   returns the parton distribution inside the proton for parton [Iparton] 
   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
              CTEQ5F4 has only 4 flavors and gluon.

*/

extern "C" {
	
	double ctq5pdf_(int*, double*, double*);
	void setctq5_(int*);

}


double cteq5_pdf(int parton, int table, double x, double scale) {
	 
	//Declaration of functions//
	
	//Declarations of variables//
	double pdfval;
	
	//Initialisation of variables//
	
	//Core//
	
	//Checks
	//10^-5 < x < 1 and 1.0 < Q < 10,000 (GeV)
	if ((x < 1E-5)||(x>1.0)) {
		std::cout << "Momentum fraction outside the range: " << x << " ...\n";
		exit(1);
	}
	if ((scale < 1.0)||(scale>1.0E4)) {
		std::cout << "Scale outside the range: " << scale << " ...\n";
		exit(1);
	}
	if ((parton > 5)||(parton < -5)) {
		std::cout << "Not a valid parton: " << parton << " ...\n";
		exit(1);
	}
	if ((table > 9)||(table < 1)) {
		std::cout << "Not a valid table: " << table << " ...\n";
		exit(1);
	}
	
	//Select the table
	setctq5_(&table);
	
	//
	pdfval=ctq5pdf_(&parton,&x,&scale);
	
	return pdfval;
	
}

//f.f.

extern "C" {
	void kkp_(int*,int*,double*,double*,double[]);
}


//subroutine kkp(ih,iset,x,qs,dh)
double kkp_ff(int hadron, int parton, int order, double x, double scale) {
	
	//Declaration of functions//
	
	//Declarations of variables//
	double ffval, fflist[11];
	int kkp_parton,kkp_hadron;
	
	//Initialisation of variables//
	
	//Core//
	
	//Checks
	if ((order != 0)&&(order != 1)) {
		std::cout << "Not a valid order: " << order << " ...\n";
		exit(1);
	}
	if ((x < 0.1)||(x>0.8)) {
		//std::cout << "Momentum fraction outside the range: " << x << " ...\n";
		//exit(1);
	}
	if ((parton < -5)||(parton > 5)) {
		std::cout << "Not a valid parton: " << parton << " ...\n";
		exit(1);
	}
	if ((scale > 100)||((scale < 1.414213562)&&(abs(parton) <= 3))||((scale < 2.9788)&&(4 == abs(parton)))||((scale < 9.46037)&&(5 == abs(parton)))) {
		//std::cout << "Scale outside the range: " << scale << " ...\n";
		//exit(1);
	}

	//
	kkp_hadron=hadron;

	//kkp_parton
	switch(parton) {
			
		case -5:
			kkp_parton=10;
			break;
		case -4:
			kkp_parton=8;
			break;
		case -3:
			kkp_parton=6;
			break;
		case -2:
			kkp_parton=4;
			break;
		case -1:
			kkp_parton=2;
			break;	
		case 0:
			kkp_parton=0;
			break;
		case 1:
			kkp_parton=1;
			break;
		case 2:
			kkp_parton=3;
			break;
		case 3:
			kkp_parton=5;
			break;
		case 4:
			kkp_parton=7;
			break;
		case 5:
			kkp_parton=9;
			break;
			
		default:
			std::cout << "Unknown parton in lhapdf(). Aborting...\n";
			exit(1);
			break;

	}
	
	kkp_(&kkp_hadron,&order,&x,&scale,fflist);
	
	ffval=fflist[kkp_parton];
	
	return ffval;

}


//photon f. f.

extern "C" {

	void fonfra_(double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
}

/*
      SUBROUTINE FONFRA(X,IFLAG,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP,XDBP
# ,XDBBP,XDTP,XDTBP,XDGP)
*/

double photons_ff(int parton, int set, double x, double scale) {
	
	//Declaration of functions//
	
	//Declarations of variables//
	double xu,xub,xd,xdb,xs,xsb;
	double xc,xcb,xb,xbb,xt,xtb,xg;
	double ffval,scale2;
	
	//Initialisation of variables//

	
	//Core//
	
	//Checks
	if ((set != 1)&&(set != 2)) {
		std::cout << "Not a valid set: " << set << " ...\n";
		exit(1);
	}
	if ((x <= 0)||(x > 1)) {
		std::cout << "Momentum fraction outside the range: " << x << " ...\n";
		exit(1);
	}
	if ((parton < -5)||(parton > 5)) {
		std::cout << "Not a valid parton: " << parton << " ...\n";
		exit(1);
	}
	/*if ((scale > 100)||((scale < 1.414213562)&&(abs(parton) <= 3))||((scale < 2.9788)&&(4 == abs(parton)))||((scale < 9.46037)&&(5 == abs(parton)))) {
		//std::cout << "Scale outside the range: " << scale << " ...\n";
		//exit(1);
}*/

	scale2=pow(scale,2);
	
	//
	fonfra_(&x,&set,&scale2,&xu,&xub,&xd,&xdb,&xs,&xc,&xb,&xbb,&xt,&xtb,&xg);
	
	//kkp_parton
	switch(parton) {
			
		case -5:
			ffval=xbb;
			break;
		case -4:
			ffval=xc;
			break;
		case -3:
			ffval=xs;
			break;
		case -2:
			ffval=xdb;
			break;
		case -1:
			ffval=xub;
			break;	
		case 0:
			ffval=xg;
			break;
		case 1:
			ffval=xu;
			break;
		case 2:
			ffval=xd;
			break;
		case 3:
			ffval=xs;
			break;
		case 4:
			ffval=xc;
			break;
		case 5:
			ffval=xb;
			break;
			
		default:
			std::cout << "Unknown parton in lhapdf(). Aborting...\n";
			exit(1);
			break;

	}
	
	ffval/=x;
	
	return ffval;

}

//double lhapdf(int parton, double x, double scale) {
//	
//	//Declaration of functions//
//	
//	//Declarations of variables//
//	double ffval;
//	int lhapdf_parton;
//	
//	//Initialisation of variables//
//	
//	//Core//
//	
//	//kkp_parton
//	switch(parton) {
//			
//		case -5:
//			lhapdf_parton=-5;
//			break;
//		case -4:
//			lhapdf_parton=-4;
//			break;
//		case -3:
//			lhapdf_parton=-3;
//			break;
//		case -2:
//			lhapdf_parton=-1;
//			break;
//		case -1:
//			lhapdf_parton=-2;
//			break;	
//		case 0:
//			lhapdf_parton=0;
//			break;
//		case 1:
//			lhapdf_parton=2;
//			break;
//		case 2:
//			lhapdf_parton=1;
//			break;
//		case 3:
//			lhapdf_parton=3;
//			break;
//		case 4:
//			lhapdf_parton=4;
//			break;
//		case 5:
//			lhapdf_parton=5;
//			break;
//
//		default:
//			std::cout << "Unknown parton in lhapdf(). Aborting...\n";
//			exit(1);
//			break;
//
//	}
//	
//	ffval=LHAPDF::xfx(x, scale, lhapdf_parton)/x;
//	
//	return ffval;
//
//}


//There seem to be no equivalent to "call setlhaparm('EPS08')" in the C++ wrapper, 
//so I assume that xfxa() is using EKS98 
//(it is actually important to know, since EKS98 and EPS08 are not using the same "reference" p.d.f.)

//I also assume that xfxa() is not already taking care of isospin symmetry
//(seems sensible since we are only supplying 1 param, "a")

//u_val(neutron)=d_val(proton)=d(proton)-dbar(proton)
//d_val(neutron)=u_val(proton)=u(proton)-ubar(proton)
//u_sea(neutron)=/ u_sea(proton)=ubar(proton) [if the sea quarks are the same for neutrons and protons]
//              | or
//              \ d_sea(proton)=dbar(proton) [if they are also changed by isospin symmetry (so that u(neutron)=d(proton))]
//d_sea(neutron)=/ d_sea(proton)=dbar(proton) [idem as above]
//              | or
//              \ u_sea(proton)=ubar(proton)  [idem as above]
//s(neutron)=s_sea(neutron)=s(proton), c(neutron)=c(proton), ...

//double mod_pdf(int beam[2], int parton, double x, double scale) {
//
//	//Declaration of functions//
//	
//	//Declarations of variables//
//	int nb_neutron, nb_proton,total;
//	double res, pdflist[13];
//	double tmp1,tmp2;
//	
//	//Initialisation of variables//
//	nb_proton=beam[0];
//	nb_neutron=beam[1];
//	total=nb_neutron+nb_proton;
//	
//	//Core//
//	
//	LHAPDF::xfxa(x, scale, total, pdflist);
//
//	switch(parton) {
//		
//		case -5:
//		case -4:
//		case -3:
//		case 0:
//		case 3:
//		case 4:
//		case 5:
//			res=pdflist[parton+6]/x;
//			break;
//		case -2:
//			//Note that LHAPDF uses 1 for down and 2 for up
//			tmp1=pdflist[-1+6]/x;
//			tmp2=pdflist[-2+6]/x;
//			//tmp2=pdflist[-1+6]/x;
//			res=nb_proton/double(total)*tmp1+nb_neutron/double(total)*tmp2;
//			break;
//		case -1:
//			tmp1=pdflist[-2+6]/x;
//			tmp2=pdflist[-1+6]/x;
//			//tmp2=pdflist[-2+6]/x;
//			res=nb_proton/double(total)*tmp1+nb_neutron/double(total)*tmp2;
//			break;	
//		case 1:
//			tmp1=pdflist[2+6]/x;
//			tmp2=pdflist[1+6]/x;
//			//tmp2=(pdflist[1+6]-pdflist[-1+6]+pdflist[-2+6])/x;
//			res=nb_proton/double(total)*tmp1+nb_neutron/double(total)*tmp2;
//			break;
//		case 2:
//			tmp1=pdflist[1+6]/x;
//			tmp2=pdflist[2+6]/x;
//			//tmp2=(pdflist[2+6]-pdflist[-2+6]+pdflist[-1+6])/x;
//			res=nb_proton/double(total)*tmp1+nb_neutron/double(total)*tmp2;
//			break;
//			
//		default:
//			std::cout << "Unknown parton in mod_pdf(). Aborting...\n";
//			exit(1);
//			break;
//
//	}
//	
//	return res;
//	
//}
