KKP_LIB_DIR=./kkp/

PHOTON_FF_LIB_DIR=./photon_ff_fnal/

CTEQ5_LIB_DIR=./Ctq5Pdf/

ghloprod : functions.o expfunctions.o utils.o ghloprod.cpp
	g++ -Wall -o ghloprod -g3 -pg ${KKP_LIB_DIR}/kkp.o ${CTEQ5_LIB_DIR}//Ctq5Pdf.o ${PHOTON_FF_LIB_DIR}//fonfra_inter2.o ${PHOTON_FF_LIB_DIR}//distribPert.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setI.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setII.o  ${PHOTON_FF_LIB_DIR}//locate.o ${PHOTON_FF_LIB_DIR}//polin2.o ${PHOTON_FF_LIB_DIR}//polint.o utils.o functions.o expfunctions.o ghloprod.cpp -lgsl -lblas64 -lm -lgfortran

expfunctions.o : expfunctions.cpp
	g++ -c -g -pg expfunctions.cpp 

#main.o : ghloprod.cpp
#	g++ -c -g ghloprod.cpp -o main.o

functions.o : functions.cpp
	g++ -c -Wall -g -pg functions.cpp

utils.o : utils.cpp
	g++ -c -Wall -g -pg utils.cpp

#optim : functions_optim.o expfunctions_optim.o utils_optim.o
#	g++ -Wall -o ghloprod_optim -O3 ${KKP_LIB_DIR}/kkp_optim.o ${CTEQ5_LIB_DIR}//Ctq5Pdf_optim.o ${PHOTON_FF_LIB_DIR}//fonfra_inter2.o ${PHOTON_FF_LIB_DIR}//distribPert.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setI.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setII.o  ${PHOTON_FF_LIB_DIR}//locate.o ${PHOTON_FF_LIB_DIR}//polin2.o ${PHOTON_FF_LIB_DIR}//polint.o utils.o functions.o expfunctions.o ghloprod.cpp -lgsl -lblas -lg2c -lm 
#
#expfunctions_optim.o : expfunctions.cpp
#	g++ -c -Wall -O3 expfunctions.cpp 
#
#functions_optim.o : functions.cpp
#	g++ -c -Wall -O3 functions.cpp
#
#utils_optim.o : utils.cpp
#	g++ -c -Wall -O3 utils.cpp
