KKP_LIB_DIR=~/lib/kkp/

LHAPDF_LIB_DIR=~/lib/lhapdf/

PHOTON_FF_LIB_DIR=~/lib/photon_ff_fnal/

CTEQ5_LIB_DIR=~/lib/Ctq5Pdf/

ghloprod : functions.o expfunctions.o utils.o
	g++ -Wall -o ghloprod -g -pg ${KKP_LIB_DIR}/kkp_optim.o ${CTEQ5_LIB_DIR}//Ctq5Pdf_optim.o ${PHOTON_FF_LIB_DIR}//fonfra_inter2.o ${PHOTON_FF_LIB_DIR}//distribPert.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setI.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setII.o  ${PHOTON_FF_LIB_DIR}//locate.o ${PHOTON_FF_LIB_DIR}//polin2.o ${PHOTON_FF_LIB_DIR}//polint.o utils.o functions.o expfunctions.o src/ghloprod.cpp -lgsl -lblas -lg2c -lm -lLHAPDF -L${LHAPDF_LIB_DIR}//lib -I${LHAPDF_LIB_DIR}//include/

expfunctions.o : src/expfunctions.cpp
	g++ -c -g -pg src/expfunctions.cpp -I${LHAPDF_LIB_DIR}//include/

#main.o : src/ghloprod.cpp
#	g++ -c -g src/ghloprod.cpp -o main.o

functions.o : src/functions.cpp
	g++ -c -Wall -g -pg src/functions.cpp

utils.o : src/utils.cpp
	g++ -c -Wall -g -pg src/utils.cpp

optim : functions_optim.o expfunctions_optim.o utils_optim.o
	g++ -Wall -o ghloprod_optim -O3 ${KKP_LIB_DIR}/kkp_optim.o ${CTEQ5_LIB_DIR}//Ctq5Pdf_optim.o ${PHOTON_FF_LIB_DIR}//fonfra_inter2.o ${PHOTON_FF_LIB_DIR}//distribPert.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setI.o ${PHOTON_FF_LIB_DIR}//distribNonPert_setII.o  ${PHOTON_FF_LIB_DIR}//locate.o ${PHOTON_FF_LIB_DIR}//polin2.o ${PHOTON_FF_LIB_DIR}//polint.o utils.o functions.o expfunctions.o src/ghloprod.cpp -lgsl -lblas -lg2c -lm -lLHAPDF -L${LHAPDF_LIB_DIR}//lib -I${LHAPDF_LIB_DIR}//include/

expfunctions_optim.o : src/expfunctions.cpp
	g++ -c -Wall -O3 src/expfunctions.cpp -I${LHAPDF_LIB_DIR}//include/

functions_optim.o : src/functions.cpp
	g++ -c -Wall -O3 src/functions.cpp

utils_optim.o : src/utils.cpp
	g++ -c -Wall -O3 src/utils.cpp
