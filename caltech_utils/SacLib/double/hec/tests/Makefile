CC = /usr/bin/gcc
CFLAG =  
F90 = f90
OPTF = -O2


LINK_LIB_DIR =/opt/seismo-util/lib
LIB = ../../libDSacLib.a
all : testcorrf90 testcorrc
testcorrf90 : testcorrf90.f90 
	$(F90)   $(OPTF) -o testcorrf90 testcorrf90.f90 \
	         -L$(LINK_LIB_DIR) -lDRWFiles -lDSacio -lm $(LIB)

testcorrc: testcorrc.c  
	$(CC) -o testcorrc testcorrc.c -L$(LINK_LIB_DIR)  -lDRWFiles $(LIB) -lDSacio -lm -I$(LINK_LIB_DIR)/includes
                

clean : 
	\rm -f $(OBJS) *~  *.o testcorrc  testcorrf90

