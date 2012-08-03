CC = /usr/bin/gcc
CFLAG =  
F90 = f90
OPTF = -O2
F90FLAGS = -YEXT_NAMES=LCS -YEXT_SFX=_

# no optimization for fortran code while compilation,
# otherwise some bad things might occur!!!

# options to link with c functions -now in $(F90FLAGS)

OBJS1 = cfft.o dcops.o hilenv.o corr1.o corr2.o convolve.o
LIB_NAME = libHilEnv.a
LINK_LIB_DIR = /opt/seismo-util/lib


testcorrf90 : testcorrf90.o 
	$(F90)  $(F90FLAGS) $(OPTF) -o testcorrf90 testcorrf90.o \
	         -L$(LINK_LIB_DIR)  -lHilEnv  -lRWFiles -lMySacio -lm

testcorrc: testcorrc.c  
	$(CC) -o testcorrc testcorrc.c -L$(LINK_LIB_DIR) \
                -lMySacio -lHilEnv -lm

testhilenvf90 : testhilenvf90.o
	$(F90)  $(F90FLAGS) $(OPTF) -o testhilenvf90 testhilenvf90.o \
                 -L$(LINK_LIB_DIR) -lMySacio -lHilEnv  -lm

testhilenvc: testhilenvc.c
	$(CC) -o testhilenvc testhilenvc.c -L$(LINK_LIB_DIR) \
                -lMySacio -lHilEnv -lm

testcorrf90.o: testcorrf90.f90
	$(F90) $(F90FLAGS)  -c testcorrf90.o testcorrf90.f90

testhilenvf90.o: testhilenvf90.f90
	$(F90) $(F90FLAGS)  -c testhilenvf90.o testhilenvf90.f90


%.o : %.c
	$(CC) $(CFLAG) -c $*.c

lib : $(OBJS1)
	ar rvu $(LIB_NAME) $(OBJS1)

clean : 
	\rm -f $(OBJS) *~  *.o testcorrc  testcorrf90

cleanall : 
	\rm -f $(OBJS) *~  *.o testcorrc testcorrf90 *.a

