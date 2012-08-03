CC = gcc
CFLAG =  
F90 = f90
OPTF = -O2

# options to link with c functions
PORT = -YEXT_NAMES=LCS  -s -B108 -YCFRL=1  # for Absoft f90 compiler
#PORT = -f -s -B108 -N90  # for Absoft f77 compiler

OBJS1 = cfft.o dcops.o hilenv.o xcorr.o
OBJS2 = testinc.o testinf90.o
OBJS = $(OBJS1) $(OBJS2) 

LIB_NAME = libHilEnv.a
LINK_LIB_DIR = /home/lqy/Lib
LINK_LIB_NAME = MySacio

testinf90 : testinf90.f90 
	$(F90)  $(PORT) $(OPTF) -o testinf90 testinf90.f90 -L$(LINK_LIB_DIR) \
		 -lMySacio -lHilEnv  -lm

testinc: testinc.c  
	$(CC) -o testinc testinc.c -L$(LINK_LIB_DIR) \
                -lMySacio -lHilEnv -lm

%.o : %.c
	$(CC) $(CFLAG) -c $*.c

lib : $(OBJS1)
	ar rvu $(LIB_NAME) $(OBJS1)

clean : 
	\rm -f $(OBJS) *~  testinc testinf90



