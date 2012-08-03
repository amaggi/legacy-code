
default:sac2mt5

FLAGS=-g
CFLAGS=-c
OBJECTS=read_sac.o resp2pz.o which_pz.o write_dsn.o  sac2mt5.o

sac2mt5:${OBJECTS}
	gcc ${FLAGS} ${OBJECTS} -o sac2mt5

read_sac.o: read_sac.c
	gcc ${FLAGS} ${CFLAGS} read_sac.c

resp2pz.o: resp2pz.c
	gcc ${FLAGS} ${CFLAGS} resp2pz.c

which_pz.o: which_pz.c
	gcc ${FLAGS} ${CFLAGS} which_pz.c

write_dsn.o: write_dsn.c
	gcc ${FLAGS} ${CFLAGS} write_dsn.c

sac2mt5.o: sac2mt5.c
	gcc ${FLAGS} ${CFLAGS} sac2mt5.c

clean:
	rm *.o
