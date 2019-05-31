# CDEBUGFLAGS = -g
CDEBUGFLAGS = -O2
CC   = cc

#----------------------------------------------------------*

plp:    plp.c plp.h nr.c nr.h
	${CC} ${CFLAGS} plp.c nr.c  -lm -o plp.ex

tlp:    tlp.c tplp.c tlp.h tplp.h nr.c nr.h
	${CC} ${CFLAGS} tlp.c tplp.c nr.c -lm -o tlp.ex

stm:    stm.c stm.h nr.c nr.h
	${CC} ${CFLAGS} stm.c nr.c -lm -o stm.ex

intnor: intnor.c intnor.h nr.c nr.h
	${CC} ${CFLAGS} intnor.c  nr.c -lm -o intnor.ex

clean:
	/bin/rm -f *.o

.KEEP_STATE:
