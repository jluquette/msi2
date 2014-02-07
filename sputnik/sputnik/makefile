#
# makefile for findrepeats
#

# for gnu compiler
#CC= gcc
#CFLAGS= -D_HPUX_SOURCE

# HP compiler
CC= cc
DEBUG= -g
CFLAGS=

PROG= sputnik
PROGOBJS= sputnik.o

all:	$(PROG)

$(PROG): sputnik.o
		 $(CC) -o $(PROG) $(PROGOBJS)

clean:	
		rm -f $(PROGOBJS)	

.c.o: 
		 $(CC) $(DEBUG) $(CFLAGS) -c $<

