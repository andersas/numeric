all: integration

integration:	integration.o	adaptive.o
	gcc -lm integration.o adaptive.o -o integration

integration.o:	integration.c adaptive.h
	gcc -Wall -c integration.c

adaptive.o:	adaptive.c adaptive.h
	gcc -Wall -c adaptive.c

.PHONY: clean
clean:
	rm -f *.o integration
