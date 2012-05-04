
src = $(wildcard *.c)
obj = $(src:.c=.o)
libobj = $(filter-out main.o, $(obj))

CC = gcc
CFLAGS = -Wall -O2 -mtune=native -march=native -mfpmath=sse -msse2 -msse3 -msse4.1 -pipe -fpic
LDFLAGS = -lm -lpthread -lgsl -lgslcblas

.PHONY: all
all:	main libraries

libraries: liblinalg.a	liblinalg.so

liblinalg.a:	$(libobj)
	ar rcs liblinalg.a $(libobj)

liblinalg.so.1:	$(libobj)
	gcc -shared $(LDFLAGS) -W1,-soname,liblinalg.so.1 $(libobj) -o liblinalg.so.1
liblinalg.so:	liblinalg.so.1
	rm -f liblinalg.so
	ln -s liblinalg.so.1 liblinalg.so
	

main:	$(obj)
	$(CC) $(LDFLAGS) $(obj) -o main

# ignore headers.
$(obj):	%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean
clean:
	rm -f *.o
.PHONY: distclean
distclean:	clean
	rm -f main *.out *.png *.a *.so* *.dat

tests:
	./main > output.txt
	@echo Testing QR decomposition of 3000x3000 matrix with
	@echo gsl and my parallel QR. Output in time.*
	(time ./main gsl) 2> time.gslQR
	(time ./main my) 2> time.myQR
