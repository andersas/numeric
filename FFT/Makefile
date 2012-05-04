src = $(wildcard *.c)
obj = $(src:.c=.o)
libobj = $(filter-out main.o, $(obj))

CC = gcc
CFLAGS = -Wall -O2 -mtune=native -march=native -mfpmath=sse -msse2 -msse3 -msse4.1 -pipe -fpic
LDFLAGS = -lm

.PHONY: all
all:	main libraries

libraries: libmyfft.a libmyfft.so

libmyfft.a:	$(libobj)
	ar rcs libmyfft.a $(libobj)

libmyfft.so: libmyfft.so.1
	rm -f libmyfft.so
	ln -s libmyfft.so.1 libmyfft.so

libmyfft.so.1:	$(libobj)
	gcc -shared $(LDFLAGS) -W1,-soname,libmyfft.so.1 $(libobj) -o libmyfft.so.1

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
	rm -f main *.out *.a *.so*

.PHONY: out
out:
	@echo With at human? Yuck, no thanks!

.PHONY: examples
examples: main
	(time ./main fft > output-fft.txt) 2> time-fft.txt
	(time ./main dft > output-dft.txt) 2> time-dft.txt
	./main fft inverse > output-fft-inverse.txt
