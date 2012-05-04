src = $(wildcard *.c)
obj = $(src:.c=.o)
compobj = $(filter-out decompress.o, $(obj))
decompobj = $(filter-out compress.o, $(obj))


CC = gcc
CFLAGS = -Wall -O2 -mtune=native -march=native -mfpmath=sse -msse2 -msse3 -msse4.1 -pipe
LDFLAGS = -lm -L./ -Wl,-rpath,./ -lmyfft

.PHONY: all
all:	compress decompress

compress:	$(compobj)
	$(CC) $(LDFLAGS) $(compobj) -o compress

decompress:	$(decompobj)
	$(CC) $(LDFLAGS) $(decompobj) -o decompress

# ignore headers.
$(obj):	%.o : %.c compress.h
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean
clean:
	rm -f *.o
.PHONY: distclean
distclean:	clean
	rm -f compress decompress

.PHONY: out
out:
	@echo With at human? Yuck, no thanks!
