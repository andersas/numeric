
src = $(wildcard *.c)
obj = $(src:.c=.o)

CC = gcc
CFLAGS = -Wall -O2 -mtune=native -march=native -mfpmath=sse -msse2 -msse3 -msse4.1 -pipe
LDFLAGS = -lm -L../linalg/ -Wl,-rpath,./ -llinalg

.PHONY: all
all:	main

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
	rm -f main *.out 

.PHONY: out
out:
	@echo With at human? Yuck, no thanks!


