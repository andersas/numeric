
src = $(wildcard *.c)
obj = $(src:.c=.o)

CC = gcc
CFLAGS = -Wall
LDFLAGS = -lm -lpthread

.PHONY: all
all:	main

main:	$(obj)
	$(CC) $(LDFLAGS) $(obj) -o main

# ignore headers.
$(obj):	%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o main *.dat *.png


examples: main
	./main > out.txt
	gnuplot gnuplot
	
