.PHONY: all
all: interpolation lspline_interpolation qspline_interpolation

# See the "tests" target below

qspline_interpolation: qspline_interpolation.o qspline.o
	gcc qspline_interpolation.o qspline.o -o qspline_interpolation

qspline_interpolation.o: qspline_interpolation.c qspline.h data.h
	gcc -c qspline_interpolation.c

qspline.o: qspline.c qspline.h data.h
	gcc -c qspline.c




lspline_interpolation: lspline_interpolation.o lspline.o
	gcc lspline_interpolation.o lspline.o -o lspline_interpolation

lspline_interpolation.o: lspline_interpolation.c lspline.h data.h
	gcc -c lspline_interpolation.c

lspline.o: lspline.c lspline.h data.h
	gcc -c lspline.c

interpolation: interpolation.o lagrange.o
	gcc interpolation.o lagrange.o -o interpolation

interpolation.o: interpolation.c lagrange.h data.h
	gcc -c interpolation.c

lagrange.o: lagrange.c lagrange.h data.h
	gcc -c lagrange.c -o lagrange.o


tests: interpolation lspline_interpolation qspline_interpolation
	@ for i in $$(seq 1 9) ; do ./testcase$$i.py > $$i.dat ; done
	@ for i in $$(seq 1 9) ; do ./interpolation -20 20 0.023 < $$i.dat > $$i.res.dat ; done
	@ for i in $$(seq 1 9) ; do ./lspline_interpolation -20 20 0.023 < $$i.dat > $$i.res.lspline.dat ; done
	@ for i in $$(seq 1 9) ; do ./qspline_interpolation -20 20 0.023 < $$i.dat > $$i.res.qspline.dat ; done
	
	@ for i in $$(seq 1 9) ; do \
		echo set term png > gnuplot.script ; \
		echo set yrange [-0.5:2.1] >> gnuplot.script ; \
		echo set output \'$$i-lagrange.png\' >> gnuplot.script ; \
		echo plot \'$$i.dat\' title \'Tabulated\' pointsize 2, \'$$i.res.dat\' title \'Lagrange\' pointsize 0.1 >> gnuplot.script ; \
		echo set output \'$$i-lspline.png\' >> gnuplot.script ; \
		echo plot \'$$i.dat\' title \'Tabulated\' pointsize 2, \'$$i.res.lspline.dat\' title \'Linear spline\' pointsize 0.1 >> gnuplot.script ; \
		echo set output \'$$i-qspline.png\' >> gnuplot.script ; \
		echo plot \'$$i.dat\' title \'Tabulated\' pointsize 2, \'$$i.res.qspline.dat\' title \'Quadratic spline\' pointsize 0.1 >> gnuplot.script ; \
		gnuplot < gnuplot.script \
	  ; done
	@ rm -f gnuplot.script
	


.PHONY: clean
clean:
	rm -f *.o interpolation *.dat *.png

