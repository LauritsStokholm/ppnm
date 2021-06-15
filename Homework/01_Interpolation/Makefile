CC = gcc
CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For gnuplot
SCRIPT = log.gpi
FORMAT = svg
IMAGE_VIEWER = ristretto
PLOTTER = gnuplot
void =

TARGETS = data.txt main lspline.$(FORMAT) qspline.$(FORMAT) cspline.$(FORMAT)


default: $(TARGETS)

# Figure 2
cspline.$(FORMAT): data.txt cspline.txt cgsl.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out left' >> $(SCRIPT)
	echo 'set title "Cubic spline"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" u 1:2 w p pt 2 ps 1 t "Data"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w l dt 1 lc "black" lw 5 t "cspline_{eval}"     \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:3 w l dt 1 lc "red"   lw 5 t "cspline_{deriv}"    \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:4 w l dt 1 lc "yellow"lw 5 t "cspline_{integ}"    \$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:2 w l dt 2 lc "grey"  lw 5 t "gsl cspline_{eval}" \$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:3 w l dt 2 lc "blue"  lw 5 t "gsl cspline_{deriv}"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:4 w l dt 2 lc "green" lw 5 t "gsl cspline_{integ}"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@

# Figure 2
qspline.$(FORMAT): data.txt qspline.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out left' >> $(SCRIPT)
	echo 'set title "Quadratic spline"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" u 1:2 w p pt 2 ps 1 t "Data"  \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w l dt 3 lc "black" lw 5 t "quadratic eval"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:3 w l dt 3 lc "red"   lw 5 t "quadratic deriv"\$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:4 w l dt 3 lc "blue"  lw 5 t "quadratic integ"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@



# Figure 1
lspline.$(FORMAT): data.txt lspline.txt lgsl.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out left' >> $(SCRIPT)
	echo 'set title "Linear spline"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"$(word 1, $^)" u 1:2 w p pt 2 ps 1 t "Data"  \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:2 w l dt 1 lc "black"  lw 5 t "linterp_{eval}"     \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:3 w l dt 1 lc "red"    lw 5 t "linterp_{deriv}"    \$(void)' >> $(SCRIPT)
	echo ', "$(word 2, $^)" u 1:4 w l dt 1 lc "yellow" lw 5 t "linterp_{integ}"    \$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:2 w l dt 2 lc "grey"   lw 5 t "gsl linterp_{eval}" \$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:3 w l dt 2 lc "blue"   lw 5 t "gsl linterp_{deriv}"\$(void)' >> $(SCRIPT)
	echo ', "$(word 3, $^)" u 1:4 w l dt 2 lc "green"  lw 5 t "gsl linterp_{integ}"\$(void)' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@


lspline.txt lgsl.txt: main data.txt
	./$<

main: main.c

data.txt:
	for x in `seq 0 1 10`; do echo "$$x\t$$(($$x * $$x * $$x + $$x * $$x - $$x - 2))" >> $@; done;



.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	$(RM) *.txt *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps


