CC = gcc
CFLAGS += -Wall $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

FORMAT = svg
SCRIPT1 = log_erf.gpi
SCRIPT2 = log_bessel.gpi
PLOTTER = gnuplot
IMAGE_VIEW = ristretto

default: task_a.txt erf.$(FORMAT) bessel.$(FORMAT)


#besselfunction
bessel.$(FORMAT): j0.txt j1.txt j2.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT2)
	echo 'set output "$@"' >> $(SCRIPT2)
	echo 'set key out' >> $(SCRIPT2)
	echo 'set title "Implementation of bessel(x) comparison"' >> $(SCRIPT2)
	echo 'set grid' >> $(SCRIPT2)
	echo 'set xlabel "x"' >> $(SCRIPT2)
	echo 'set ylabel "y"' >> $(SCRIPT2)
	echo 'set xtics 5' >> $(SCRIPT2)
	echo 'set ytics 0.1' >> $(SCRIPT2)
	echo 'set xrange [0:$(MAX)]' >> $(SCRIPT2)
	echo 'array point[1]' >> $(SCRIPT2)
	echo 'plot \$(void)' >> $(SCRIPT2)
#	echo '"$(word 1, $^)" using 1:2 with points title "Data" \$(void)' >> $(SCRIPT2)
#	echo ', "$(word 1, $^)" using 1:3 with lines title "fit"' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"res0" with lines linetype 1 dashtype 1 linewidth 5 linecolor "red"    title "j0 (qags)",          \$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"res1" with lines linetype 1 dashtype 2 linewidth 5 linecolor "grey"   title "j0 (cquad)",         \$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"res2" with lines linetype 1 dashtype 3 linewidth 5 linecolor "black"  title "j0 (romberg)",       \$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"res3" with lines linetype 1 dashtype 4 linewidth 5 linecolor "yellow" title "j0 (gauss-legendre)",\$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"func"  with lines linetype 1 dashtype 5 linewidth 5 linecolor "green" title "j0(x)",\$(void)' >> $(SCRIPT2)
	echo '"$(word 2, $^)" using "x":"res0" with lines linetype 2 dashtype 1 linewidth 5 linecolor "red"    title "j1 (qags)",          \$(void)' >> $(SCRIPT2)
	echo '"$(word 2, $^)" using "x":"res1" with lines linetype 2 dashtype 2 linewidth 5 linecolor "grey"   title "j1 (cquad)",         \$(void)' >> $(SCRIPT2)
	echo '"$(word 2, $^)" using "x":"res2" with lines linetype 2 dashtype 3 linewidth 5 linecolor "black"  title "j1 (romberg)",       \$(void)' >> $(SCRIPT2)
	echo '"$(word 2, $^)" using "x":"res3" with lines linetype 2 dashtype 4 linewidth 5 linecolor "yellow" title "j1 (gauss-legendre)",\$(void)' >> $(SCRIPT2)
	echo '"$(word 2, $^)" using "x":"func"  with lines linetype 2 dashtype 5 linewidth 5 linecolor "green" title "j1(x)",\$(void)' >> $(SCRIPT2)
	echo '"$(word 3, $^)" using "x":"res0" with lines linetype 3 dashtype 1 linewidth 5 linecolor "red"    title "j2 (qags)",          \$(void)' >> $(SCRIPT2)
	echo '"$(word 3, $^)" using "x":"res1" with lines linetype 3 dashtype 2 linewidth 5 linecolor "grey"   title "j2 (cquad)",         \$(void)' >> $(SCRIPT2)
	echo '"$(word 3, $^)" using "x":"res2" with lines linetype 3 dashtype 3 linewidth 5 linecolor "black"  title "j2 (romberg)",       \$(void)' >> $(SCRIPT2)
	echo '"$(word 3, $^)" using "x":"res3" with lines linetype 3 dashtype 4 linewidth 5 linecolor "yellow" title "j2 (gauss-legendre)",\$(void)' >> $(SCRIPT2)
	echo '"$(word 3, $^)" using "x":"func"  with lines linetype 3 dashtype 5 linewidth 5 linecolor "green" title "j2(x)",\$(void)' >> $(SCRIPT2)
	$(PLOTTER) $(SCRIPT2)
	$(IMAGE_VIEW) $@


# errorfunction
erf.$(FORMAT): task_b.txt
	echo 'set terminal $(FORMAT)' > $(SCRIPT1)
	echo 'set output "$@"' >> $(SCRIPT1)
	echo 'set key out' >> $(SCRIPT1)
	echo 'set title "Implementation of erf(x) comparison"' >> $(SCRIPT1)
	echo 'set grid' >> $(SCRIPT1)
	echo 'set xlabel "x"' >> $(SCRIPT1)
	echo 'set ylabel "y"' >> $(SCRIPT1)
	echo 'set xtics 1' >> $(SCRIPT1)
	echo 'set ytics 1' >> $(SCRIPT1)
	echo 'set xrange [0:$(MAX)]' >> $(SCRIPT1)
	echo 'array point[1]' >> $(SCRIPT1)
	echo 'plot \$(void)' >> $(SCRIPT1)
#	echo '"$(word 1, $^)" using 1:2 with points title "Data" \$(void)' >> $(SCRIPT1)
#	echo ', "$(word 1, $^)" using 1:3 with lines title "fit"' >> $(SCRIPT1)
	echo '"$<" using "x":"res0" with lines linetype 1 dashtype 1 linewidth 5 linecolor "red"    title "qags",          \$(void)' >> $(SCRIPT1)
	echo '"$<" using "x":"res1" with lines linetype 1 dashtype 2 linewidth 5 linecolor "grey"   title "cquad",         \$(void)' >> $(SCRIPT1)
	echo '"$<" using "x":"res2" with lines linetype 1 dashtype 3 linewidth 5 linecolor "black"  title "romberg",       \$(void)' >> $(SCRIPT1)
	echo '"$<" using "x":"res3" with lines linetype 1 dashtype 4 linewidth 5 linecolor "yellow" title "gauss-legendre",\$(void)' >> $(SCRIPT1)
	echo '"$<" using "x":"func"  with lines linetype 1 dashtype 5 linewidth 5 linecolor "green" title "erf(x)",\$(void)' >> $(SCRIPT1)
	$(PLOTTER) $(SCRIPT1)
	$(IMAGE_VIEW) $@



task_a.txt: task_a
	./$<
	cat $@

task_b.txt: task_b
	./$<

j0.txt j1.txt j2.txt: task_c
	./$<

task_a: task_a.c
task_b: task_b.c
task_c: task_c.c


.PHONEY: clean
clean:
	$(RM) *.txt *.$(FORMAT) *.gpi
	find ./ -type f -executable -delete

