# Compiler
CC = gcc
CFLAGS += -O -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For big projects:
# VPATH is included for all files
#VPATH = ~~@~~
## vpath is file specific
#vpath .mca ..~//~


#==============================================================================#
# For gnuplot
SCRIPT1 = log_erf.gpi
SCRIPT2 = log_my_tgamma.gpi
SCRIPT3 = log_my_lgamma.gpi
SCRIPT4 = log_my_cgamma.gpi
FORMAT = eps
PLOTTER = gnuplot
void =

#==============================================================================#
TARGETS = erf.eps my_tgamma.eps my_lgamma.eps my_cgamma.eps
# Remember: $^ (prerequisite) and $@ (target)
#==============================================================================#

.PHONEY: all
all: $(TARGETS)

#==============================================================================#
# COMPLEX GAMMA
my_cgamma.$(FORMAT): my_cgamma_out.txt Makefile my_cgamma_in.txt
	echo 'set terminal $(FORMAT)' > $(SCRIPT4)
	echo 'set output "$@"' >> $(SCRIPT4)
# Image format
	echo 'set title "Absolute value of complex Gamma function"' >> $(SCRIPT4)
	echo 'set nokey; set grid' >> $(SCRIPT4)
	echo 'set xlabel "Re(z)"; set ylabel "Im(z)"; set zlabel "|G(z)|"' >> $(SCRIPT4)
	echo 'set xrange [-4.5:4.5]; set yrange [-1.5:1.5]; set zrange [-0.5:6]' >> $(SCRIPT4)
	echo 'set xtics 1; set ytics 1; set ztics 1 ' >> $(SCRIPT4)
	echo 'set view 75, 350, 1, 1' >> $(SCRIPT4)
# Contour settings
# hidden3d: cant be seen through - removes lines behind other
	#echo 'set hidden3d' >> $(SCRIPT4)
	#echo 'set dgrid3d spline' >> $(SCRIPT4)
#	echo 'set contour surface' >> $(SCRIPT4)
#	echo 'set cntrparam levels incr -0.5, 0.25, 5' >> $(SCRIPT4)
	echo 'splot "$<" u "x":"y":"z" w l lc "black"' >> $(SCRIPT4)
	$(PLOTTER) $(SCRIPT4)
	gv $@
#==============================================================================#
# LOGGAMMAFUNCTION
my_lgamma.$(FORMAT): my_lgamma_out.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT3)
	echo 'set output "$@"' >> $(SCRIPT3)
	echo 'set key out' >> $(SCRIPT3)
	echo 'set title "Comparison of Logarithmic Gammafunctions"' >> $(SCRIPT3)
	echo 'set grid' >> $(SCRIPT3)
	echo 'set xlabel "x"' >> $(SCRIPT3)
	echo 'set ylabel "y"' >> $(SCRIPT3)
	echo 'set xtics 1' >> $(SCRIPT3)
	echo 'set ytics 1' >> $(SCRIPT3)
#	echo 'set yrange [-7:7]' >> $(SCRIPT3)
	echo 'set xrange [-1:5]' >> $(SCRIPT3)
	echo 'plot \$(void)' >> $(SCRIPT3)
	echo '"$(word 1, $^)" using "x":"math" with lines linetype 1 dashtype 1 linewidth 5 linecolor "red"   title "math.h",        \$(void)' >> $(SCRIPT3)
	echo '"$(word 1, $^)" using "x":"gsl"  with lines linetype 1 dashtype 3 linewidth 5 linecolor "grey"  title "gsl",           \$(void)' >> $(SCRIPT3)
	echo '"$(word 1, $^)" using "x":"my"   with lines linetype 1 dashtype 3 linewidth 5 linecolor "black" title "self-defined",  \$(void)' >> $(SCRIPT3)
	$(PLOTTER) $(SCRIPT3)
	gv $@

#==============================================================================#
# GAMMAFUNCTION
my_tgamma.$(FORMAT): my_tgamma_out.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT2)
	echo 'set output "$@"' >> $(SCRIPT2)
	echo 'set key out' >> $(SCRIPT2)
	echo 'set title "Comparison of Gammafunctions"' >> $(SCRIPT2)
	echo 'set grid' >> $(SCRIPT2)
	echo 'set xlabel "x"' >> $(SCRIPT2)
	echo 'set ylabel "y"' >> $(SCRIPT2)
	echo 'set xtics 1' >> $(SCRIPT2)
	echo 'set ytics 1' >> $(SCRIPT2)
	echo 'set yrange [-7:7]' >> $(SCRIPT2)
	echo 'set xrange [-5:5]' >> $(SCRIPT2)
	echo 'array point[1]' >> $(SCRIPT2)
	echo 'plot \$(void)' >> $(SCRIPT2)
#	echo '"$(word 1, $^)" using 1:2 with points title "Data" \$(void)' >> $(SCRIPT2)
#	echo ', "$(word 1, $^)" using 1:3 with lines title "fit"' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"math" with lines linetype 1 dashtype 1 linewidth 5 linecolor "red"   title "math.h",        \$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"gsl"  with lines linetype 1 dashtype 2 linewidth 5 linecolor "grey"  title "gsl",           \$(void)' >> $(SCRIPT2)
	echo '"$(word 1, $^)" using "x":"my"   with lines linetype 1 dashtype 3 linewidth 5 linecolor "black" title "self-defined",  \$(void)' >> $(SCRIPT2)
	echo 'point using (1):(1) pt 2 ps 1 lc "black" t "n!", \$(void)' >> $(SCRIPT2)
	echo 'point using (2):(1) pt 2 ps 1 lc "black" notitle, \$(void)' >> $(SCRIPT2)
	echo 'point using (3):(2) pt 2 ps 1 lc "black" notitle, \$(void)' >> $(SCRIPT2)
	echo 'point using (4):(6) pt 2 ps 1 lc "black" notitle, \$(void)' >> $(SCRIPT2)
	$(PLOTTER) $(SCRIPT2)
	gv $@
#==============================================================================#

my_tgamma_out.txt my_lgamma_out.txt my_cgamma_out.txt: gamma my_tgamma_in.txt my_lgamma_in.txt my_cgamma_in.txt
	./$<

gamma: gamma.c

my_cgamma_in.txt:
	for i in $$(seq -4.2 0.025 4.2);\
		do\
		for j in $$(seq -2 0.025 2);\
			do\
			echo "$$i $$j" >> $@;\
			done;\
		done

my_lgamma_in.txt:
	seq 0.0001 0.001 5 >> $@

my_tgamma_in.txt:
	seq -5.0001 0.001 -4.1 > $@
	seq -4.0001 0.001 -3.01 >> $@
	seq -3.001  0.01 -2.2 >> $@
	seq -2.001  0.01 -1.1 >> $@
	seq -1.001  0.1 -0.01 >> $@
	seq  0.001  0.1 5 >> $@


#==============================================================================#
# ERRORFUNCTION
# Figure 1
erf.$(FORMAT): erf_output.txt Makefile
	echo 'set terminal $(FORMAT)' > $(SCRIPT1)
	echo 'set output "$@"' >> $(SCRIPT1)
	echo 'set key out' >> $(SCRIPT1)
	echo 'set title "Comparison of Errorfunctions"' >> $(SCRIPT1)
	echo 'set grid' >> $(SCRIPT1)
	echo 'set xlabel "x"' >> $(SCRIPT1)
	echo 'set ylabel "y"' >> $(SCRIPT1)
	echo 'set xrange [0:2]' >> $(SCRIPT1)
	echo 'set yrange [0:1.5]' >> $(SCRIPT1)
	echo 'plot \$(void)' >> $(SCRIPT1)
#	echo '"$(word 1, $^)" using 1:2 with points title "Data" \$(void)' >> $(SCRIPT1)
#	echo ', "$(word 1, $^)" using 1:3 with lines title "fit"' >> $(SCRIPT1)
	echo '"$(word 1, $^)" using "x":"matherf" with lines linetype 1 dashtype 1 linewidth 5 linecolor "red"   title "math.h",\$(void)' >> $(SCRIPT1)
	echo '"$(word 1, $^)" using "x":"gslerf"  with lines linetype 1 dashtype 2 linewidth 5 linecolor "grey"  title "gsl", \$(void)' >> $(SCRIPT1)
	echo '"$(word 1, $^)" using "x":"myerf"   with lines linetype 1 dashtype 3 linewidth 5 linecolor "black" title "self-defined"'             >> $(SCRIPT1)
	$(PLOTTER) $(SCRIPT1)
	gv $@

# Output of main
erf_output.txt: erf erf_input.txt
	./$(word 1, $^) < $(word 2, $^)

erf: erf.c

erf_input.txt:
	seq 0 0.01 4 > $@
#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
	$(RM) $(SCRIPT1) $(SCRIPT2) *.dat *.svg *.png *.pdf log* *.eps *.txt


