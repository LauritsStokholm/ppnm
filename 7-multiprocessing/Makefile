# Compiler
CC = gcc
CFLAGS += -pthread
LDLIBS += -lm -pthread

SCRIPT = log_convergence.gpi
FORMAT = svg
PLOTTER = gnuplot
IMAGE_VIEWER = ristretto


convergence.$(FORMAT): convergence_input.txt convergence_input_omp.txt Makefile
	echo '\
	set terminal "$(FORMAT)";\
	set output "$@";\
	set title "Monte Carlo calculation of pi";\
	set logscale x; set logscale y;\
	set xlabel "Number of points; N"; set ylabel "Absolute error";\
	a = 1;\
	b = 1;\
	f(x) = a / sqrt(x);\
	g(x) = b / sqrt(x);\
	fit f(x) "$<" using 1:2 via a;\
	fit g(x) "$(word 2, $^)" using 1:2 via b;\
	plot "$<" u 1:2 w lp lt 1 lc "black"\
	, "$(word 2, $^)" u 1:2 w lp lt 2 lc "black"\
	, f(x) with lines title "$$1/\sqrt{N}$$"\
	, g(x) with lines title "$$1/\sqrt{N}$$ (OMP)"\
	' > $(SCRIPT)
	$(PLOTTER) $(SCRIPT)
	$(IMAGE_VIEWER) $@


#==============================================================================#

ns = 1 1.5 2 4 8
convergence_input.txt: main
	./$< -n 100000 > $@
	for i in $$(seq 3e5 2e5 1e6);\
		do\
		./$< -n $$i >> $@;\
		done
	for i in $$(seq 1e6 2e6 1e7);\
		do\
		./$< -n $$i >> $@;\
		done
	for i in $$(seq 1e7 2e7 1e8);\
		do\
		./$< -n $$i >> $@;\
		done

convergence_input_omp.txt: omp
	./$< -n 100000 > $@
	for i in $$(seq 3e5 2e5 1e6);\
		do\
		./$< -n $$i >> $@;\
		done
	for i in $$(seq 1e6 2e6 1e7);\
		do\
		./$< -n $$i >> $@;\
		done
	for i in $$(seq 1e7 2e7 1e8);\
		do\
		./$< -n $$i >> $@;\
		done


main: main.c

omp: omp.c


#==============================================================================#

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	$(RM) *.txt *.log *.gpi *.$(FORMAT)


