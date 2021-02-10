# Compiler
CC = gcc
CFLAGS += -Wall -std=gnu1x# $$(gsl-config --cflags)
LDLIBS += -lm #$$(gsl-config --libs)

# For big projects:
# VPATH is included for all files
#VPATH = ~~@~~
## vpath is file specific
#vpath .mca ..~//~


#==============================================================================#
# For gnuplot
#SCRIPT = log.gpi
#FORMAT = eps
#PLOTTER = gnuplot
#void =

#==============================================================================#
TARGETS = output.txt
# Remember: $^ (prerequisite) and $@ (target)

#.PHONEY: all
all: $(TARGETS)
#==============================================================================#

# Output of main
output.txt: main
	./$(word 1, $^)> $@
	cat $@

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	rm -f $(TARGETS)
	find ./ -type f -executable -delete
#	$(RM) *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps


