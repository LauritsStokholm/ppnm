# Makefiles are not yet introduced for the course, but is however a great way
# to simplify and reproduce work. I have created this simple makefile just for
# easy replication, and to have a cleaner folder.


# First of all, makefile is clever, and it does not need every detail 
# A first try, would look like this:

# Recipe for objective file
helloworld: helloworld.c
	gcc helloworld.c -o helloworld
	./helloworld

# Recipe to clean up folder:
.PHONEY: clean
clean:
	rm helloworld
