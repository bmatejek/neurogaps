#
# Application name 
#

NAME=dijkstraview



#
# Source files
#

CCSRCS=$(NAME).cpp

#
# Libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lNeuron -lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut -lX11
else
PKG_LIBS=-lNeuron -lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut
endif


#
# Include standard makefile
#

include ../../makefiles/Makefile.apps















