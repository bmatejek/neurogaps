#
# Application name 
#

NAME=graphcut



#
# Source files
#

CCSRCS=$(NAME).cpp

#
# Libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lgco -lX11 -lNeuron -lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut
else
PKG_LIBS=-lgco -lNeuron -lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut
endif


#
# Include standard makefile
#

include ../../makefiles/Makefile.apps















