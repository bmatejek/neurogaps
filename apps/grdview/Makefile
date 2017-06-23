#
# Application name 
#

NAME=grdview



#
# Source files.
#

CCSRCS=$(NAME).cpp



#
# Libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut -lX11
else
PKG_LIBS=-lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lfglut -ljpeg -lpng
endif



#
# Include standard makefile
#

include ../../makefiles/Makefile.apps

