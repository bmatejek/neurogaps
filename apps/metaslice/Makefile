#
# Application name and list of source files.
#

NAME=metaslice
CCSRCS=$(NAME).cpp



#
# Dependency libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut -lX11
else
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lfglut -ljpeg -lpng
endif

#
# R3 application makefile
#

include ../../makefiles/Makefile.apps


