#
# Application name and list of source files.
#

NAME=refinemeta
CCSRCS=$(NAME).cpp R3GridArray.cpp



#
# Dependency libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut -lX11
else
PKG_LIBS=-lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lfglut -ljpeg -lpng
endif

#
# R3 application makefile
#

include ../../makefiles/Makefile.apps


