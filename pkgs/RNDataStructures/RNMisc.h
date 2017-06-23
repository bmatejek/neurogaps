// Include file for miscellaneous useful functions

#ifndef __RN__MISC__H__
#define __RN__MISC__H__


#include "RNDataStructures.h"



// function list
void RNProgressBar(int index, int nindices);
void RNDeflateIntegerArray(int *entires, int nentries);
void RNBestFitLine(RNScalar *x, RNScalar *y, int n, RNScalar& alpha, RNScalar& beta, RNScalar& RSquared);

#endif
