/* Include file for GAPS data structures module */

#ifndef __RN__DATA_STRUCTURES__H__
#define __RN__DATA_STRUCTURES__H__

// define a custom assertion
#define rn_assertion(__x__) \
    if (!(__x__)) { fprintf(stderr, "Assertion error %s at line %d in file %s", #__x__, __LINE__, __FILE__); exit(-1); }


/* data structure include files */

#include "R3Graphics/R3Graphics.h"
#include "RNMinBinaryHeap.h"
#include "RNVector.h"
#include "RNMetaRawFile.h"
#include "RNHistogram.h"
#include "RNMisc.h"
#include "RNDistribution.h"
#include "RNSkeleton.h"
#include "RNPyPlot.h"



/* Initialization functions */

int RNInitBasics(void);
void RNStopBasics(void);


#endif







