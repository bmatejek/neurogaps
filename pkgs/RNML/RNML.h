// Include file for the machine learning class

#ifndef __RNML__H__
#define __RNML__H__



// external dependencies

#include <vector>
#include <string>
#include <cv.h>
#include <ml.h>



//////////////////////////////////////////////////////////////////////
//// PREDICTION FUNCTIONS
//////////////////////////////////////////////////////////////////////



// dependency includes

#include "R3Shapes/R3Shapes.h"
#include "RNDataStructures/RNDataStructures.h"


// class declarations

class RNDataset;
class RNClassifier;
class RNRandomForest;



// rnml includes

#include "RNClassifier.h"
#include "RNRandomForest.h"
#include "RNDataset.h"



// initialization functions

int RNMLInit(void);
void RNMLStop(void);



// initialization functions

int RNMLInit(void);
void RNMLStop(void);



#endif
