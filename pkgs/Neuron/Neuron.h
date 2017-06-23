// Include file for Neuron class

#ifndef __NEURON__H__
#define __NEURON__H__



// dependency includes

#include "R3Shapes/R3Shapes.h"
#include "RNDataStructures/RNDataStructures.h"



// class declarations

class NeuronVolume;
class NeuronReconstruction;
class NeuronHumanLabel;
class NeuronPrediction;
class NeuronSupervoxel;
class NeuronCellular;
class NeuronExtracellular;
class NeuronPredictionBoundary;
class NeuronBoundary;
class NeuronVoxel;
class NeuronData;



// neuron includes

#include "NeuronData.h"
#include "NeuronPredictionBoundary.h"
#include "NeuronBoundary.h"
#include "NeuronVoxel.h"
#include "NeuronVolume.h"
#include "NeuronReconstruction.h"
#include "NeuronHumanLabel.h"
#include "NeuronPrediction.h"
#include "NeuronSupervoxel.h"
#include "NeuronCellular.h"
#include "NeuronExtracellular.h"



// initialization functions

int NeuronInit(void);
void NeuronStop(void);



// debugging constants

#define NEURON_DEBUG TRUE



// useful enums for binary smoothing

enum { MEAN_BOUNDARY_METRIC, MIN_BOUNDARY_METRIC, RANDOM_WALK_METRIC, TRUTH_METRIC };



#endif