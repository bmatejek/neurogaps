// Source file for the machine learning classifier class




////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////

#include "RNML.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(RNClassifier);


////////////////////////////////////////////////////////////////////
// Constructors/destructors 
////////////////////////////////////////////////////////////////////

RNClassifier::
RNClassifier(void)
{
}



RNClassifier::
~RNClassifier(void)
{
}



//////////////////////////////////////////////////////////////////
// Classification functions
//////////////////////////////////////////////////////////////////

unsigned int RNClassifier::
Classify(const std::vector<float>& attributes) const
{
   /* overridden */
   return 0.0;
}



float RNClassifier::
PredictBinaryProbability(const std::vector<float>& attributes) const
{
   /* overridden */
   return 0.0;
}



std::vector<unsigned int> RNClassifier::
Classify(RNDataset& testing_data) const
{
   /* overridden */
   return std::vector<unsigned int>();
}



std::vector<float> RNClassifier::
PredictBinaryProbabilities(RNDataset& testing_data) const
{
   /* overridden */
   return std::vector<float>();
}



//////////////////////////////////////////////////////////////////
// Input/output functions
//////////////////////////////////////////////////////////////////

unsigned int **RNClassifier::
PrintPerformance(RNDataset& dataset) const
{
   /* overridden */
   return NULL;
}



void RNClassifier::
LoadClassifier(const char input_filename[4096])
{
   /* overridden */
}



void RNClassifier::
SaveClassifier(const char output_filename[4096]) const
{
   /* overridden */
}
