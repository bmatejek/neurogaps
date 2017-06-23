// Include file for histogram

#ifndef __RN__HISTOGRAM__H__
#define __RN__HISTOGRAM__H__


#include "RNDataStructures.h"



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

class RNHistogram {
public:
   RNHistogram(void);
   ~RNHistogram(void);

   // manipulation functions
   void AddScalarQuantity(RNScalar value);
   void SortVector(void);

   // histogram functions
   int Histogram(int nbins, int *bins, RNScalar *boundaries, const char *output_filename = NULL);

private:
   RNVector<RNScalar> vector;
   RNBoolean sorted;
};

#endif