// Source file for RNScalar distributions



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RNDataStructures/RNDataStructures.h"
#include <algorithm>



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

RNDistribution::
RNDistribution(std::vector<RNScalar> &data) :
data(data)
{
   unsigned int nentries = data.size();
   if (nentries > 0) {
      sort(data.begin(), data.end());
      // update all of the distribution variables
      maximum = data[nentries - 1];
      minimum = data[0];
      if (nentries % 2 == 0)
         median = 0.5 * (data[nentries / 2 - 1] + data[nentries / 2]);
      else
         median = data[nentries / 2];
      for (int dim = 1; dim < 5; ++dim) {
         quintile[dim - 1] = data[(int)((RNScalar)(nentries / 5) * dim)];
      }

      // calculate the mean
      mean = 0.0;
      for (unsigned int ie = 0; ie < nentries; ++ie)
         mean += data[ie];
      mean /= nentries;

      if (nentries == 1) {
         stddev = 0.0;
         skew = 0.0;
         kurtosis = 0.0;
      }
      else {
         stddev = 0.0;
         skew = 0.0;
         kurtosis = 0.0;
         for (unsigned int ie = 0; ie < nentries; ++ie) {
            stddev += (data[ie] - mean) * (data[ie] - mean);
            skew += (data[ie] - mean) * (data[ie] - mean) * (data[ie] - mean);
            kurtosis += (data[ie] - mean) * (data[ie] - mean) * (data[ie] - mean) * (data[ie] - mean);
         }
         stddev /= (nentries - 1);
         stddev = sqrt(stddev);
         if (stddev > 0.0) {
            skew /= nentries;
            skew /= (stddev * stddev * stddev);
            kurtosis /= nentries;
            kurtosis /= (stddev * stddev * stddev * stddev);
         }
         else {
            skew = 0.0;
            kurtosis = 0.0;
         }
      }
   }
   else {
      kurtosis = 0.0;
      maximum = 0.0;
      mean = 0.0;
      median = 0.0;
      minimum = 0.0;
      for (int dim = 0; dim < 4; ++dim)
         quintile[dim] = 0.0; 
      stddev = 0.0;
      skew = 0.0;
   }
}



RNDistribution::
~RNDistribution(void)
{
}



