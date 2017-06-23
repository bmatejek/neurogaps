// Include file for RNScalar distributions

#ifndef __RN_DISTRIBUTION_H__
#define __RN_DISTRIBUTION_H__

#include <vector>



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

class RNDistribution {
public:
   // constructors destructors
   RNDistribution(std::vector<RNScalar> &data);
   ~RNDistribution(void);

   // access functions
   RNScalar Kurtosis(void) const;
   RNScalar Maximum(void) const;
   RNScalar Mean(void) const;
   RNScalar Median(void) const;
   RNScalar Minimum(void) const;
   RNScalar Quintile(int dim) const;
   RNScalar Skew(void) const;
   RNScalar StdDev(void) const;
   int NPoints(void) const;


private:
   // instance variables
   std::vector<RNScalar> &data;
   RNScalar kurtosis;
   RNScalar maximum;
   RNScalar mean;
   RNScalar median;
   RNScalar minimum;
   RNScalar quintile[4];
   RNScalar skew;
   RNScalar stddev;
};



/* inline functions */

inline RNScalar RNDistribution::
Kurtosis(void) const
{
   // return the kurtosis
   return kurtosis;
}



inline RNScalar RNDistribution::
Maximum(void) const
{
   // return the maximum
   return maximum;
}



inline RNScalar RNDistribution::
Mean(void) const
{
   // return the mean
   return mean;
}



inline RNScalar RNDistribution::
Median(void) const
{
   // return the median 
   return median;
}



inline RNScalar RNDistribution::
Minimum(void) const
{
   // return the minimum
   return minimum;
}



inline RNScalar RNDistribution::
Quintile(int dim) const
{
   rn_assertion((1 <= dim) && (dim < 5));
   // return the dim quintile rank
   return quintile[dim - 1];
}



inline RNScalar RNDistribution::
StdDev(void) const
{
   // return the standard deviation
   return stddev;
}



inline RNScalar RNDistribution::
Skew(void) const
{
   // return the skew
   return skew;
}



inline int RNDistribution::
NPoints(void) const
{
   // return the number of points in the distribution
   return data.size();
}



#endif