// Include file for the neuron prediction boundary class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronPredictionBoundary {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronPredictionBoundary(void);
   virtual ~NeuronPredictionBoundary(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // data access functions
   NeuronData *Data(void) const;
   int DataIndex(void) const;

   // boundary access functions
   int NAffinities(void) const;
   RNScalar Affinity(int affinity_index) const;

   // index function
   int PredictionOneIndex(void) const;
   int PredictionTwoIndex(void) const;

   // return the volumes
   NeuronPrediction *PredictionOne(void) const;
   NeuronPrediction *PredictionTwo(void) const;
   NeuronPrediction *OtherPrediction(const NeuronPrediction *prediction) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // return statistics
   RNScalar Maximum(void) const;
   RNScalar Mean(void) const;
   RNScalar Median(void) const;
   RNScalar Minimum(void) const;
   RNScalar StdDev(void) const;
   RNScalar Skew(void) const;
   RNScalar BoundaryRank(void) const;
   RNScalar BoundaryScaledRanking(void) const;
   RNScalar CombinedVoxels(void) const;
   RNScalar CombinedVoxelsProportion(void) const;


   // sophisticated features
   RNScalar RandomWalk(RNBoolean drunkards_walk) const;

   // get user data
   void *UserData(void) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // only used by NeuronData
   void CreateNAffinities(int naffinities);
   void UpdateVoxelAffinity(RNScalar affinity, int index);
   void UpdateStatistics(void);

   // set user data
   void SetUserData(void *user_data);


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////

public:
   // I/O functions
   RNBoolean AreAffinitiesResident(void) const;
   int ReadAffinities(void);
   int ReadAffinities(FILE *fp, RNBoolean seek = TRUE);
   int WriteAffinities(FILE *fp, RNBoolean seek = TRUE);
   int ReleaseAffinities(void);

   RNBoolean AreVoxelsResident(void) const;
   int ReadVoxels(void);
   int ReadVoxels(FILE *fp, RNBoolean seek = TRUE);
   int ReleaseVoxels(void);

   // file functions
   unsigned long long FileOffset(void) const;
   void SetFileOffset(unsigned long long file_offset);
   void SetNAffinities(int naffinities);

protected:
   friend class NeuronData;
   friend class NeuronPrediction;
   NeuronData *data;
   int data_index;
   int naffinities;
   unsigned long long file_offset;
   unsigned int read_voxel_count;
   unsigned int read_affinities_count;
   RNScalar *affinities;
   NeuronPrediction *prediction_one;
   NeuronPrediction *prediction_two;
   int prediction_one_index;
   int prediction_two_index;
   RNScalar maximum;
   RNScalar mean;
   RNScalar median;
   RNScalar minimum;
   RNScalar stddev;
   RNScalar skew;
   void *user_data;
};


/////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/* access functions */

inline NeuronData *NeuronPredictionBoundary::
Data(void) const
{
   // return neuron data
   return data;
}



inline int NeuronPredictionBoundary::
DataIndex(void) const
{
   // return the index of this boundary
   return data_index;
}



inline int NeuronPredictionBoundary::
NAffinities(void) const
{
   // return number of affinities
   return naffinities;
}



inline RNScalar NeuronPredictionBoundary::
Affinity(int affinity_index) const
{
   rn_assertion((0 <= affinity_index) && (affinity_index < naffinities));
   rn_assertion(affinities != NULL);
   // return kth affinity
   return affinities[affinity_index];
}



inline int NeuronPredictionBoundary::
PredictionOneIndex(void) const
{
   // return the index of the first volume
   return prediction_one_index;
}



inline int NeuronPredictionBoundary::
PredictionTwoIndex(void) const
{
   // return the index of the second volume
   return prediction_two_index;
}



inline NeuronPrediction *NeuronPredictionBoundary::
PredictionOne(void) const
{
   // return the first volume
   return prediction_one;
}



inline NeuronPrediction *NeuronPredictionBoundary::
PredictionTwo(void) const
{
   // return the second volume
   return prediction_two;
}



inline NeuronPrediction *NeuronPredictionBoundary::
OtherPrediction(const NeuronPrediction *prediction) const
{
   rn_assertion(prediction != NULL);
   // return other volume
   if (prediction_one == prediction) return prediction_two;
   else if (prediction_two == prediction) return prediction_one;
   else return NULL;
}



/* property functions */

inline RNScalar NeuronPredictionBoundary::
Maximum(void) const
{
   if (maximum == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return the maximum
   return maximum;
}



inline RNScalar NeuronPredictionBoundary::
Mean(void) const
{
   if (mean == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return the average
   return mean;
}



inline RNScalar NeuronPredictionBoundary::
Median(void) const
{
   if (median == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return the median
   return median;
}



inline RNScalar NeuronPredictionBoundary::
Minimum(void) const
{
   if (minimum == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return the minimum
   return minimum;
}



inline RNScalar NeuronPredictionBoundary::
StdDev(void) const
{
   if (stddev == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return the standard deviation
   return stddev;
}



inline RNScalar NeuronPredictionBoundary::
Skew(void) const
{
   if (skew == RN_UNKNOWN) {
      ((NeuronPredictionBoundary *) this)->UpdateStatistics();
   }
   // return skew of affinities
   return skew;
}



inline void *NeuronPredictionBoundary::
UserData(void) const
{
   // return user data
   return user_data;
}



/* manipulation functions */

inline void NeuronPredictionBoundary::
SetUserData(void *user_data)
{
   // set user data
   this->user_data = user_data;
}



/* internal functions */

inline unsigned long long NeuronPredictionBoundary::
FileOffset(void) const
{
   // return offset of affinities in .neuron file
   return file_offset;
}



inline void NeuronPredictionBoundary::
SetFileOffset(unsigned long long file_offset)
{
   // set offset of affininties in .neuron file
   this->file_offset = file_offset;
}



inline void NeuronPredictionBoundary::
SetNAffinities(int naffinities)
{
   // set number of affinities
   this->naffinities = naffinities;
}


/* I/O functions */

inline RNBoolean NeuronPredictionBoundary::
AreAffinitiesResident(void) const
{
   // return if affinities are resident
   return read_affinities_count;
}



inline RNBoolean NeuronPredictionBoundary::
AreVoxelsResident(void) const
{
   // return if voxels are resident
   return (read_voxel_count > 0);
}

