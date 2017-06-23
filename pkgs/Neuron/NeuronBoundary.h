// Include file for the neuron boundary class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronBoundary {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronBoundary(void);
   virtual ~NeuronBoundary(void);


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
   int SupervoxelOneIndex(void) const;
   int SupervoxelTwoIndex(void) const;

   // return the volumes
   NeuronSupervoxel *SupervoxelOne(void) const;
   NeuronSupervoxel *SupervoxelTwo(void) const;
   NeuronSupervoxel *OtherSupervoxel(const NeuronSupervoxel *supervoxel) const;


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
   friend class NeuronSupervoxel;
   NeuronData *data;
   int data_index;
   int naffinities;
   unsigned long long file_offset;
   unsigned int read_voxel_count;
   unsigned int read_affinities_count;
   RNScalar *affinities;
   NeuronSupervoxel *supervoxel_one;
   NeuronSupervoxel *supervoxel_two;
   int supervoxel_one_index;
   int supervoxel_two_index;
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

inline NeuronData *NeuronBoundary::
Data(void) const
{
   // return neuron data
   return data;
}



inline int NeuronBoundary::
DataIndex(void) const
{
   // return the index of this boundary
   return data_index;
}



inline int NeuronBoundary::
NAffinities(void) const
{
   // return number of affinities
   return naffinities;
}



inline RNScalar NeuronBoundary::
Affinity(int affinity_index) const
{
   rn_assertion((0 <= affinity_index) && (affinity_index < naffinities));
   rn_assertion(affinities != NULL);
   // return kth affinity
   return affinities[affinity_index];
}



inline int NeuronBoundary::
SupervoxelOneIndex(void) const
{
   // return the index of the first volume
   return supervoxel_one_index;
}



inline int NeuronBoundary::
SupervoxelTwoIndex(void) const
{
   // return the index of the second volume
   return supervoxel_two_index;
}



inline NeuronSupervoxel *NeuronBoundary::
SupervoxelOne(void) const
{
   // return the first volume
   return supervoxel_one;
}



inline NeuronSupervoxel *NeuronBoundary::
SupervoxelTwo(void) const
{
   // return the second volume
   return supervoxel_two;
}



inline NeuronSupervoxel *NeuronBoundary::
OtherSupervoxel(const NeuronSupervoxel *supervoxel) const
{
   rn_assertion(supervoxel != NULL);
   // return other volume
   if (supervoxel_one == supervoxel) return supervoxel_two;
   else if (supervoxel_two == supervoxel) return supervoxel_one;
   else return NULL;
}



/* property functions */

inline RNScalar NeuronBoundary::
Maximum(void) const
{
   if (maximum == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return the maximum
   return maximum;
}



inline RNScalar NeuronBoundary::
Mean(void) const
{
   if (mean == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return the average
   return mean;
}



inline RNScalar NeuronBoundary::
Median(void) const
{
   if (median == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return the median
   return median;
}



inline RNScalar NeuronBoundary::
Minimum(void) const
{
   if (minimum == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return the minimum
   return minimum;
}



inline RNScalar NeuronBoundary::
StdDev(void) const
{
   if (stddev == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return the standard deviation
   return stddev;
}



inline RNScalar NeuronBoundary::
Skew(void) const
{
   if (skew == RN_UNKNOWN) {
      ((NeuronBoundary *) this)->UpdateStatistics();
   }
   // return skew of affinities
   return skew;
}



inline void *NeuronBoundary::
UserData(void) const
{
   // return user data
   return user_data;
}



/* manipulation functions */

inline void NeuronBoundary::
SetUserData(void *user_data)
{
   // set user data
   this->user_data = user_data;
}



/* internal functions */

inline unsigned long long NeuronBoundary::
FileOffset(void) const
{
   // return offset of affinities in .neuron file
   return file_offset;
}



inline void NeuronBoundary::
SetFileOffset(unsigned long long file_offset)
{
   // set offset of affininties in .neuron file
   this->file_offset = file_offset;
}



inline void NeuronBoundary::
SetNAffinities(int naffinities)
{
   // set number of affinities
   this->naffinities = naffinities;
}


/* I/O functions */

inline RNBoolean NeuronBoundary::
AreAffinitiesResident(void) const
{
   // return if affinities are resident
   return read_affinities_count;
}



inline RNBoolean NeuronBoundary::
AreVoxelsResident(void) const
{
   // return if voxels are resident
   return (read_voxel_count > 0);
}

