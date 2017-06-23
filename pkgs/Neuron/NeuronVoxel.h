// Include file for the neuron voxel class



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

class NeuronVoxel {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronVoxel(void);
   virtual ~NeuronVoxel(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // data access functions
   NeuronData *Data(void) const;
   int DataIndex(void) const;

   // volume access functions
   NeuronSupervoxel *Supervoxel(void) const;
   int SupervoxelIndex(void) const;

   // human label access functions
   NeuronHumanLabel *HumanLabel(void) const;
   int HumanLabelIndex(void) const;

   // prediction access functions
   NeuronPrediction *Prediction(void) const;
   int PredictionIndex(void) const;

   // voxel neighbor access functions
   int NNeighbors(void) const;
   NeuronVoxel *Neighbor(int neighbor_index) const;
   NeuronVoxel *PreviousVoxel(int dim) const;
   NeuronVoxel *NextVoxel(int dim) const;

   // indirect affinity access functions
   RNScalar AverageAffinity(void) const;
   RNScalar AffinityToNeighbor(NeuronVoxel *neighbor) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // position access functions
   R3Point Position(void) const;
   int Coordinate(int dim) const;
   int XCoordinate(void) const;
   int YCoordinate(void) const;
   int ZCoordinate(void) const;

   // affinity access functions
   RNScalar Affinity(int dim) const;
   RNScalar XAffinity(void) const;
   RNScalar YAffinity(void) const;
   RNScalar ZAffinity(void) const;

   // image access functions
   RNScalar Image(void) const;

   // get user data
   void *UserData(void) const;

   // cellular/extracellular label functions
   RNBoolean IsCellular(void) const;
   RNBoolean IsExtracellular(void) const;

   // boundary property functions
   RNBoolean IsOnBoundary(void) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // data manipulation functions
   void SetData(NeuronData *data);

   // position manipulation functions
   void SetPosition(int coordinates[3]);
   void SetPosition(int ix, int iy, int iz);

   // affinity manipulation functions
   void SetAffinities(RNScalar affinities[3]);
   void SetAffinities(RNScalar xaffinity, RNScalar yaffinity, RNScalar zaffinity);

   // image manipulation functions
   void SetImage(RNScalar image);

   // set user data
   void SetUserData(void *user_data);


   ///////////////////////////
   //// DISPLAY FUNCTIONS ////
   ///////////////////////////
public:
   // draw function
   virtual void Draw(void) const;
   virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////
public:
   // I/O functions
   int ReadVoxel(FILE *fp, NeuronData *data);
   int WriteVoxel(FILE *fp);

protected:
   friend class NeuronData;
   NeuronData *data;
   NeuronSupervoxel *supervoxel;
   NeuronHumanLabel *human_label;
   NeuronPrediction *prediction;
   int supervoxel_index;
   int human_label_index;
   int prediction_index;
   int coordinates[3];
   RNScalar affinities[3];
   RNScalar image;
   void *user_data;
};



/////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/* access functions */

inline NeuronData *NeuronVoxel::
Data(void) const
{
   // return data 
   return data;
}



inline int NeuronVoxel::
DataIndex(void) const
{
   rn_assertion(data != NULL);
   // return index in data
   return data->IndicesToIndex(coordinates[RN_X], coordinates[RN_Y], coordinates[RN_Z]);
}



inline NeuronSupervoxel *NeuronVoxel::
Supervoxel(void) const
{
   // return supervoxel that this belongs to
   return supervoxel;
}



inline int NeuronVoxel::
SupervoxelIndex(void) const
{
   // return index in supervoxel
   return supervoxel_index;
}



inline NeuronHumanLabel *NeuronVoxel::
HumanLabel(void) const
{
   // return human label
   return human_label;
}



inline int NeuronVoxel::
HumanLabelIndex(void) const
{
   // return index in human label
   return human_label_index;
}



inline NeuronPrediction *NeuronVoxel::
Prediction(void) const
{
   // return prediction
   return prediction;
}



inline int NeuronVoxel::
PredictionIndex(void) const
{
   // return index in prediction
   return prediction_index;
}



/* property functions */

inline R3Point NeuronVoxel::
Position(void) const
{
   // return position of voxel
   return R3Point(coordinates[RN_X], coordinates[RN_Y], coordinates[RN_Z]);
}



inline int NeuronVoxel::
Coordinate(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return dim coordinate
   return coordinates[dim];
}



inline int NeuronVoxel::
XCoordinate(void) const
{
   // return x coordinate
   return coordinates[RN_X];
}



inline int NeuronVoxel::
YCoordinate(void) const
{
   // return y coordinate
   return coordinates[RN_Y];
}



inline int NeuronVoxel::
ZCoordinate(void) const
{
   // return z coordinate
   return coordinates[RN_Z];
}



inline RNScalar NeuronVoxel::
Affinity(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return dim affinity
   return affinities[dim];
}



inline RNScalar NeuronVoxel::
XAffinity(void) const
{
   // return x affinity
   return affinities[RN_X];
}



inline RNScalar NeuronVoxel::
YAffinity(void) const
{
   // return y affinity
   return affinities[RN_Y];
}



inline RNScalar NeuronVoxel::
ZAffinity(void) const
{
   // return z affinity
   return affinities[RN_Z];
}



inline RNScalar NeuronVoxel::
Image(void) const
{
   // return image
   return image;
}



inline void *NeuronVoxel::
UserData(void) const
{
   // return user data
   return user_data;
}



/* manipulation functions */

inline void NeuronVoxel::
SetData(NeuronData *data)
{
   // set data
   this->data = data;
}



inline void NeuronVoxel::
SetPosition(int coordinates[3])
{
   // set coordinates of voxel
   this->coordinates[RN_X] = coordinates[RN_X];
   this->coordinates[RN_Y] = coordinates[RN_Y];
   this->coordinates[RN_Z] = coordinates[RN_Z];
}



inline void NeuronVoxel::
SetPosition(int ix, int iy, int iz)
{
   // set coordinates of voxel
   this->coordinates[RN_X] = ix;
   this->coordinates[RN_Y] = iy;
   this->coordinates[RN_Z] = iz;
}



inline void NeuronVoxel::
SetAffinities(RNScalar affinities[3])
{
   // set affinities of voxel
   this->affinities[RN_X] = affinities[RN_X];
   this->affinities[RN_Y] = affinities[RN_Y];
   this->affinities[RN_Z] = affinities[RN_Z];
}



inline void NeuronVoxel::
SetAffinities(RNScalar xaffinity, RNScalar yaffinity, RNScalar zaffinity)
{
   // set affinities of voxel
   this->affinities[RN_X] = xaffinity;
   this->affinities[RN_Y] = yaffinity;
   this->affinities[RN_Z] = zaffinity;
}



inline void NeuronVoxel::
SetImage(RNScalar image)
{
   // set the image value
   this->image = image;
}



inline void NeuronVoxel::
SetUserData(void *user_data)
{
   // set user data
   this->user_data = user_data;
}