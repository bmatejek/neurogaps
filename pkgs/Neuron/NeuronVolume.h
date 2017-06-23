// Include file for the neuron volume class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronVolume {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronVolume(void);
   virtual ~NeuronVolume(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // data access functions
   NeuronData *Data(void) const;
   int DataIndex(void) const;

   // voxel access functions
   int NVoxels(void) const;
   int LocalIndex(int global_index) const;
   int GlobalIndex(int local_index) const;
   virtual NeuronVoxel *LocalVoxel(int local_index) const;
   virtual NeuronVoxel *GlobalVoxel(int global_index) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // data location functions
   RNBoolean IsOnBoundary(void) const;
   R3Point Centroid(void) const;

   // spatial property functions
   const R3Box& BBox(void) const;

   // get user data
   void *UserData(void) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // voxel manipulation functions
   void SetNVoxels(int nvoxels);
   void SetVoxelMapping(int local_index, int global_index);

   // boundary instance variable
   void SetBoundaryFlag(RNBoolean on_boundary);

   // spatial property manipulation
   void SetBBox(const R3Box& input_bbox);

   // set user data
   void SetUserData(void *input_user_data);


   ///////////////////////////
   //// DISPLAY FUNCTIONS ////
   ///////////////////////////
public:
   // draw functions
   virtual void Draw(void) const;

   // print functions
   virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////
protected:
   // update functions
   virtual void UpdateBBox(void);
   void InvalidateBBox(void);
   RNBoolean DoesBBoxNeedUpdate(void) const;

public:
   // I/O functions
   RNBoolean AreVoxelsResident(void) const;
   unsigned int ReadCount(void) const;
   virtual int ReadVoxels(void);
   virtual int ReadVoxels(FILE *fp, RNBoolean seek = TRUE);
   virtual int ReleaseVoxels(void);

   // class type definitions
   RN_CLASS_TYPE_DECLARATIONS(NeuronVolume);

protected:
   friend class NeuronData;
   NeuronData *data;
   int data_index;
   RNBoolean boundary_flag;
   unsigned int read_voxel_count;
   int *voxel_mapping;
   int nvoxels;
   R3Box bbox;
   void *user_data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

/* access functions */

inline NeuronData *NeuronVolume::
Data(void) const
{
   // return neuron data
   return data;
}



inline int NeuronVolume::
DataIndex(void) const
{
   // return data index
   return data_index;
}



inline int NeuronVolume::
NVoxels(void) const
{
   // return number of voxels
   return nvoxels;
}



/* property functions */

inline RNBoolean NeuronVolume::
IsOnBoundary(void) const
{
   // return if volume is on boundary
   return boundary_flag;
}



inline R3Point NeuronVolume::
Centroid(void) const
{
   // return center of bounding box
   return bbox.Centroid();
}



inline const R3Box& NeuronVolume::
BBox(void) const
{
   // return bounding box
   return bbox;
}



inline void *NeuronVolume::
UserData(void) const
{
   // return user data
   return user_data;
}



/* manipulation functions */

inline void NeuronVolume::
SetNVoxels(int nvoxels)
{
   if (voxel_mapping) delete[] voxel_mapping;

   // set number of voxels
   this->nvoxels = nvoxels;

   // allocate memory for voxel mapping
   voxel_mapping = new int[nvoxels];
}



inline void NeuronVolume::
SetVoxelMapping(int local_index, int global_index)
{
   rn_assertion(voxel_mapping != NULL);
   rn_assertion((0 <= local_index) && (local_index < nvoxels));
   rn_assertion((0 <= global_index) && (global_index < data->NVoxels()));
   // set voxel mapping of local index to global index
   voxel_mapping[local_index] = global_index;
}



inline void NeuronVolume::
SetBoundaryFlag(RNBoolean on_boundary)
{
   // set the boundary flag
   this->boundary_flag = on_boundary;
}


inline void NeuronVolume::
SetBBox(const R3Box& bbox)
{
   // set bounding box
   this->bbox = bbox;
}



inline void NeuronVolume::
SetUserData(void *user_data)
{
   // set user data
   this->user_data = user_data;
}



/* I/O functions */

inline RNBoolean NeuronVolume::
AreVoxelsResident(void) const
{
   // return if voxels are resident
   return (read_voxel_count > 0);
}



inline unsigned int NeuronVolume::
ReadCount(void) const
{
   // return the read count
   return read_voxel_count;
}