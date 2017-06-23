// Include file for the neuron reconstruction class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronReconstruction : public NeuronVolume {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronReconstruction(void);
   ~NeuronReconstruction(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // supervoxel access functions
   int NSupervoxels(void) const;
   NeuronSupervoxel *Supervoxel(int supervoxel_index) const;

   // voxel access functions
   virtual NeuronVoxel *LocalVoxel(int local_index) const;
   virtual NeuronVoxel *GlobalVoxel(int global_index) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // supervoxel manipulation functions
   virtual void InsertSupervoxel(NeuronSupervoxel *supervoxel);
   virtual void RemoveSupervoxel(NeuronSupervoxel *supervoxel);

   // create voxel mapping
   void CreateVoxelMapping(int nvoxels);


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


public:
   // I/O functions
   virtual int ReadVoxels(void);
   virtual int ReadVoxels(FILE *fp, RNBoolean seek = TRUE);
   virtual int ReleaseVoxels(void);

   // class type definitions
   RN_CLASS_TYPE_DECLARATIONS(NeuronReconstruction);

protected:
   friend class NeuronData;
   RNArray<NeuronSupervoxel *> supervoxels;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

/* access functions */

inline int NeuronReconstruction::
NSupervoxels(void) const
{
   // return number of supervoxels
   return supervoxels.NEntries();
}



inline NeuronSupervoxel *NeuronReconstruction::
Supervoxel(int supervoxel_index) const
{
   rn_assertion((0 <= supervoxel_index) && (supervoxel_index < supervoxels.NEntries()));
   // return kth supervoxel
   return supervoxels.Kth(supervoxel_index);
}



