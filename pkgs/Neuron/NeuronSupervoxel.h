// Include file for the neuron supervoxel class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronSupervoxel : public NeuronVolume {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronSupervoxel(void);
   virtual ~NeuronSupervoxel(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // supervoxel access functions
   virtual int SupervoxelIndex(void) const;

   // boundary access functions
   int NBoundaries(void) const;
   NeuronBoundary *Boundary(int boundary_index) const;
   NeuronBoundary *Boundary(NeuronSupervoxel *neighbor) const;

   // human label access functions
   int NHumanLabels(void) const;
   NeuronHumanLabel *HumanLabel(int human_label_index) const;

   // summary access functions
   int MajorityHumanLabelIndex(void) const;
   NeuronHumanLabel *MajorityHumanLabel(void) const;

   int CenterVoxelIndex(void) const;
   NeuronVoxel *CenterVoxel(void) const;

   int CenterHumanLabelIndex(void) const;
   NeuronHumanLabel *CenterHumanLabel(void) const;

   // prediction access functions
   int NPredictions(void) const;
   NeuronPrediction *Prediction(int prediction_index) const;

   // voxel access functions
   virtual NeuronVoxel *LocalVoxel(int local_index) const;
   virtual NeuronVoxel *GlobalVoxel(int global_index) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // cellular/extracellular property
   virtual RNBoolean IsCellular(void) const;
   virtual RNBoolean IsExtracellular(void) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // boundary manipulation functions
   void InsertBoundary(NeuronBoundary *boundary);
   void RemoveBoundary(NeuronBoundary *boundary);

   // human label manipulation functions
   int InsertHumanLabel(NeuronHumanLabel *human_label);
   void RemoveHumanLabel(NeuronHumanLabel *human_label);

   // prediction manipulation functions
   int InsertPrediction(NeuronPrediction *prediction);
   void RemovePrediction(NeuronPrediction *prediction);

   // voxel creation functions
   void CreateVoxels(int nvoxels);


   ///////////////////////////
   //// DISPLAY FUNCTIONS ////
   ///////////////////////////
public:
   virtual void Draw(void) const;
   virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////
public:
   // update functions
   virtual void UpdateBBox(void);

   // I/O functions
   virtual int ReadVoxels(void);
   virtual int ReadVoxels(FILE *fp, RNBoolean seek = TRUE);
   virtual int WriteVoxels(FILE *fp, RNBoolean seek = TRUE);
   virtual int ReleaseVoxels(void);

protected:
   // file functions
   unsigned long long FileOffset(void) const;
   void SetFileOffset(unsigned long long file_offset);

public:
   // class type declarations
   RN_CLASS_TYPE_DECLARATIONS(NeuronSupervoxel);

protected:
   friend class NeuronData;
   friend class NeuronHumanLabel;
   friend class NeuronPrediction;
   unsigned long long file_offset;
   NeuronVoxel *voxels;
   RNArray<NeuronHumanLabel *> human_labels;
   RNArray<NeuronPrediction *> predictions;
   RNArray<NeuronBoundary *> boundaries;
   int majority_human_label_index;
   int center_human_label_index;
   int center_voxel_index;
};



/////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/* access functions */

inline int NeuronSupervoxel::
SupervoxelIndex(void) const
{
   // overridden
   return -1;
}



inline int NeuronSupervoxel::
NBoundaries(void) const
{
   // return number of boundaries
   return boundaries.NEntries();
}



inline NeuronBoundary *NeuronSupervoxel::
Boundary(int boundary_index) const
{
   rn_assertion((0 <= boundary_index) && (boundary_index < boundaries.NEntries()));
   // return kth boundary
   return boundaries.Kth(boundary_index);
}



inline int NeuronSupervoxel::
NHumanLabels(void) const
{
   // return number of human labels
   return human_labels.NEntries();
}



inline NeuronHumanLabel *NeuronSupervoxel::
HumanLabel(int human_label_index) const
{
   rn_assertion((0 <= human_label_index) && (human_label_index < human_labels.NEntries()));
   // return kth human label
   return human_labels.Kth(human_label_index);
}



inline int NeuronSupervoxel::
MajorityHumanLabelIndex(void) const
{
   // return the index of the most prominent human label
   return majority_human_label_index;
}



inline NeuronHumanLabel *NeuronSupervoxel::
MajorityHumanLabel(void) const
{
   // return the most prominent human label
   if (majority_human_label_index == -1) return NULL;
   else return data->HumanLabel(majority_human_label_index);
}



inline int NeuronSupervoxel::
CenterVoxelIndex(void) const
{
   // return the index of the center voxel
   return center_voxel_index;
}



inline NeuronVoxel *NeuronSupervoxel::
CenterVoxel(void) const
{
   rn_assertion(AreVoxelsResident());   
   // return the center voxel
   return LocalVoxel(center_voxel_index);
}



inline int NeuronSupervoxel::
CenterHumanLabelIndex(void) const
{
   // return the index of the human label
   return center_human_label_index;
}



inline NeuronHumanLabel *NeuronSupervoxel::
CenterHumanLabel(void) const
{
   // return the human label
   if (center_human_label_index == -1) return NULL;
   return data->HumanLabel(center_human_label_index);
}



inline int NeuronSupervoxel::
NPredictions(void) const
{
   // return number of predictions
   return predictions.NEntries();
}



inline NeuronPrediction *NeuronSupervoxel::
Prediction(int prediction_index) const
{
   rn_assertion((0 <= prediction_index) && (prediction_index < predictions.NEntries()));
   // return kth prediction
   return predictions.Kth(prediction_index);
}



inline NeuronVoxel *NeuronSupervoxel::
LocalVoxel(int local_index) const
{
   rn_assertion((0 <= local_index) && (local_index < nvoxels));
   rn_assertion(voxels);
   // return voxels
   return &voxels[local_index];
}



/* property functions */

inline RNBoolean NeuronSupervoxel::
IsCellular(void) const
{
   // overridden
   return FALSE;
}



inline RNBoolean NeuronSupervoxel::
IsExtracellular(void) const
{
   // overridden
   return FALSE;
}



/* file functions */

inline unsigned long long NeuronSupervoxel::
FileOffset(void) const
{
   // return offset of voxels in .neuron file
   return file_offset;
}



inline void NeuronSupervoxel::
SetFileOffset(unsigned long long file_offset)
{
   // set offset of voxels in .neuron file
   this->file_offset = file_offset;
}
