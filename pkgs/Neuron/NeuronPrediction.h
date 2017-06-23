// Include file for the neuron prediction class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronPrediction : public NeuronReconstruction {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronPrediction(void);
   virtual ~NeuronPrediction(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // boundary access functions
   int NBoundaries(void) const;
   NeuronPredictionBoundary *Boundary(int boundary_index) const;
   NeuronPredictionBoundary *Boundary(NeuronPrediction *neighbor) const;

   // summary access functions
   int MajorityHumanLabelIndex(void) const;
   NeuronHumanLabel *MajorityHumanLabel(void) const;

   int CenterVoxelIndex(void) const;
   NeuronVoxel *CenterVoxel(void) const;

   int CenterHumanLabelIndex(void) const;
   NeuronHumanLabel *CenterHumanLabel(void) const;

   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // boundary manipulation functions
   void InsertBoundary(NeuronPredictionBoundary *boundary);
   void RemoveBoundary(NeuronPredictionBoundary *boundary);

   // supervoxel manipulation functions
   virtual void InsertSupervoxel(NeuronSupervoxel *supervoxel);
   virtual void RemoveSupervoxel(NeuronSupervoxel *supervoxel);


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
public:
   // class type definitions
   RN_CLASS_TYPE_DECLARATIONS(NeuronPrediction);

protected:
   friend class NeuronData;
   RNArray<NeuronPredictionBoundary *> boundaries;
   int majority_human_label_index;
   int center_human_label_index;
   int center_voxel_index;
};



inline int NeuronPrediction::
NBoundaries(void) const
{
   // return number of boundaries
   return boundaries.NEntries();
}



inline NeuronPredictionBoundary *NeuronPrediction::
Boundary(int boundary_index) const
{
   rn_assertion((0 <= boundary_index) && (boundary_index < boundaries.NEntries()));
   // return kth boundary
   return boundaries.Kth(boundary_index);
}


inline int NeuronPrediction::
CenterVoxelIndex(void) const
{
   // return the index of the center voxel
   return center_voxel_index;
}



inline NeuronVoxel *NeuronPrediction::
CenterVoxel(void) const
{
   rn_assertion(AreVoxelsResident());
   // return the center voxel
   return LocalVoxel(center_voxel_index);
}



inline int NeuronPrediction::
MajorityHumanLabelIndex(void) const
{
   // return the index of the most prominent human label
   return majority_human_label_index;
}



inline NeuronHumanLabel *NeuronPrediction::
MajorityHumanLabel(void) const
{
   // return the most prominent human label
   if (majority_human_label_index == -1) return NULL;
   else return data->HumanLabel(majority_human_label_index);
}


inline int NeuronPrediction::
CenterHumanLabelIndex(void) const
{
   // return the index of the human label
   return center_human_label_index;
}



inline NeuronHumanLabel *NeuronPrediction::
CenterHumanLabel(void) const
{
   // return the human label
   if (center_human_label_index == -1) return NULL;
   return data->HumanLabel(center_human_label_index);
}



