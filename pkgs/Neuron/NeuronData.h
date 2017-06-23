// Include file for the neuron data class


// useful constants



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

class NeuronData {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronData(void);
   virtual ~NeuronData(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // voxel access functions
   int NVoxels(void) const;
   NeuronVoxel *Voxel(int ix, int iy, int iz) const;
   NeuronVoxel *Voxel(int voxel_index) const;

   // human label access functions
   int NHumanLabels(void) const;
   NeuronHumanLabel *HumanLabel(int human_label_index) const;

   // prediction access functions
   int NPredictions(void) const;
   NeuronPrediction *Prediction(int prediction_index) const;

   // supervoxel access functions
   int NSupervoxels(void) const;
   NeuronSupervoxel *Supervoxel(int supervoxel_index) const;

   // cellular access functions
   int NCellulars(void) const;
   NeuronCellular *Cellular(int cellular_index) const;

   // extracellular access functions
   int NExtracellulars(void) const;
   NeuronExtracellular *Extracellular(int extracellular_index) const;

   // prediction boundary access functions
   int NPredictionBoundaries(void) const;
   NeuronPredictionBoundary *PredictionBoundary(int prediction_boundary_index) const;

   // boundary access functions
   int NBoundaries(void) const;
   NeuronBoundary *Boundary(int boundary_index) const;

public:
   // index/indices conversion functions
   int IndicesToIndex(const int coordinates[3]) const;
   int IndicesToIndex(int ix, int iy, int iz) const;
   void IndicesToIndex(const int coordinates[3], int& index) const;
   void IndicesToIndex(int ix, int iy, int iz, int& index) const;
   void IndexToIndices(int index, int(&coordinates)[3]) const;
   void IndexToIndices(int index, int& ix, int& iy, int& iz) const;

   // voxel mapping functions
   int VoxelMapping(int voxel_index) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // resolution property functions
   int XResolution(void) const;
   int YResolution(void) const;
   int ZResolution(void) const;
   int Resolution(int dim) const;

   // scaling property functions
   RNScalar XScale(void) const;
   RNScalar YScale(void) const;
   RNScalar ZScale(void) const;
   RNScalar Scale(int dim) const;

   // get bounding box
   const R3Box& GridBox(void);
   R3Box WorldBox(void);

   // get affine transformation
   const R3Affine& Transformation(void) const;

   // world position functions
   R3Point GridPosition(const R3Point& world_position) const;
   R3Point WorldPosition(const R3Point& grid_position) const;
   R3Point WorldPosition(int x, int y, int z) const;

   // get user data
   void *UserData(void) const;

   // affinity property functions
   int AffinityStatistics(RNScalar(&mean)[4], RNScalar(&stddev)[4]) const;

   // return the rand error at voxel level
   int NTestMetrics(void) const;
   const char *TestMetricName(int index) const;
   int TestMetric(int *voxel_proposals, RNScalar *metrics, int print_verbose = TRUE) const;

   // return segmentation errors at supervoxel level
   int NSegmentationMetrics(void) const;
   const char *SegmentationMetricName(int index) const;
   int SegmentationMetric(int *cellular_proposals, RNScalar *metrics) const;


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////

   // cellular manipulation
   void InsertCellular(NeuronCellular *cellular);
   void RemoveCellular(NeuronCellular *cellular);

   // extracellular manipulation
   void InsertExtracellular(NeuronExtracellular *extracellular);
   void RemoveExtracellular(NeuronExtracellular *extracellular);

   // human label manipulation
   void InsertHumanLabel(NeuronHumanLabel *human_label);
   void RemoveHumanLabel(NeuronHumanLabel *human_label);

   // prediction manipulation
   void InsertPrediction(NeuronPrediction *prediction);
   void RemovePrediction(NeuronPrediction *prediction);

   // prediction boundary manipulation
   void InsertPredictionBoundary(NeuronPredictionBoundary *prediction_boundary);
   void RemovePredictionBoundary(NeuronPredictionBoundary *prediction_boundary);

   // boundary manipulation
   void InsertBoundary(NeuronBoundary *boundary);
   void RemoveBoundary(NeuronBoundary *boundary);

   // affinity manipulation
   int NormalizeAffinities(void);
   int FilterBoundaryAffinities(int metric);
   int FilterAffinity(const char *filter_type, const char *extension = ".txt");


   ///////////////////////////
   //// DISPLAY FUNCTIONS ////
   ///////////////////////////




   ///////////////////////
   //// I/O FUNCTIONS ////
   ///////////////////////

   // grid functions
   int ReadAffinities(R3Grid *affinities_grids[3]) const;
   R3Grid *ReadTruthTxtFile(const char human_labels_txt_filename[4096]) const;
   int ReadHumanLabels(R3Grid *human_label_grid) const;
   int ReadImages(R3Grid *image_grid) const;
   int ReadMachineLabels(R3Grid *machine_labels_grid) const;
   int ReadPredictions(R3Grid *predictions_grid) const;

   // read/write functions
   virtual int ReadFile(const char *input_filename, RNBoolean read_voxels = FALSE, RNBoolean read_affinities = FALSE);
   virtual int ReadASCIIFile(const char *ascii_filename, RNBoolean read_voxels = FALSE, RNBoolean read_affinities = FALSE);
   virtual int ReadNeuronFile(const char *neuron_filename, RNBoolean read_voxels = FALSE, RNBoolean read_affinities = FALSE);
   virtual int WriteFile(const char *output_filename);
   virtual int WriteASCIIFile(const char *ascii_filename);
   virtual int WriteNeuronFile(const char *neuron_filename);


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////

public:
   // update functions
   void UpdateBBox(void);
   void InvalidateBBox(void);
   RNBoolean DoesBBoxNeedUpdate(void) const;
   void SetBBox(const R3Box& bbox);
   int ReadPredictions(const char meta_filename[4096]);

   // ASCII to Neuron helper functions
   void LabelExtracellulars(R3Grid *machine_labels, int ncellulars) const;
   void CreateSupervoxels(R3Grid *machine_labels);
   void SetAffinities(R3Grid *affinties[3]);
   void SetImageValues(R3Grid *image);
   void CreateHumanLabels(R3Grid *human_labels);
   void CreatePredictions(R3Grid *predictions);
   void CreatePredictionMapping(void);
   void CreatePredictionBoundaries(void);
   void CreateBoundaries(void);
   void CreateHumanLabelMapping(void);
   void CreateSupervoxelMapping(void);

   // filename functions
   const char *FilePath(void) const;
   const char *Filename(void) const;
   const char *AffinitiesFilename(void) const;
   const char *HumanLabelsFilename(void) const;
   const char *ImageFilename(void) const;
   const char *MachineLabelsFilename(void) const;

   // I/O functions
   RNBoolean AreVoxelsResident(void) const;
   unsigned int ReadCount(void) const;
   int ReadVoxels(void);
   int ReleaseVoxels(void);

protected:
   RNArray<NeuronHumanLabel *> human_labels;
   RNArray<NeuronPrediction *> predictions;
   RNArray<NeuronCellular *> cellulars;
   RNArray<NeuronExtracellular *> extracellulars;
   RNArray<NeuronPredictionBoundary *> prediction_boundaries;
   RNArray<NeuronBoundary *> boundaries;
   unsigned int read_voxel_count;
   int *voxel_mapping;
   int resolution[3];
   RNScalar scaling[3];
   R3Affine transformation;
   int data_size;
   int data_sheet_size;
   int data_row_size;
   char file_path[4096];
   char filename[4096];
   char affinities_filename[4096];
   char human_labels_filename[4096];
   char image_filename[4096];
   char machine_labels_filename[4096];
   R3Box bbox;
   void *user_data;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

/* access functions */

inline int NeuronData::
NHumanLabels(void) const
{
   // return number of human labels
   return human_labels.NEntries();
}



inline NeuronHumanLabel *NeuronData::
HumanLabel(int human_label_index) const
{
   rn_assertion((0 <= human_label_index) && (human_label_index < human_labels.NEntries()));
   // return kth human label
   return human_labels.Kth(human_label_index);
}



inline int NeuronData::
NPredictions(void) const
{
   // return number of predictions
   return predictions.NEntries();
}



inline NeuronPrediction *NeuronData::
Prediction(int prediction_index) const
{
   rn_assertion((0 <= prediction_index) && (prediction_index < predictions.NEntries()));
   // return kth prediction
   return predictions.Kth(prediction_index);
}



inline int NeuronData::
NSupervoxels(void) const
{
   // return number of supervoxels
   return cellulars.NEntries() + extracellulars.NEntries();
}



inline NeuronSupervoxel *NeuronData::
Supervoxel(int supervoxel_index) const
{
   rn_assertion((0 <= supervoxel_index) && (supervoxel_index < cellulars.NEntries() + extracellulars.NEntries()));
   // return Kth supervoxel
   if (supervoxel_index < cellulars.NEntries()) return (NeuronSupervoxel *)cellulars.Kth(supervoxel_index);
   else return (NeuronSupervoxel *)extracellulars.Kth(supervoxel_index - cellulars.NEntries());
}



inline int NeuronData::
NCellulars(void) const
{
   // return number of cellulars
   return cellulars.NEntries();
}



inline NeuronCellular *NeuronData::
Cellular(int cellular_index) const
{
   rn_assertion((0 <= cellular_index) && (cellular_index < cellulars.NEntries()));
   // return kth cellular
   return cellulars.Kth(cellular_index);
}



inline int NeuronData::
NExtracellulars(void) const
{
   // return number of extracellulars
   return extracellulars.NEntries();
}



inline NeuronExtracellular *NeuronData::
Extracellular(int extracellular_index) const
{
   rn_assertion((0 <= extracellular_index) && (extracellular_index < extracellulars.NEntries()));
   // return kth extracellular
   return extracellulars.Kth(extracellular_index);
}



inline int NeuronData::
NBoundaries(void) const
{
   // return number of boundaries
   return boundaries.NEntries();
}



inline NeuronBoundary *NeuronData::
Boundary(int boundary_index) const
{
   rn_assertion((0 <= boundary_index) && (boundary_index < boundaries.NEntries()));
   // return kth boundary
   return boundaries.Kth(boundary_index);
}


inline int NeuronData::
NPredictionBoundaries(void) const
{
   // return number of boundaries
   return prediction_boundaries.NEntries();
}



inline NeuronPredictionBoundary *NeuronData::
PredictionBoundary(int prediction_boundary_index) const
{
   rn_assertion((0 <= prediction_boundary_index) && (prediction_boundary_index < prediction_boundaries.NEntries()));
   // return kth boundary
   return prediction_boundaries.Kth(prediction_boundary_index);
}




inline int NeuronData::
IndicesToIndex(const int coordinates[3]) const
{
   rn_assertion((0 <= coordinates[RN_X]) && (coordinates[RN_X] < resolution[RN_X]));
   rn_assertion((0 <= coordinates[RN_Y]) && (coordinates[RN_Y] < resolution[RN_Y]));
   rn_assertion((0 <= coordinates[RN_Z]) && (coordinates[RN_Z] < resolution[RN_Z]));
   // return the index cooresponding to this point
   return coordinates[RN_Z] * data_sheet_size + coordinates[RN_Y] * data_row_size + coordinates[RN_X];
}



inline int NeuronData::
IndicesToIndex(int ix, int iy, int iz) const
{
   rn_assertion((0 <= ix) && (ix < resolution[RN_X]));
   rn_assertion((0 <= iy) && (iy < resolution[RN_Y]));
   rn_assertion((0 <= iz) && (iz < resolution[RN_Z]));
   // return the index cooresponding to these coordinates
   return iz * data_sheet_size + iy * data_row_size + ix;
}



inline void NeuronData::
IndicesToIndex(const int coordinates[3], int& index) const
{
   rn_assertion((0 <= coordinates[RN_X]) && (coordinates[RN_X] < resolution[RN_X]));
   rn_assertion((0 <= coordinates[RN_Y]) && (coordinates[RN_Y] < resolution[RN_Y]));
   rn_assertion((0 <= coordinates[RN_Z]) && (coordinates[RN_Z] < resolution[RN_Z]));
   // set the index
   index = coordinates[RN_Z] * data_sheet_size + coordinates[RN_Y] * data_row_size + coordinates[RN_X];
}



inline void NeuronData::
IndicesToIndex(int ix, int iy, int iz, int& index) const
{
   rn_assertion((0 <= ix) && (ix < resolution[RN_X]));
   rn_assertion((0 <= iy) && (iy < resolution[RN_Y]));
   rn_assertion((0 <= iz) && (iz < resolution[RN_Z]));
   // set the index
   index = iz * data_sheet_size + iy * data_row_size + ix;
}



inline void NeuronData::
IndexToIndices(int index, int(&coordinates)[3]) const
{
   rn_assertion((0 <= index) && (index < data_size));
   // set indices at index
   coordinates[RN_Z] = index / data_sheet_size;
   coordinates[RN_Y] = (index - coordinates[RN_Z] * data_sheet_size) / data_row_size;
   coordinates[RN_X] = index % data_row_size;
}



inline void NeuronData::
IndexToIndices(int index, int& ix, int& iy, int& iz) const
{
   rn_assertion((0 <= index) && (index < data_size));
   // set indices at index
   iz = index / data_sheet_size;
   iy = (index - iz * data_sheet_size) / data_row_size;
   ix = index % data_row_size;
}



inline int NeuronData::
VoxelMapping(int voxel_index) const
{
   rn_assertion(voxel_mapping != NULL);
   rn_assertion((0 <= voxel_index) && (voxel_index < data_size));
   // return voxel mapping for this index
   return voxel_mapping[voxel_index];
}



/* property functions */

inline int NeuronData::
XResolution(void) const
{
   // return x resolution
   return resolution[RN_X];
}



inline int NeuronData::
YResolution(void) const
{
   // return y resolution
   return resolution[RN_Y];
}



inline int NeuronData::
ZResolution(void) const
{
   // return z resolution
   return resolution[RN_Z];
}



inline int NeuronData::
Resolution(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return dim resolution
   return resolution[dim];
}



inline RNScalar NeuronData::
XScale(void) const
{
   // return x scale
   return scaling[RN_X];
}



inline RNScalar NeuronData::
YScale(void) const
{
   // return y scale
   return scaling[RN_Y];
}



inline RNScalar NeuronData::
ZScale(void) const
{
   // return z scale
   return scaling[RN_Z];
}



inline RNScalar NeuronData::
Scale(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return dim scaling
   return scaling[dim];
}



inline const R3Box& NeuronData::
GridBox(void)
{
   // return bounding box
   if (bbox.XMin() == FLT_MAX)
      this->UpdateBBox();
   return bbox;
}



inline R3Box NeuronData::
WorldBox(void)
{
   // return bounding box in world coordinates
   R3Box world_box = GridBox();
   world_box.Transform(Transformation());
   return world_box;
}



inline const R3Affine& NeuronData::
Transformation(void) const
{
   // return affine transformation
   return transformation;
}



inline R3Point NeuronData::
GridPosition(const R3Point& world_position) const
{
   // return the grid position from the world position
   R3Point grid_position = R3Point(world_position);
   grid_position.InverseTransform(transformation);
   return grid_position;
}



inline R3Point NeuronData::
WorldPosition(const R3Point& grid_position) const
{
   // return world position from the grid position
   R3Point world_position = R3Point(grid_position);
   world_position.Transform(transformation);
   return world_position;
}



inline R3Point NeuronData::
WorldPosition(int x, int y, int z) const
{
   // return world position from voxel indices
   return WorldPosition(R3Point(x, y, z));
}



inline void *NeuronData::
UserData(void) const
{
   // return user data
   return user_data;
}



/* internal functions */

inline const char *NeuronData::
FilePath(void) const
{
   // return directory
   return file_path;
}



inline const char *NeuronData::
Filename(void) const
{
   // return filename
   return filename;
}



inline const char *NeuronData::
AffinitiesFilename(void) const
{
   // return affinities filename
   return affinities_filename;
}



inline const char *NeuronData::
HumanLabelsFilename(void) const
{
   // return human labels filename
   return human_labels_filename;
}



inline const char *NeuronData::
ImageFilename(void) const
{
   // return image filename
   return image_filename;
}



inline const char *NeuronData::
MachineLabelsFilename(void) const
{
   // return machine labels filename
   return machine_labels_filename;
}



/* I/O functions */

inline RNBoolean NeuronData::
AreVoxelsResident(void) const
{
   // return if voxels are resident
   return (read_voxel_count > 0);
}



inline unsigned int NeuronData::
ReadCount(void) const
{
   // return the read count
   return read_voxel_count;
}
