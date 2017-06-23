// Source file for the neuron voxel class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



////////////////////////////////////////////////////////////////////////
// Useful constants
////////////////////////////////////////////////////////////////////////

#define VOXEL_FILLER 1
#define NNEIGHBORS 6



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronVoxel::
NeuronVoxel(void) :
data(NULL),
supervoxel(NULL),
human_label(NULL),
prediction(NULL),
supervoxel_index(-1),
human_label_index(-1),
prediction_index(-1),
user_data(NULL)
{
   for (int dim = 0; dim <= 2; ++dim) {
      coordinates[dim] = -1;
      affinities[dim] = -1;
   }
}



NeuronVoxel::
~NeuronVoxel(void)
{
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int NeuronVoxel::
NNeighbors(void) const
{
   // return the number of possible neighbors
   return NNEIGHBORS;
}



NeuronVoxel *NeuronVoxel::
Neighbor(int neighbor_index) const {
   rn_assertion((0 <= neighbor_index) && (neighbor_index < NNEIGHBORS));
   // get the offset from this voxel
   int dx = 0, dy = 0, dz = 0;
   if (neighbor_index == 0) dx = -1;
   else if (neighbor_index == 1) dx = 1;
   else if (neighbor_index == 2) dy = -1;
   else if (neighbor_index == 3) dy = 1;
   else if (neighbor_index == 4) dz = -1;
   else if (neighbor_index == 5) dz = 1;

   // get the neighbor coordinates
   int neighbor_x = coordinates[RN_X] + dx;
   int neighbor_y = coordinates[RN_Y] + dy;
   int neighbor_z = coordinates[RN_Z] + dz;

   // see if this neighbor is within bounds
   if (neighbor_x < 0) return NULL;
   else if (neighbor_x > data->XResolution() - 1) return NULL;
   else if (neighbor_y < 0) return NULL;
   else if (neighbor_y > data->YResolution() - 1) return NULL;
   else if (neighbor_z < 0) return NULL;
   else if (neighbor_z > data->ZResolution() - 1) return NULL;

   // return the voxel offset from this coordinate
   else return data->Voxel(neighbor_x, neighbor_y, neighbor_z);
}



NeuronVoxel *NeuronVoxel::
PreviousVoxel(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return the previous voxel along dimension
   return Neighbor(2 * dim);
}



NeuronVoxel *NeuronVoxel::
NextVoxel(int dim) const
{
   rn_assertion((0 <= dim) && (dim <= 2));
   // return the next voxel along dimension
   return Neighbor(1 + 2 * dim);
}



RNScalar NeuronVoxel::
AverageAffinity(void) const
{
   int nneighbors = 0;
   RNScalar average = 0.0;
   for (int in = 0; in < NNEIGHBORS; ++in) {
      NeuronVoxel *neighbor = Neighbor(in);
      if (!neighbor) continue;
      nneighbors++;
      average += AffinityToNeighbor(neighbor);
   }

   // return the average affinity
   return average / nneighbors;
}



RNScalar NeuronVoxel::
AffinityToNeighbor(NeuronVoxel *neighbor) const
{
   rn_assertion(neighbor != NULL);
   // get change between neighbor and this
   int dx = neighbor->coordinates[RN_X] - coordinates[RN_X];
   int dy = neighbor->coordinates[RN_Y] - coordinates[RN_Y];
   int dz = neighbor->coordinates[RN_Z] - coordinates[RN_Z];

   // make sure neighbor and this differ by one dimension
   rn_assertion(abs(dx) + abs(dy) + abs(dz) == 1);

   // return either neighbor or this affinity
   if (dx == -1) return affinities[RN_X];
   else if (dx == 1) return neighbor->affinities[RN_X];
   else if (dy == -1) return affinities[RN_Y];
   else if (dy == 1) return neighbor->affinities[RN_Y];
   else if (dz == -1) return affinities[RN_Z];
   else return neighbor->affinities[RN_Z];
}



////////////////////////////////////////////////////////////////////////
// Property functions
////////////////////////////////////////////////////////////////////////

RNBoolean NeuronVoxel::
IsCellular(void) const
{
   // return if voxel is cellular
   return supervoxel->IsCellular();
}



RNBoolean NeuronVoxel::
IsExtracellular(void) const
{
   // return if voxel is extracellular
   return supervoxel->IsExtracellular();
}



RNBoolean NeuronVoxel::
IsOnBoundary(void) const
{
   rn_assertion(data != NULL);
   // check all boundaries
   if (coordinates[RN_X] == 0) return TRUE;
   else if (coordinates[RN_X] == data->XResolution() - 1) return TRUE;
   else if (coordinates[RN_Y] == 0) return TRUE;
   else if (coordinates[RN_Y] == data->YResolution() - 1) return TRUE;
   else if (coordinates[RN_Z] == 0) return TRUE;
   else if (coordinates[RN_Z] == data->ZResolution() - 1) return TRUE;
   else return FALSE;
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronVoxel::
Draw(void) const
{
   rn_assertion((coordinates[RN_X] != -1) && (coordinates[RN_Y] != -1) && (coordinates[RN_Z] != -1));
   // draw the point
   R3LoadPoint(R3Point(coordinates[RN_X], coordinates[RN_Y], coordinates[RN_Z]));
}



void NeuronVoxel::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   rn_assertion((coordinates[RN_X] != -1) && (coordinates[RN_Y] != -1) && (coordinates[RN_Z] != -1));
   // check fp
   if (!fp) fp = stdout;

   // print prefix
   if (prefix) fprintf(fp, "%s", prefix);

   // print voxel coordinates
   fprintf(fp, "(%d, %d, %d)", coordinates[RN_X], coordinates[RN_Y], coordinates[RN_Z]);

   // print suffix
   if (suffix)  fprintf(fp, "%s", suffix);
   fprintf(fp, "\n");
}



////////////////////////////////////////////////////////////////////////
//Input/Output functions
////////////////////////////////////////////////////////////////////////

int NeuronVoxel::
ReadVoxel(FILE *fp, NeuronData *data)
{
   // just checking
   rn_assertion(fp != NULL);
   rn_assertion(data != NULL);

   // convenient variables
   int dummy = 0;

   // read the important voxel attributes
   int supervoxel_data_index;
   int human_label_data_index;
   int prediction_data_index;
   fread(&supervoxel_index, sizeof(int), 1, fp);
   fread(&supervoxel_data_index, sizeof(int), 1, fp);
   fread(&human_label_index, sizeof(int), 1, fp);
   fread(&human_label_data_index, sizeof(int), 1, fp);
   fread(&prediction_index, sizeof(int), 1, fp);
   fread(&prediction_data_index, sizeof(int), 1, fp);
   fread(&(coordinates[0]), sizeof(int), 3, fp);
   fread(&(affinities[0]), sizeof(RNScalar), 3, fp);
   fread(&image, sizeof(RNScalar), 1, fp);
   for (int j = 0; j < VOXEL_FILLER; ++j)
      fread(&dummy, sizeof(int), 1, fp);

   // set data
   this->data = data;

   // set supervoxel, human label, and prediction
   this->supervoxel = data->Supervoxel(supervoxel_data_index);
   if (human_label_data_index != -1)
      this->human_label = data->HumanLabel(human_label_data_index);
   if (prediction_data_index != -1)
      this->prediction = data->Prediction(prediction_data_index);

   // return success
   return 1;
}



int NeuronVoxel::
WriteVoxel(FILE *fp)
{
   // just checking
   rn_assertion(fp != NULL);

   // convenient variables
   int dummy = 0;

   // write the important voxel attributes
   rn_assertion(supervoxel != NULL);
   int supervoxel_data_index = supervoxel->SupervoxelIndex();

   int human_label_data_index;
   if (human_label)
      human_label_data_index = human_label->DataIndex();
   else
      human_label_data_index = -1;

   int prediction_data_index;
   if (prediction)
      prediction_data_index = prediction->DataIndex();
   else
      prediction_data_index = -1;

   fwrite(&supervoxel_index, sizeof(int), 1, fp);
   fwrite(&supervoxel_data_index, sizeof(int), 1, fp);
   fwrite(&human_label_index, sizeof(int), 1, fp);
   fwrite(&human_label_data_index, sizeof(int), 1, fp);
   fwrite(&prediction_index, sizeof(int), 1, fp);
   fwrite(&prediction_data_index, sizeof(int), 1, fp);
   fwrite(&(coordinates[0]), sizeof(int), 3, fp);
   fwrite(&(affinities[0]), sizeof(RNScalar), 3, fp);
   fwrite(&image, sizeof(RNScalar), 1, fp);
   for (int j = 0; j < VOXEL_FILLER; ++j)
      fwrite(&dummy, sizeof(int), 1, fp);

   // return success
   return 1;
}
