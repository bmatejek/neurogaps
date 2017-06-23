// Source file for the neuron volume class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronVolume);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronVolume::
NeuronVolume(void) :
data(NULL),
data_index(-1),
boundary_flag(FALSE),
read_voxel_count(0),
voxel_mapping(NULL),
nvoxels(0),
bbox(R3null_box),
user_data(NULL)
{
}



NeuronVolume::
~NeuronVolume(void)
{
   // remove from data
   if (voxel_mapping) delete[] voxel_mapping;
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int NeuronVolume::
LocalIndex(int global_index) const
{
   // just checking
   rn_assertion(voxel_mapping != NULL);

   // do a binary search over voxel mapping
   int local_index = nvoxels / 2;
   int search_low = 0;
   int search_high = nvoxels - 1;

   while (voxel_mapping[local_index] != global_index) {
      if (local_index == search_low) break;
      local_index = (search_high - search_low) / 2 + search_low;
      if (voxel_mapping[local_index] < global_index)
         search_low = local_index + 1;
      else
         search_high = local_index - 1;
   }

   // make sure that there is a match
   rn_assertion(voxel_mapping[local_index] == global_index);
   return local_index;
}



int NeuronVolume::
GlobalIndex(int local_index) const
{
   rn_assertion((0 <= local_index) && (local_index < nvoxels));
   rn_assertion(voxel_mapping);
   // return voxel mapping at local index
   return voxel_mapping[local_index];
}



NeuronVoxel *NeuronVolume::
LocalVoxel(int local_index) const
{
   // overridden
   return NULL;
}



NeuronVoxel *NeuronVolume::
GlobalVoxel(int global_index) const
{
   // overridden
   return NULL;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronVolume::
Draw(void) const
{
   // overridden
}



void NeuronVolume::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // overridden
}



////////////////////////////////////////////////////////////////////////
// Internal functions
////////////////////////////////////////////////////////////////////////

void NeuronVolume::
UpdateBBox(void)
{
   // overridden
}



void NeuronVolume::
InvalidateBBox(void)
{
   // invalidate bounding box
   bbox = R3null_box;
   if (data) data->InvalidateBBox();
}



RNBoolean NeuronVolume::
DoesBBoxNeedUpdate(void) const
{
   if (bbox.XMin() == FLT_MAX) return TRUE;
   else return FALSE;
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int NeuronVolume::
ReadVoxels(void)
{
   // overridden
   return 0;
}



int NeuronVolume::
ReadVoxels(FILE *fp, RNBoolean seek)
{
   // overridden
   return 0;
}



int NeuronVolume::
ReleaseVoxels(void)
{
   // overridden
   return 0;
}
