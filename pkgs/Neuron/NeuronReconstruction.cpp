// Source file for the neuron reconstruction class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronReconstruction);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronReconstruction::
NeuronReconstruction(void) :
supervoxels()
{
}



NeuronReconstruction::
~NeuronReconstruction(void)
{
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

NeuronVoxel *NeuronReconstruction::
LocalVoxel(int local_index) const
{
   rn_assertion(voxel_mapping);
   // return voxel at local index
   int global_voxel = voxel_mapping[local_index];
   return data->Voxel(global_voxel);
}



NeuronVoxel *NeuronReconstruction::
GlobalVoxel(int global_index) const
{
   // return global voxel
   return data->Voxel(global_index);
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void NeuronReconstruction::
InsertSupervoxel(NeuronSupervoxel *supervoxel)
{
   // overridden
}



void NeuronReconstruction::
RemoveSupervoxel(NeuronSupervoxel *supervoxel)
{
   // overriden
}



void NeuronReconstruction::
CreateVoxelMapping(int nvoxels)
{
   rn_assertion(voxel_mapping == NULL);
   rn_assertion(nvoxels != 0);
   // create voxel mapping array
   voxel_mapping = new int[nvoxels];
   this->nvoxels = nvoxels;
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronReconstruction::
Draw(void) const
{
   // overridden
}



void NeuronReconstruction::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // overridden
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void NeuronReconstruction::
UpdateBBox(void)
{
   rn_assertion(AreVoxelsResident());

   // update bounding box
   bbox = R3null_box;
   for (int iv = 0; iv < nvoxels; ++iv) {
      NeuronVoxel *voxel = LocalVoxel(iv);
      bbox.Union(voxel->Position());
   }
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int NeuronReconstruction::
ReadVoxels(void)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;

   // read all supervoxel voxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (!Supervoxel(is)->ReadVoxels()) return 0;
   }

   // return success
   return 1;
}



int NeuronReconstruction::
ReadVoxels(FILE *fp, RNBoolean seek)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;

   // read all supervoxel voxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (!Supervoxel(is)->ReadVoxels(fp, seek)) return 0;
   }

   // return success
   return 1;
}



int NeuronReconstruction::
ReleaseVoxels(void)
{
   // check read counter
   if (--read_voxel_count > 0) return 1;

   // release all of the supervoxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (!Supervoxel(is)->ReleaseVoxels()) return 0;
   }

   // return success
   return 1;
}