// Source file for the neuron cellular class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronCellular);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronCellular::
NeuronCellular(void)
{
}



NeuronCellular::
~NeuronCellular(void)
{
   data->RemoveCellular(this);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int NeuronCellular::
SupervoxelIndex(void) const
{
   if (!data) return -1;
   // return offset supervoxel index
   return DataIndex();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronCellular::
Draw(void) const
{
   glBegin(GL_POINTS);
   /* TODO add color for cellular */
   if (read_voxel_count) {
      // draw supervoxels
      for (int iv = 0; iv < NVoxels(); ++iv) {
         R3LoadPoint(LocalVoxel(iv)->Position());
      }
   }
   else {
      // draw center of cellular
      R3LoadPoint(bbox.Centroid());
   }
   glEnd();
}



void NeuronCellular::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // check fp
   if (!fp) fp = stdout;

   // print cellular
   if (prefix) fprintf(fp, "%s", prefix);
   fprintf(fp, "Cellular %d:", data_index);
   if (suffix) fprintf(fp, "%s", suffix);
   fprintf(fp, "\n");
}