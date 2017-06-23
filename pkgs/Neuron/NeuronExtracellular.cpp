// Source file for the neuron extracellular class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronExtracellular);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronExtracellular::
NeuronExtracellular(void)
{
}



NeuronExtracellular::
~NeuronExtracellular(void)
{
   data->RemoveExtracellular(this);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int NeuronExtracellular::
SupervoxelIndex(void) const 
{
   if (!data) return -1;
   // return offset supervoxel index
   return data->NCellulars() + DataIndex();
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronExtracellular::
Draw(void) const
{
   glBegin(GL_POINTS);
   /* TODO add color for extracellular */
   if (read_voxel_count) {
      // draw supervoxels
      for (int iv = 0; iv < NVoxels(); ++iv) {
         R3LoadPoint(LocalVoxel(iv)->Position());
      }
   }
   else {
      // draw center of extracellular
      R3LoadPoint(bbox.Centroid());
   }
   glEnd();
}



void NeuronExtracellular::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // check fp
   if (!fp) fp = stdout;

   // print extracellular
   if (prefix) fprintf(fp, "%s", prefix);
   fprintf(fp, "Extracellular %d:", data_index);
   if (suffix) fprintf(fp, "%s", suffix);
   fprintf(fp, "\n");
}