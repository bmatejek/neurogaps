// Source file for the neuron human label class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronHumanLabel);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronHumanLabel::
NeuronHumanLabel(void)
{
}



NeuronHumanLabel::
~NeuronHumanLabel(void)
{
   data->RemoveHumanLabel(this);
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void NeuronHumanLabel::
InsertSupervoxel(NeuronSupervoxel *supervoxel)
{
   if (!supervoxel->InsertHumanLabel(this)) return;
   supervoxels.Insert(supervoxel);
}



void NeuronHumanLabel::
RemoveSupervoxel(NeuronSupervoxel *supervoxel)
{
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (supervoxel == Supervoxel(is)) {
         supervoxel->RemoveHumanLabel(this);
         supervoxels.RemoveKth(is);
         return;
      }
   }
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronHumanLabel::
Draw(void) const
{
   glBegin(GL_POINTS);
   /* TODO add color for human label */
   if (read_voxel_count) {
      // draw supervoxels
      for (int is = 0; is < NSupervoxels(); ++is) {
         Supervoxel(is)->Draw();
      }
   }
   else {
      // draw center of human label
      R3LoadPoint(bbox.Centroid());
   }
   glEnd();
}



void NeuronHumanLabel::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // check fp
   if (!fp) fp = stdout;

   // print human label
   if (prefix) fprintf(fp, "%s", prefix);
   fprintf(fp, "Human Label %d:", data_index);
   if (suffix) fprintf(fp, "%s", suffix);
   fprintf(fp, "\n");

   // add indentation to prefix
   char indented_prefix[1024];
   sprintf(indented_prefix, "%s  ", prefix);

   // print supervoxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      Supervoxel(is)->Print(fp, indented_prefix, suffix);
   }
}