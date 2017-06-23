// Source file for the neuron prediction class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronPrediction);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronPrediction::
NeuronPrediction(void) :
boundaries()
{
}



NeuronPrediction::
~NeuronPrediction(void)
{
   data->RemovePrediction(this);
}


////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

NeuronPredictionBoundary *NeuronPrediction::
Boundary(NeuronPrediction *neighbor) const
{
   rn_assertion(neighbor != NULL);

   // linear time search of boundaries
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      NeuronPredictionBoundary *candidate = boundaries[ib];
      if (candidate->OtherPrediction(this) == neighbor) return candidate;
   }

   // neighbor not found
   return NULL;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void NeuronPrediction::
InsertSupervoxel(NeuronSupervoxel *supervoxel)
{
   if (!supervoxel->InsertPrediction(this)) return;
   supervoxels.Insert(supervoxel);
}



void NeuronPrediction::
RemoveSupervoxel(NeuronSupervoxel *supervoxel)
{
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (supervoxel == Supervoxel(is)) {
         supervoxel->RemovePrediction(this);
         supervoxels.RemoveKth(is);
         return;
      }
   }
}



void NeuronPrediction::
InsertBoundary(NeuronPredictionBoundary *boundary)
{
   rn_assertion(!boundary->PredictionOne() || !boundary->PredictionTwo());

   // insert boundary
   if (!boundary->PredictionOne()) {
      boundary->prediction_one = this;
      boundary->prediction_one_index = boundaries.NEntries();
   }
   else {
      boundary->prediction_two = this;
      boundary->prediction_two_index = boundaries.NEntries();
   }
   boundaries.Insert(boundary);
}



void NeuronPrediction::
RemoveBoundary(NeuronPredictionBoundary *boundary)
{
   if (boundary->prediction_one == this) {
      rn_assertion(boundary->prediction_one_index >= 0);

      // remove boundary
      RNArrayEntry *entry = boundaries.KthEntry(boundary->prediction_one_index);
      NeuronPredictionBoundary *tail = boundaries.Tail();
      if (tail->prediction_one == this)
         tail->prediction_one_index = boundary->prediction_one_index;
      else
         tail->prediction_two_index = boundary->prediction_one_index;
      boundaries.EntryContents(entry) = tail;
      boundaries.RemoveTail();
      boundary->prediction_one_index = -1;
      boundary->prediction_one = NULL;
   }
   else if (boundary->prediction_two == this) {
      rn_assertion(boundary->prediction_two_index >= 0);

      // remove boundary
      RNArrayEntry *entry = boundaries.KthEntry(boundary->prediction_two_index);
      NeuronPredictionBoundary *tail = boundaries.Tail();
      if (tail->prediction_one == this)
         tail->prediction_one_index = boundary->prediction_two_index;
      else
         tail->prediction_two_index = boundary->prediction_two_index;
      boundaries.EntryContents(entry) = tail;
      boundaries.RemoveTail();
      boundary->prediction_two_index = -1;
      boundary->prediction_two = NULL;
   }
   else {
      /* should never reach here */
      rn_assertion(FALSE);
   }
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronPrediction::
Draw(void) const
{
   glBegin(GL_POINTS);
   /* TODO add color for prediction */
   if (read_voxel_count) {
      // draw supervoxels
      for (int is = 0; is < NSupervoxels(); ++is) {
         Supervoxel(is)->Draw();
      }
   }
   else {
      // draw center of prediction
      R3LoadPoint(bbox.Centroid());
   }
   glEnd();
}



void NeuronPrediction::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // check fp
   if (!fp) fp = stdout;

   // print prediction
   if (prefix) fprintf(fp, "%s", prefix);
   fprintf(fp, "Prediction %d:", data_index);
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