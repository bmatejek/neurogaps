// Source file for the neuron supervoxel class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(NeuronSupervoxel);



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronSupervoxel::
NeuronSupervoxel(void) :
file_offset(0),
voxels(NULL),
human_labels(),
predictions(),
boundaries()
{
}



NeuronSupervoxel::
~NeuronSupervoxel(void)
{
   if (voxels) delete[] voxels;
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

NeuronVoxel *NeuronSupervoxel::
GlobalVoxel(int global_index) const
{
   // just checking
   rn_assertion(voxels != NULL);
   rn_assertion(voxel_mapping != NULL);

   // do a binary search over voxel mapping
   return &voxels[LocalIndex(global_index)];
}



NeuronBoundary *NeuronSupervoxel::
Boundary(NeuronSupervoxel *neighbor) const
{
   rn_assertion(neighbor != NULL);

   // linear time search of boundaries
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      NeuronBoundary *candidate = boundaries[ib];
      if (candidate->OtherSupervoxel(this) == neighbor) return candidate;
   }

   // neighbor not found
   return NULL;
}





////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

int NeuronSupervoxel::
InsertHumanLabel(NeuronHumanLabel *human_label)
{
   // insert this human label into array
   for (int ih = 0; ih < NHumanLabels(); ++ih) {
      if (human_label == HumanLabel(ih)) return 0;
   }

   human_labels.Insert(human_label);

   // return success
   return 1;
}



void NeuronSupervoxel::
RemoveHumanLabel(NeuronHumanLabel *human_label)
{
   // remove this human label from array
   for (int ih = 0; ih < NHumanLabels(); ++ih) {
      if (human_label == HumanLabel(ih)) {
         human_labels.RemoveKth(ih);
         break;
      }
   }
}



int NeuronSupervoxel::
InsertPrediction(NeuronPrediction *prediction)
{
   // insert this prediction into array
   for (int ip = 0; ip < NPredictions(); ++ip) {
      if (prediction == Prediction(ip)) return 0;
   }

   predictions.Insert(prediction);

   // return success
   return 1;
}



void NeuronSupervoxel::
RemovePrediction(NeuronPrediction *prediction)
{
   // remove this prediction from array
   for (int ip = 0; ip < NPredictions(); ++ip) {
      if (prediction == Prediction(ip)) {
         predictions.RemoveKth(ip);
         break;
      }
   }
}



void NeuronSupervoxel::
InsertBoundary(NeuronBoundary *boundary)
{
   rn_assertion(!boundary->SupervoxelOne() || !boundary->SupervoxelTwo());

   // insert boundary
   if (!boundary->SupervoxelOne()) {
      boundary->supervoxel_one = this;
      boundary->supervoxel_one_index = boundaries.NEntries();
   }
   else {
      boundary->supervoxel_two = this;
      boundary->supervoxel_two_index = boundaries.NEntries();
   }
   boundaries.Insert(boundary);
}



void NeuronSupervoxel::
RemoveBoundary(NeuronBoundary *boundary)
{
   if (boundary->supervoxel_one == this) {
      rn_assertion(boundary->supervoxel_one_index >= 0);

      // remove boundary
      RNArrayEntry *entry = boundaries.KthEntry(boundary->supervoxel_one_index);
      NeuronBoundary *tail = boundaries.Tail();
      if (tail->supervoxel_one == this)
         tail->supervoxel_one_index = boundary->supervoxel_one_index;
      else
         tail->supervoxel_two_index = boundary->supervoxel_one_index;
      boundaries.EntryContents(entry) = tail;
      boundaries.RemoveTail();
      boundary->supervoxel_one_index = -1;
      boundary->supervoxel_one = NULL;
   }
   else if (boundary->supervoxel_two == this) {
      rn_assertion(boundary->supervoxel_two_index >= 0);

      // remove boundary
      RNArrayEntry *entry = boundaries.KthEntry(boundary->supervoxel_two_index);
      NeuronBoundary *tail = boundaries.Tail();
      if (tail->supervoxel_one == this)
         tail->supervoxel_one_index = boundary->supervoxel_two_index;
      else
         tail->supervoxel_two_index = boundary->supervoxel_two_index;
      boundaries.EntryContents(entry) = tail;
      boundaries.RemoveTail();
      boundary->supervoxel_two_index = -1;
      boundary->supervoxel_two = NULL;
   }
   else {
      /* should never reach here */
      rn_assertion(FALSE);
   }
}



void NeuronSupervoxel::
CreateVoxels(int nvoxels)
{
   rn_assertion(voxels == NULL);
   this->nvoxels = nvoxels;
   
   // create voxels
   this->voxels = new NeuronVoxel[nvoxels];
   rn_assertion(voxels != NULL);

   // create voxel mapping
   this->voxel_mapping = new int[nvoxels];
   rn_assertion(voxel_mapping != NULL);

   // should never remove this from memory
   read_voxel_count = 1;
}



////////////////////////////////////////////////////////////////////////
// Display functions
////////////////////////////////////////////////////////////////////////

void NeuronSupervoxel::
Draw(void) const
{
   // overridden
}



void NeuronSupervoxel::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
   // overridden
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void NeuronSupervoxel::
UpdateBBox(void)
{
   if (!voxels) { bbox = R3null_box; }
   else {
      rn_assertion(AreVoxelsResident());
      bbox = R3null_box;
      for (int iv = 0; iv < nvoxels; ++iv) {
         NeuronVoxel *voxel = LocalVoxel(iv);
         bbox.Union(voxel->Position());
      }
   }
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int NeuronSupervoxel::
ReadVoxels(void)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;
   if (voxels) delete[] voxels;
   if (voxel_mapping) delete[] voxel_mapping;

   // get the data filename
   const char *path = data->FilePath();
   const char *file = data->Filename();
   char filename[4096];
   sprintf(filename, "%s/%s", path, file);

   // open file
   FILE *fp = fopen(filename, "rb");
   if (!fp) {
      fprintf(stderr, "Unable to open Neuron data file: %s\n", filename);
      return 0;
   }

   // seek within file
   if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
      fprintf(stderr, "Unable to seek to volume position\n");
      return 0;
   }

   // read number of voxels
   if (fread(&nvoxels, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to read number of voxels\n");
      return 0;
   }

   // check number of voxels
   if (nvoxels > 0) {
      // allocate voxels
      voxels = new NeuronVoxel[nvoxels];
      if (!voxels) {
         fprintf(stderr, "Unable to allocate voxels\n");
         return 0;
      }
      voxel_mapping = new int[nvoxels];
      if (!voxel_mapping) {
         fprintf(stderr, "Unable to allocate voxels\n");
         return 0;
      }
      // read voxels
      for (int iv = 0; iv < nvoxels; ++iv) {
         voxels[iv].ReadVoxel(fp, data);
         SetVoxelMapping(iv, voxels[iv].DataIndex());
      }

      // update bounding box
      bbox = R3null_box;
      UpdateBBox();
   }
   else {
      // no voxels
      voxels = NULL;
      bbox = R3null_box;
   }

   // close the file
   fclose(fp);

   // return success
   return 1;
}



int NeuronSupervoxel::
ReadVoxels(FILE *fp, RNBoolean seek)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;
   rn_assertion(fp != NULL);
   if (voxels) delete[] voxels;
   if (voxel_mapping) delete[] voxel_mapping;

   // seek within file
   if (seek && (file_offset > 0)) {
      if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
         fprintf(stderr, "Unable to seek to volume position\n");
         return 0;
      }
   }
   else {
      // update file offset
      file_offset = RNFileTell(fp);
   }

   // read number of voxels
   if (fread(&nvoxels, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to read number of voxels\n");
      return 0;
   }

   // check number of voxels
   if (nvoxels > 0) {
      // allocate voxels
      voxels = new NeuronVoxel[nvoxels];
      if (!voxels) {
         fprintf(stderr, "Unable to allocate voxels\n");
         return 0;
      }
      voxel_mapping = new int[nvoxels];
      if (!voxel_mapping) {
         fprintf(stderr, "Unable to allocate voxels\n");
         return 0;
      }
      // read voxels
      for (int iv = 0; iv < nvoxels; ++iv) {
         voxels[iv].ReadVoxel(fp, data);
         SetVoxelMapping(iv, voxels[iv].DataIndex());
      }

      // update bounding box
      bbox = R3null_box;
      UpdateBBox();
   }
   else {
      // no voxels
      voxels = NULL;
      bbox = R3null_box;
   }

   // return success
   return 1;
}



int NeuronSupervoxel::
WriteVoxels(FILE *fp, RNBoolean seek)
{
   // just cheking
   rn_assertion(fp != NULL);
   rn_assertion(voxels != NULL);
   rn_assertion(voxel_mapping != NULL);

   // seek within file
   if (seek && (file_offset > 0)) {
      // seek to file offset
      if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
         fprintf(stderr, "Unable to seek to volume position\n");
         return 0;
      }
   }
   else {
      // update file offset
      file_offset = RNFileTell(fp);
   }

   // write number of voxels
   if (fwrite(&nvoxels, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to write number of voxels\n");
      return 0;
   }

   // write voxels
   if (nvoxels > 0) {
      for (int iv = 0; iv < nvoxels; ++iv) {
         if (!voxels[iv].WriteVoxel(fp)) {
            fprintf(stderr, "Unable to write voxels\n");
            return 0;
         }
      }
   }

   // return success
   return 1;
}



int NeuronSupervoxel::
ReleaseVoxels(void)
{
   // check read counter
   if (--read_voxel_count > 0) return 1;
   rn_assertion(voxels);

   // release voxels
   delete[] voxels;
   voxels = NULL;

   // delete voxel mapping
   delete[] voxel_mapping;
   voxel_mapping = NULL;

   // return success
   return 1;
}
