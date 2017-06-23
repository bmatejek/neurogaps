// Source file for the neuron boundary class

/////////////////////////////////////////////////////////////////////
// Include files
/////////////////////////////////////////////////////////////////////

#include "Neuron.h"
#include <algorithm>



// useful cache directories

static const char *features_directory = "postprocess";
static RNScalar *random_walk_features = NULL;
static char random_walk_features_data_name[4096];
static RNScalar *drunkard_walk_features = NULL;
static char drunkard_walk_features_data_name[4096];



/////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
/////////////////////////////////////////////////////////////////////

NeuronPredictionBoundary::
NeuronPredictionBoundary(void) :
data(NULL),
data_index(-1),
naffinities(0),
file_offset(0),
read_voxel_count(0),
read_affinities_count(0),
affinities(NULL),
prediction_one(NULL),
prediction_two(NULL),
prediction_one_index(-1),
prediction_two_index(-1),
maximum(RN_UNKNOWN),
mean(RN_UNKNOWN),
median(RN_UNKNOWN),
minimum(RN_UNKNOWN),
stddev(RN_UNKNOWN),
skew(RN_UNKNOWN),
user_data(NULL)
{
}



NeuronPredictionBoundary::
~NeuronPredictionBoundary(void)
{
   // delete the affinities array if it exists
   if (affinities) delete[] affinities;

   // remove the boundary from data
   if (data) data->RemovePredictionBoundary(this);
}



/////////////////////////////////////////////////////////////////////
// Manipulation functions
/////////////////////////////////////////////////////////////////////

void NeuronPredictionBoundary::
CreateNAffinities(int naffinities)
{
   rn_assertion(affinities == NULL);
   // create affinities matrix
   this->naffinities = naffinities;
   this->affinities = new RNScalar[naffinities];
   read_affinities_count++;
}



void NeuronPredictionBoundary::
UpdateVoxelAffinity(RNScalar affinity, int index)
{
   rn_assertion(affinities != NULL);
   rn_assertion((0 <= index) && (index < naffinities));
   // add the voxel affinity
   affinities[index] = affinity;
}



void NeuronPredictionBoundary::
UpdateStatistics(void)
{
   // read affinities if not in memory
   RNBoolean read_affinities = FALSE;
   if (!affinities) { ReadAffinities(); read_affinities = TRUE; }

   // sort the affinities
   std::sort(affinities, affinities + naffinities);

   // update all of the statistics on the boundary
   mean = 0.0;
   minimum = affinities[0];
   maximum = affinities[naffinities - 1];
   for (int ia = 0; ia < naffinities; ++ia) {
      mean += affinities[ia];
   }
   mean /= naffinities;

   // update the median value
   if (naffinities % 2 == 0) {
      median = 0.5 * (affinities[naffinities / 2 - 1] + affinities[naffinities / 2]);
   }
   else {
      median = affinities[naffinities / 2];
   }

   // calculate the standard deviation
   stddev = 0.0;
   for (int ia = 0; ia < naffinities; ++ia) {
      stddev += (affinities[ia] - mean) * (affinities[ia] - mean);
   }
   if (naffinities > 1)
      stddev = sqrt(stddev / (naffinities - 1));
   else
      stddev = 0.0;

   // calculate the skew of the affinities
   skew = 0.0;
   for (int ia = 0; ia < naffinities; ++ia) {
      skew += (affinities[ia] - mean) * (affinities[ia] - mean) * (affinities[ia] - mean);
   }
   skew /= naffinities;
   if (stddev > 0.0)
      skew /= (stddev * stddev * stddev);
   else
      skew = 0.0;

   // release affinities if previous read
   if (read_affinities) ReleaseAffinities();
}



/////////////////////////////////////////////////////////////////////
// Property functions
/////////////////////////////////////////////////////////////////////

RNScalar NeuronPredictionBoundary::
RandomWalk(RNBoolean drunkards_walk) const
{
   // get root filename
   char root_filename[4096];
   strncpy(root_filename, data->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   if (drunkards_walk && drunkard_walk_features) {
      if (!strcmp(drunkard_walk_features_data_name, root_filename)) {
         return drunkard_walk_features[this->DataIndex()];
      }
      else {
         delete[] drunkard_walk_features;
      }
   }
   else if (!drunkards_walk && random_walk_features) {
      if (!strcmp(random_walk_features_data_name, root_filename)) {
         return random_walk_features[this->DataIndex()];
      }
      else {
         delete[] random_walk_features;
      }
   }

   // see if the cache file exists
   char random_walk_cache_filename[4096];
   if (drunkards_walk)
      sprintf(random_walk_cache_filename, "%s/%s_drunkards_walk.feature", features_directory, root_filename);
   else
      sprintf(random_walk_cache_filename, "%s/%s_random_walk.feature", features_directory, root_filename);

   // open up the cache file if it exists
   FILE *cache_fp = fopen(random_walk_cache_filename, "rb");
   if (cache_fp) {
      int nboundaries;
      fread(&nboundaries, sizeof(int), 1, cache_fp);
      rn_assertion(nboundaries == data->NPredictionBoundaries());

      RNScalar value;
      if (drunkards_walk) {
         drunkard_walk_features = new RNScalar[data->NPredictionBoundaries()];
         strncpy(drunkard_walk_features_data_name, root_filename, 4096);

         // read into memory
         fread(drunkard_walk_features, sizeof(RNScalar), data->NPredictionBoundaries(), cache_fp);

         value = drunkard_walk_features[this->DataIndex()];
      }
      else {
         random_walk_features = new RNScalar[data->NPredictionBoundaries()];
         strncpy(random_walk_features_data_name, root_filename, 4096);

         // read into memory
         fread(random_walk_features, sizeof(RNScalar), data->NPredictionBoundaries(), cache_fp);

         value = random_walk_features[this->DataIndex()];
      }

      // close the cache
      fclose(cache_fp);

      // return the value of the drunkards walk
      return value;
   }

   // make sure voxels are resident
   rn_assertion(data->AreVoxelsResident());

   // get two voxels within the supervoxels
   NeuronVoxel *voxel_one = prediction_one->CenterVoxel();
   NeuronVoxel *voxel_two = prediction_two->CenterVoxel();

   // coordinates of the goal voxels
   int goalx = voxel_two->XCoordinate();
   int goaly = voxel_two->YCoordinate();
   int goalz = voxel_two->ZCoordinate();

   // number of random iterations
   int niterations = 250;
   RNScalar max_step_sum = 0.0;

   for (int it = 0; it < niterations; ++it) {
      // run a randomized walk between voxel one and voxel two
      int ix = voxel_one->XCoordinate();
      int iy = voxel_one->YCoordinate();
      int iz = voxel_one->ZCoordinate();

      RNScalar max_step = 0.0;
      while (ix != goalx || iy != goaly || iz != goalz) {
         // get the current voxel
         NeuronVoxel *voxel = data->Voxel(ix, iy, iz);

         // get distance left to travel
         int dx = goalx - ix;
         int dy = goaly - iy;
         int dz = goalz - iz;

         RNScalar xdirectional_value = 0.0;
         if (dx > 0) { xdirectional_value = voxel->AffinityToNeighbor(voxel->NextVoxel(RN_X)); }
         else if (dx < 0) { xdirectional_value = voxel->AffinityToNeighbor(voxel->PreviousVoxel(RN_X)); }

         RNScalar ydirectional_value = 0.0;
         if (dy > 0) { ydirectional_value = voxel->AffinityToNeighbor(voxel->NextVoxel(RN_Y)); }
         else if (dy < 0) { ydirectional_value = voxel->AffinityToNeighbor(voxel->PreviousVoxel(RN_Y)); }

         RNScalar zdirectional_value = 0.0;
         if (dz > 0) { zdirectional_value = voxel->AffinityToNeighbor(voxel->NextVoxel(RN_Z)); }
         else if (dz < 0) { zdirectional_value = voxel->AffinityToNeighbor(voxel->PreviousVoxel(RN_Z)); }

         // get a random number
         RNScalar random_number = RNRandomScalar();
         int direction_of_motion;

         // drunkard prefers no direction
         if (drunkards_walk) {
            RNScalar sum_values = (dx != 0) + (dy != 0) + (dz != 0);
            RNScalar probabilities[3] = { (dx != 0) / sum_values, (dy != 0) / sum_values, (dz != 0) / sum_values };

            // get direction of random walk
            if (random_number < probabilities[RN_X]) direction_of_motion = RN_X;
            else if (random_number < probabilities[RN_X] + probabilities[RN_Y]) direction_of_motion = RN_Y;
            else direction_of_motion = RN_Z;
         }
         // favor areas of high affinity
         else {
            // find the probability of moving in each direction
            RNScalar sum_values = xdirectional_value + ydirectional_value + zdirectional_value;
            RNScalar probabilities[3] = { xdirectional_value / sum_values, ydirectional_value / sum_values, zdirectional_value / sum_values };

            // get direction of random walk
            if (random_number < probabilities[RN_X]) direction_of_motion = RN_X;
            else if (random_number < probabilities[RN_X] + probabilities[RN_Y]) direction_of_motion = RN_Y;
            else direction_of_motion = RN_Z;
         }

         if (direction_of_motion == RN_X) {
            if (dx < 0) --ix;
            else ++ix;
            RNScalar tmp_distance = (1.0 - xdirectional_value);
            if (tmp_distance > max_step) max_step = tmp_distance;
         }
         else if (direction_of_motion == RN_Y) {
            if (dy < 0) --iy;
            else ++iy;
            RNScalar tmp_distance = (1.0 - ydirectional_value);
            if (tmp_distance > max_step) max_step = tmp_distance;
         }
         else {
            if (dz < 0) --iz;
            else ++iz;
            RNScalar tmp_distance = (1.0 - zdirectional_value);
            if (tmp_distance > max_step) max_step = tmp_distance;
         }
      }
      max_step_sum += max_step;
   }

   return max_step_sum / niterations;
}



RNScalar NeuronPredictionBoundary::
BoundaryRank(void) const
{
   int greater_boundaries = 0;
   for (int ib = 0; ib < prediction_one->NBoundaries(); ++ib) {
      NeuronPredictionBoundary *boundary = prediction_one->Boundary(ib);
      if (boundary == this) continue;

      RNScalar boundary_mean = boundary->Mean();
      if (boundary_mean > mean) greater_boundaries++;
   }

   for (int ib = 0; ib < prediction_two->NBoundaries(); ++ib) {
      NeuronPredictionBoundary *boundary = prediction_two->Boundary(ib);
      if (boundary == this) continue;

      RNScalar boundary_mean = boundary->Mean();
      if (boundary_mean > mean) greater_boundaries++;
   }

   // return the boundary ranking
   return (RNScalar)(greater_boundaries + 1);
}



RNScalar NeuronPredictionBoundary::
BoundaryScaledRanking(void) const
{
   int greater_boundaries = 0;
   for (int ib = 0; ib < prediction_one->NBoundaries(); ++ib) {
      NeuronPredictionBoundary *boundary = prediction_one->Boundary(ib);
      if (boundary == this) continue;

      RNScalar boundary_mean = boundary->Mean();
      if (boundary_mean > mean) greater_boundaries++;
   }

   for (int ib = 0; ib < prediction_two->NBoundaries(); ++ib) {
      NeuronPredictionBoundary *boundary = prediction_two->Boundary(ib);
      if (boundary == this) continue;

      RNScalar boundary_mean = boundary->Mean();
      if (boundary_mean > mean) greater_boundaries++;
   }

   int ncompeting_boundaries = prediction_one->NBoundaries() + prediction_two->NBoundaries() - 1;

   return (ncompeting_boundaries - greater_boundaries) / (RNScalar)ncompeting_boundaries;
}



RNScalar NeuronPredictionBoundary::
CombinedVoxels(void) const
{
   // return the number of voxels combined between the two volumes
   return (RNScalar)(prediction_one->NVoxels() + prediction_two->NVoxels());
}



RNScalar NeuronPredictionBoundary::
CombinedVoxelsProportion(void) const
{
   // return the proportion of voxels in the smaller volume
   if (prediction_one->NVoxels() < prediction_two->NVoxels())
      return prediction_one->NVoxels() / (RNScalar)(prediction_one->NVoxels() + prediction_two->NVoxels());
   else
      return prediction_two->NVoxels() / (RNScalar)(prediction_one->NVoxels() + prediction_two->NVoxels());
}



////////////////////////////////////////////////////////////////////////
// Input/Output functions
////////////////////////////////////////////////////////////////////////

int NeuronPredictionBoundary::
ReadAffinities(void)
{
   // check read counter
   if (++read_affinities_count > 1) return 1;
   rn_assertion(affinities == NULL);

   // get the data filename
   const char *filename = data->Filename();
   if (!filename) return 0;

   // open file
   FILE *fp = fopen(filename, "rb");
   if (!fp) {
      fprintf(stderr, "Unable to open Neuron data file: %s\n", filename);
      return 0;
   }

   // seek within file
   if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
      fprintf(stderr, "Unable to seek to affinity position\n");
      return 0;
   }

   // read number of affinities
   if (fread(&naffinities, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to read number of affinities\n");
      return 0;
   }

   // allocate affinities
   if (naffinities > 0) {
      // allocate affinities
      affinities = new RNScalar[naffinities];
      if (!affinities) {
         fprintf(stderr, "Unable to allocate affinities\n");
         return 0;
      }
      // read affinities
      if (fread(&affinities[0], sizeof(RNScalar), naffinities, fp) != (unsigned int)naffinities) {
         fprintf(stderr, "Unable to read affinities\n");
         return 0;
      }
   }
   else {
      // no affinitites
      affinities = NULL;
   }

   // close the file
   fclose(fp);

   // return success
   return 1;
}



int NeuronPredictionBoundary::
ReadAffinities(FILE *fp, RNBoolean seek)
{
   // check read counter
   if (++read_affinities_count > 1) return 1;
   rn_assertion(affinities == NULL);
   rn_assertion(fp != NULL);

   // seek within file
   if (seek && (file_offset > 0)) {
      if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
         fprintf(stderr, "Unable to seek to boundary position\n");
         return 0;
      }
   }
   else {
      // update file offset
      file_offset = RNFileTell(fp);
   }

   // read number of affinities
   if (fread(&naffinities, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to read number of affinities\n");
      return 0;
   }

   // check number of affinities
   if (naffinities > 0) {
      // allocate affinities
      affinities = new RNScalar[naffinities];
      if (!affinities) {
         fprintf(stderr, "Unable to allocate affinities\n");
         return 0;
      }
      // read affinities
      if (fread(&affinities[0], sizeof(RNScalar), naffinities, fp) != (unsigned int)naffinities) {
         fprintf(stderr, "Unable to read affinities\n");
         return 0;
      }
   }
   else {
      // no affinities
      affinities = NULL;
   }

   // return success
   return 1;
}



int NeuronPredictionBoundary::
WriteAffinities(FILE *fp, RNBoolean seek)
{
   // just checking
   rn_assertion(fp != NULL);

   // seek within file
   if (seek && (file_offset > 0)) {
      // seek to file offset
      if (!RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET)) {
         fprintf(stderr, "Unable to seek to boundary position\n");
         return 0;
      }
   }
   else {
      // update file offset
      file_offset = RNFileTell(fp);
   }

   // write number of affinities
   if (fwrite(&naffinities, sizeof(int), 1, fp) != (unsigned int)1) {
      fprintf(stderr, "Unable to write number of affinities\n");
      return 0;
   }

   // write affinities
   if (naffinities > 0) {
      if (fwrite(&affinities[0], sizeof(RNScalar), naffinities, fp) != (unsigned int)naffinities) {
         fprintf(stderr, "Unable to write affinities\n");
         return 0;
      }
   }

   // return success
   return 1;
}



int NeuronPredictionBoundary::
ReleaseAffinities(void)
{
   // check read counter
   if (--read_affinities_count > 0) return 1;
   rn_assertion(affinities != NULL);

   // release  affinities
   delete[] affinities;
   affinities = NULL;

   // return success
   return 1;
}



int NeuronPredictionBoundary::
ReadVoxels(void)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;

   if (!prediction_one->ReadVoxels()) return 0;
   if (!prediction_two->ReadVoxels()) return 0;

   // return success
   return 1;
}



int NeuronPredictionBoundary::
ReadVoxels(FILE *fp, RNBoolean seek)
{
   // just checking
   rn_assertion(fp != NULL);

   // check read counter
   if (++read_voxel_count > 1) return 1;

   if (!prediction_one->ReadVoxels(fp, seek)) return 0;
   if (!prediction_two->ReadVoxels(fp, seek)) return 0;

   // return success
   return 1;
}



int NeuronPredictionBoundary::
ReleaseVoxels(void)
{
   // check read counter
   if (--read_voxel_count > 0) return 1;

   if (!prediction_one->ReleaseVoxels()) return 0;
   if (!prediction_two->ReleaseVoxels()) return 0;

   // return success
   return 1;
}
