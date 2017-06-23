// Source file for the neuron data class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "Neuron.h"
#include <vector>
#include <stdio.h>
#include <algorithm>



////////////////////////////////////////////////////////////////////////
// Useful constants
////////////////////////////////////////////////////////////////////////

#define HEADER_FILLER 12
#define DATA_FILLER 12
#define BOUNDARY_FILLER 10
#define CELLULAR_FILLER 4
#define EXTRACELLULAR_FILLER 4
#define PREDICTION_FILLER 4
#define HUMAN_LABEL_FILLER 4
#define NTEST_METRICS 12
#define NSEGMENTATION_METRICS 8

enum { MAXIMUM_METRIC, MEAN_METRIC, MEDIAN_METRIC, MINIMUM_METRIC };



static const char *test_metric_names[NTEST_METRICS] = {
   "Rand Error Full", "Rand Error Merge", "Rand Error Split",
   "Rand F-Score Full", "Rand F-Score Merge", "Rand F-Score Split", 
   "Variation of F-Score Full", "Variation F-Score Merge", "Variation F-Score Split", 
   "Variation of Information Full", "Variation of Information Merge", "Variation of Information Split"
};



static const char *segmentation_metric_names[NSEGMENTATION_METRICS] = {
   "Correct Merges", "Bad Merges", "Bad Splits", "Good Splits",
   "Merge Precision", "Merge Recall", "Split Precision", "Split Recall"
};



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

NeuronData::
NeuronData(void) :
human_labels(),
predictions(),
cellulars(),
extracellulars(),
prediction_boundaries(),
boundaries(),
read_voxel_count(0),
voxel_mapping(NULL),
transformation(R3null_affine),
data_size(-1),
data_sheet_size(-1),
data_row_size(-1),
bbox(R3null_box),
user_data(NULL)
{
   // set strings to NULL
   file_path[0] = '\0';
   filename[0] = '\0';
   affinities_filename[0] = '\0';
   human_labels_filename[0] = '\0';
   image_filename[0] = '\0';
   machine_labels_filename[0] = '\0';

   // set dimensions of resolution and scaling
   for (int dim = 0; dim <= 2; ++dim) {
      resolution[dim] = -1;
      scaling[dim] = -1;
   }
}



NeuronData::
~NeuronData(void)
{
   if (voxel_mapping) delete[] voxel_mapping;

   // delete all volumes
   while (NCellulars() > 0)
      delete Cellular(0);
   while (NExtracellulars() > 0)
      delete Extracellular(0);
   while (NHumanLabels() > 0)
      delete HumanLabel(0);
   while (NPredictions() > 0)
      delete Prediction(0);
   while (NBoundaries() > 0)
      delete Boundary(0);
}



////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
NVoxels(void) const
{
   return data_size;
}



NeuronVoxel *NeuronData::
Voxel(int ix, int iy, int iz) const
{
   // return voxel at (ix, iy, iz)
   return Voxel(IndicesToIndex(ix, iy, iz));
}



NeuronVoxel *NeuronData::
Voxel(int voxel_index) const
{
   rn_assertion((0 <= voxel_index) && (voxel_index < data_size));
   rn_assertion(voxel_mapping);

   // get voxel from supervoxel
   int supervoxel_index = voxel_mapping[voxel_index];
   NeuronSupervoxel *supervoxel = Supervoxel(supervoxel_index);
   if (supervoxel->AreVoxelsResident()) return supervoxel->GlobalVoxel(voxel_index);
   else return NULL;
}



int NeuronData::
NTestMetrics(void) const
{
   // return the number of test metrics
   return NTEST_METRICS;
}



const char *NeuronData::
TestMetricName(int index) const
{
   rn_assertion((0 <= index) && (index < NTEST_METRICS));
   // return this metric name
   return test_metric_names[index];
}



int NeuronData::
NSegmentationMetrics(void) const
{
   // return the number of segmentation metrics
   return NSEGMENTATION_METRICS;
}



const char *NeuronData::
SegmentationMetricName(int index) const
{
   rn_assertion((0 <= index) && (index < NSEGMENTATION_METRICS));
   // return this metric name
   return segmentation_metric_names[index];
}



////////////////////////////////////////////////////////////////////////
// Property functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
AffinityStatistics(RNScalar(&mean)[4], RNScalar(&stddev)[4]) const
{
   // make sure voxels are resident
   rn_assertion(AreVoxelsResident());

   // set base mean
   const int RN_ALL = 3;
   mean[RN_X] = 0.0;
   mean[RN_Y] = 0.0;
   mean[RN_Z] = 0.0;
   mean[RN_ALL] = 0.0;

   // go through each dimension and voxel
   for (int ix = 0; ix < resolution[RN_X]; ++ix) {
      for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
         for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);
            for (int dim = 0; dim <= 2; ++dim) {
               // add to mean variables            
               mean[dim] += voxel->affinities[dim];
               mean[RN_ALL] += voxel->affinities[dim];
            }
         }
      }
   }

   // divide sums
   mean[RN_X] /= NVoxels();
   mean[RN_Y] /= NVoxels();
   mean[RN_Z] /= NVoxels();
   mean[RN_ALL] /= (RN_ALL * NVoxels());

   // go through every voxel to calculate standard deviation
   for (int ix = 0; ix < resolution[RN_X]; ++ix) {
      for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
         for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);
            for (int dim = 0; dim <= 2; ++dim) {
               // add to stddev variables
               stddev[dim] += (voxel->affinities[dim] - mean[dim]) * (voxel->affinities[dim] - mean[dim]);
               stddev[RN_ALL] += (voxel->affinities[dim] - mean[RN_ALL]) * (voxel->affinities[dim] - mean[RN_ALL]);
            }
         }
      }
   }

   // divide sums
   stddev[RN_X] = sqrt(stddev[RN_X] / NVoxels());
   stddev[RN_Y] = sqrt(stddev[RN_Y] / NVoxels());
   stddev[RN_Z] = sqrt(stddev[RN_Z] / NVoxels());
   stddev[RN_ALL] = sqrt(stddev[RN_ALL] / (RN_ALL * NVoxels()));

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void NeuronData::
InsertCellular(NeuronCellular *cellular)
{
   // just checking
   rn_assertion(cellular->data_index == -1);
   rn_assertion(cellular->data == NULL);

   // insert supervoxel
   cellular->data = this;
   cellular->data_index = cellulars.NEntries();
   cellulars.Insert(cellular);
}



void NeuronData::
RemoveCellular(NeuronCellular *cellular)
{
   // just checking
   rn_assertion(cellular->data_index >= 0);
   rn_assertion(cellular->data_index < cellulars.NEntries());
   rn_assertion(cellular->data == this);

   // remove supervoxel
   RNArrayEntry *entry = cellulars.KthEntry(cellular->data_index);
   NeuronCellular *tail = cellulars.Tail();
   tail->data_index = cellular->data_index;
   cellulars.EntryContents(entry) = tail;
   cellulars.RemoveTail();
   cellular->data_index = -1;
   cellular->data = NULL;
}



void NeuronData::
InsertExtracellular(NeuronExtracellular *extracellular)
{
   // just checking
   rn_assertion(extracellular->data_index == -1);
   rn_assertion(extracellular->data == NULL);

   // insert extracellular
   extracellular->data = this;
   extracellular->data_index = extracellulars.NEntries();
   extracellulars.Insert(extracellular);
}



void NeuronData::
RemoveExtracellular(NeuronExtracellular *extracellular)
{
   // just checking
   rn_assertion(extracellular->data_index >= 0);
   rn_assertion(extracellular->data_index < extracellulars.NEntries());
   rn_assertion(extracellular->data == this);

   // remove extracellular
   RNArrayEntry *entry = extracellulars.KthEntry(extracellular->data_index);
   NeuronExtracellular *tail = extracellulars.Tail();
   tail->data_index = extracellular->data_index;
   extracellulars.EntryContents(entry) = tail;
   extracellulars.RemoveTail();
   extracellular->data_index = -1;
   extracellular->data = NULL;
}



void NeuronData::
InsertHumanLabel(NeuronHumanLabel *human_label)
{
   // just checking
   rn_assertion(human_label->data_index == -1);
   rn_assertion(human_label->data == NULL);

   // insert truth
   human_label->data = this;
   human_label->data_index = human_labels.NEntries();
   human_labels.Insert(human_label);
}



void NeuronData::
RemoveHumanLabel(NeuronHumanLabel *human_label)
{
   // just checking
   rn_assertion(human_label->data_index >= 0);
   rn_assertion(human_label->data_index < human_labels.NEntries());
   rn_assertion(human_label->data == this);

   // remove truth
   RNArrayEntry *entry = human_labels.KthEntry(human_label->data_index);
   NeuronHumanLabel *tail = human_labels.Tail();
   tail->data_index = human_label->data_index;
   human_labels.EntryContents(entry) = tail;
   human_labels.RemoveTail();
   human_label->data_index = -1;
   human_label->data = NULL;
}



void NeuronData::
InsertPrediction(NeuronPrediction *prediction)
{
   // just checking
   rn_assertion(prediction->data_index == -1);
   rn_assertion(prediction->data == NULL);

   // insert truth
   prediction->data = this;
   prediction->data_index = predictions.NEntries();
   predictions.Insert(prediction);
}



void NeuronData::
RemovePrediction(NeuronPrediction *prediction)
{
   // just checking
   rn_assertion(prediction->data_index >= 0);
   rn_assertion(prediction->data_index < predictions.NEntries());
   rn_assertion(prediction->data == this);

   // remove truth
   RNArrayEntry *entry = predictions.KthEntry(prediction->data_index);
   NeuronPrediction *tail = predictions.Tail();
   tail->data_index = prediction->data_index;
   predictions.EntryContents(entry) = tail;
   predictions.RemoveTail();
   prediction->data_index = -1;
   prediction->data = NULL;
}



void NeuronData::
InsertBoundary(NeuronBoundary *boundary)
{
   // just checking
   rn_assertion(boundary->data_index == -1);
   rn_assertion(boundary->data == NULL);

   // insert truth
   boundary->data = this;
   boundary->data_index = boundaries.NEntries();
   boundaries.Insert(boundary);
}



void NeuronData::
RemoveBoundary(NeuronBoundary *boundary)
{
   // just checking
   rn_assertion(boundary->data_index >= 0);
   rn_assertion(boundary->data_index < boundaries.NEntries());
   rn_assertion(boundary->data == this);

   // remove boundary
   RNArrayEntry *entry = boundaries.KthEntry(boundary->data_index);
   NeuronBoundary *tail = boundaries.Tail();
   tail->data_index = boundary->data_index;
   boundaries.EntryContents(entry) = tail;
   boundaries.RemoveTail();
   boundary->data_index = -1;
   boundary->data = NULL;
}



void NeuronData::
InsertPredictionBoundary(NeuronPredictionBoundary *prediction_boundary)
{
   // just checking
   rn_assertion(prediction_boundary->data_index == -1);
   rn_assertion(prediction_boundary->data == NULL);

   // insert truth
   prediction_boundary->data = this;
   prediction_boundary->data_index = prediction_boundaries.NEntries();
   prediction_boundaries.Insert(prediction_boundary);
}



void NeuronData::
RemovePredictionBoundary(NeuronPredictionBoundary *prediction_boundary)
{
   // just checking
   rn_assertion(prediction_boundary->data_index >= 0);
   rn_assertion(prediction_boundary->data_index < prediction_boundaries.NEntries());
   rn_assertion(prediction_boundary->data == this);

   // remove boundary
   RNArrayEntry *entry = prediction_boundaries.KthEntry(prediction_boundary->data_index);
   NeuronPredictionBoundary *tail = prediction_boundaries.Tail();
   tail->data_index = prediction_boundary->data_index;
   prediction_boundaries.EntryContents(entry) = tail;
   prediction_boundaries.RemoveTail();
   prediction_boundary->data_index = -1;
   prediction_boundary->data = NULL;
}



int NeuronData::
NormalizeAffinities(void)
{
   // create mean, stddev and define
   RNScalar mean[4];
   RNScalar stddev[4];
   const int ALL = 3;

   // get affinity statistics
   this->AffinityStatistics(mean, stddev);

   // make sure that voxels are resident
   rn_assertion(AreVoxelsResident());

   for (int ix = 0; ix < resolution[RN_X]; ++ix) {
      for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
         for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);
            for (int dim = 0; dim <= 2; ++dim) {
               voxel->affinities[dim] = ((voxel->affinities[dim] - mean[dim]) / stddev[dim]) * stddev[ALL] + mean[ALL];
            }
         }
      }
   }

   // return success
   return 1;
}



int NeuronData::
FilterBoundaryAffinities(int metric)
{
   // go through every voxel and its neighbors
   for (int iv = 0; iv < NVoxels(); ++iv) {
      NeuronVoxel *voxel = Voxel(iv);
      NeuronSupervoxel *supervoxel = voxel->Supervoxel();
      for (int in = 0; in < voxel->NNeighbors() / 2; ++in) {
         // get the neighbor volume
         NeuronVoxel *neighbor_voxel = voxel->NextVoxel(in);
         if (!neighbor_voxel) continue;

         NeuronSupervoxel *neighbor_supervoxel = neighbor_voxel->Supervoxel();
         if (supervoxel == neighbor_supervoxel) continue;

         // find this boundary
         NeuronBoundary *boundary = supervoxel->Boundary(neighbor_supervoxel);
         rn_assertion(boundary != NULL);
         RNScalar affinity;
         if (metric == MAXIMUM_METRIC) affinity = boundary->Maximum();
         else if (metric == MEAN_METRIC) affinity = boundary->Mean();
         else if (metric == MEDIAN_METRIC) affinity = boundary->Median();
         else if (metric == MINIMUM_METRIC) affinity = boundary->Minimum();
         else { fprintf(stderr, "Unknown metric: %d\n", metric); return 0; }

         // update the affinity
         voxel->affinities[in] = affinity;
      }
   }

   // return success
   return 1;
}



int NeuronData::
FilterAffinity(const char *filter_type, const char *extension)
{
   rn_assertion(AreVoxelsResident());

   // run an affinity filtering algorithm
   if (!strcmp(filter_type, "maximum")) { if (!FilterBoundaryAffinities(MAXIMUM_METRIC)) return 0; }
   else if (!strcmp(filter_type, "mean")) { if (!FilterBoundaryAffinities(MEAN_METRIC)) return 0; }
   else if (!strcmp(filter_type, "median")) { if (!FilterBoundaryAffinities(MEDIAN_METRIC)) return 0; }
   else if (!strcmp(filter_type, "minimum")) { if (!FilterBoundaryAffinities(MINIMUM_METRIC)) return 0; }
   else if (!strcmp(filter_type, "normalize")) { if (!NormalizeAffinities()) return 0; }
   else { fprintf(stderr, "Unrecognized filter type: %s\n", filter_type); return 0; }

   // create output filename
   char root_filename[4096];
   strncpy(root_filename, filename, 4096);
   char *endp = strchr(root_filename, '.');
   *endp = '\0';

   // create output filename
   char output_filename[4096];

   if (!strcmp(filter_type, "normalize")) sprintf(output_filename, "%s_normalize", root_filename);
   else if (!strcmp(filter_type, "minimum")) sprintf(output_filename, "%s_minimum", root_filename);
   else if (!strcmp(filter_type, "maximum")) sprintf(output_filename, "%s_maximum", root_filename);
   else if (!strcmp(filter_type, "mean")) sprintf(output_filename, "%s_mean", root_filename);
   else if (!strcmp(filter_type, "median")) sprintf(output_filename, "%s_median", root_filename);
   else if (!strcmp(filter_type, "gaussian")) sprintf(output_filename, "%s_gaussian", root_filename);
   else if (!strcmp(filter_type, "diffusion")) sprintf(output_filename, "%s_diffusion", root_filename);
   else { fprintf(stderr, "Unrecognized filter type: %s\n", filter_type); return 0; }

   // write the file
   if (!strncmp(extension, ".txt", 4)) {
      // create new affinity root filename
      char affinity_root[4096];
      strncpy(affinity_root, affinities_filename, 4096);

      // output depends on name
      if (!strcmp(filter_type, "normalize")) sprintf(affinity_root, "%s_normalize", affinity_root);
      else if (!strcmp(filter_type, "minimum")) sprintf(affinity_root, "%s_minimum", affinity_root);
      else if (!strcmp(filter_type, "maximum")) sprintf(affinity_root, "%s_maximum", affinity_root);
      else if (!strcmp(filter_type, "mean")) sprintf(affinity_root, "%s_mean", affinity_root);
      else if (!strcmp(filter_type, "median")) sprintf(affinity_root, "%s_median", affinity_root);
      else if (!strcmp(filter_type, "gaussian")) sprintf(affinity_root, "%s_gaussian", affinity_root);
      else if (!strcmp(filter_type, "diffusion")) sprintf(affinity_root, "%s_diffusion", affinity_root);
      else { fprintf(stderr, "Unrecognized filter type: %s\n", filter_type); return 0; }

      char affinities_root_filename[4096];
      sprintf(affinities_root_filename, "%s/%s", file_path, affinity_root);

      // update affinity filename
      sprintf(affinities_filename, "%s", affinity_root);

      R3Grid *affinities[3] = { NULL, NULL, NULL };
      for (int ia = 0; ia < 3; ++ia) {
         affinities[ia] = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
      }

      for (int iv = 0; iv < NVoxels(); ++iv) {
         NeuronVoxel *voxel = Voxel(iv);
         int ix, iy, iz;
         IndexToIndices(iv, ix, iy, iz);

         for (int dim = 0; dim <= 2; ++dim) {
            affinities[dim]->SetGridValue(ix, iy, iz, voxel->affinities[dim]);
         }
      }

      // save the file
      RNMeta affinities_meta = RNMeta("Float32", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 3);

      if (!RNWriteNeuronMetaRawFile(affinities_root_filename, affinities_meta, affinities)) return 0;

      // delete affinities grid
      for (int ia = 0; ia < 3; ++ia) {
         delete affinities[ia];
      }

      char output_txt_filename[4096];
      sprintf(output_txt_filename, "%s/%s.txt", file_path, output_filename);

      // open file
      FILE *fp = fopen(output_txt_filename, "w");
      if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_txt_filename); return 0; }

      fprintf(fp, "%s.meta\n", affinities_filename);
      fprintf(fp, "%s.meta\n", human_labels_filename);
      fprintf(fp, "%s.meta\n", image_filename);
      fprintf(fp, "%s.meta\n", machine_labels_filename);
      fprintf(fp, "Scale(%f,%f,%f)", scaling[RN_X], scaling[RN_Y], scaling[RN_Z]);

      // close file
      fclose(fp);
   }
   else { fprintf(stderr, "Unrecognized extension: %s\n", extension); return 0; }

   // update filename
   sprintf(filename, "%s%s", output_filename, extension);

   // return success
   return 1;
}



int NeuronData::
TestMetric(int *voxel_proposals, RNScalar *metrics, int print_verbose) const
{
   rn_assertion(voxel_proposals != NULL);
   rn_assertion(metrics != NULL);

   // find the number of proposals
   int nproposals = 0;
   int nvoxels = NVoxels();
   for (int iv = 0; iv < nvoxels; ++iv)
      if (voxel_proposals[iv] > nproposals) nproposals = voxel_proposals[iv];

   // get an array of voxel truths
   int nhumanlabels = NHumanLabels();
   int *voxel_truths = new int[nvoxels];
   if (print_verbose) { printf("Creating voxel truth mapping...\n  "); fflush(stdout); }
   for (int iv = 0; iv < nvoxels; ++iv) {
      if (print_verbose) RNProgressBar(iv, nvoxels);
      NeuronVoxel *voxel = Voxel(iv);
      if (voxel->HumanLabel()) voxel_truths[iv] = voxel->HumanLabel()->DataIndex() + 1;
      else voxel_truths[iv] = 0;
   }
   if (print_verbose) printf("\ndone!\n");

   // size of each segment
   unsigned long long *s = new unsigned long long[nproposals + 1];
   for (int ip = 0; ip < nproposals + 1; ++ip)
      s[ip] = 0;
   unsigned long long *t = new unsigned long long[nhumanlabels + 1];
   for (int ih = 0; ih < nhumanlabels + 1; ++ih)
      t[ih] = 0;
   if (print_verbose) { printf("Creating s and t arrays...\n  "); fflush(stdout); }
   for (int iv = 0; iv < nvoxels; ++iv) {
      if (print_verbose) RNProgressBar(iv, nvoxels);
      if (voxel_truths[iv] == 0) continue;
      // just checking
      rn_assertion((0 <= voxel_proposals[iv]) && (voxel_proposals[iv] < nproposals + 1));
      rn_assertion((0 <= voxel_truths[iv]) && (voxel_truths[iv] < nhumanlabels + 1));
      s[voxel_proposals[iv]]++;
      t[voxel_truths[iv]]++;
   }
   if (print_verbose) printf("\ndone!\n");

   // populate overlap matrix
   unsigned long long **c = new unsigned long long *[nproposals + 1];
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      c[ip] = new unsigned long long[nhumanlabels + 1];
      for (int ih = 0; ih < nhumanlabels + 1; ++ih)
         c[ip][ih] = 0;
   }
   if (print_verbose) { printf("Creating overlap matrix...\n  "); fflush(stdout); }
   for (int iv = 0; iv < nvoxels; ++iv) {
      if (print_verbose) RNProgressBar(iv, nvoxels);
      if (voxel_truths[iv] == 0) continue;
      rn_assertion((0 <= voxel_proposals[iv]) && (voxel_proposals[iv] < nproposals + 1));
      rn_assertion((0 <= voxel_truths[iv]) && (voxel_truths[iv] < nhumanlabels + 1));
      c[voxel_proposals[iv]][voxel_truths[iv]]++;
   }
   if (print_verbose) printf("\ndone!\n");

   unsigned long long ncellular_voxels = 0;
   for (int iv = 0; iv < nvoxels; ++iv) {
      if (voxel_truths[iv] != 0) ncellular_voxels++;
   }

   // create normalized arrays
   RNScalar *sp = new RNScalar[nproposals + 1];
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      sp[ip] = (RNScalar)s[ip] / ncellular_voxels;
   }
   RNScalar *tp = new RNScalar[nhumanlabels + 1];
   for (int ih = 0; ih < nhumanlabels + 1; ++ih) {
      tp[ih] = (RNScalar)t[ih] / ncellular_voxels;
   }
   RNScalar **p = new RNScalar *[nproposals + 1];
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      p[ip] = new RNScalar[nhumanlabels + 1];
      for (int ih = 0; ih < nhumanlabels + 1; ++ih) {
         p[ip][ih] = c[ip][ih] / (RNScalar)ncellular_voxels;
      }
   }

   // make sure that certain parameters are maintained
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      unsigned long long scheck = 0;
      RNScalar spcheck = 0;
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         scheck += c[ip][ih];
         spcheck += p[ip][ih];
      }
      rn_assertion(scheck == s[ip]);
      rn_assertion(abs(spcheck == sp[ip]) < RN_MATH_PRECISION);
   }
   for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
      unsigned long long tcheck = 0;
      RNScalar tpcheck = 0;
      for (int ip = 0; ip < nproposals + 1; ++ip) {
         tcheck += c[ip][ih];
         tpcheck += p[ip][ih];
      }
      rn_assertion(tcheck == t[ih]);
      rn_assertion(abs(tpcheck == tp[ih]) < RN_MATH_PRECISION);
   }

   // TP voxel pairs mapped to the same segment within S and T
   // FP voxel pairs mapped to the same segment within S but not T
   // FN voxel pairs mapped to the same segment of T but not S
   // TN voxel pairs mapped to different segments of both S and T
   unsigned long long TP = 0;
   unsigned long long TP_FP = 0;
   unsigned long long TP_FN = 0;
   unsigned long long prev = 0;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         TP += (c[ip][ih] * (c[ip][ih] - 1)) / 2;
         if (TP < prev) { fprintf(stderr, "Error!\n"); exit(-1); }
         prev = TP;
      }
   }
   for (int ip = 0; ip < nproposals + 1; ++ip)
      TP_FP += (s[ip] * (s[ip] - 1)) / 2;
   for (int ih = 1; ih < nhumanlabels + 1; ++ih)
      TP_FN += (t[ih] * (t[ih] - 1)) / 2;

   unsigned long long FP = TP_FP - TP;
   unsigned long long FN = TP_FN - TP;
   unsigned long long Nchoose2 = (ncellular_voxels * (ncellular_voxels - 1)) / 2;
   unsigned long long TN = Nchoose2 - TP - FP - FN;

   // get random index
   RNScalar RI = (TP + TN) / (RNScalar)Nchoose2;
   RNScalar RE = 1 - RI;
   RNScalar REsplit = FN / (RNScalar)Nchoose2;
   RNScalar REmerge = FP / (RNScalar)Nchoose2;
   
   if (print_verbose) {
      printf("Rand Error Full: %lf\n", RE);
      printf("Rand Error Merge: %lf\n", REmerge);
      printf("Rand Error Split: %lf\n", REsplit);
   }
   metrics[0] = RE;
   metrics[1] = REmerge;
   metrics[2] = REsplit;

   RNScalar RFSsplit = 0.0;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         RFSsplit += p[ip][ih] * p[ip][ih];
      }
   }
   RNScalar RFSsplit_denominator = 0.0;
   for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
      RFSsplit_denominator += tp[ih] * tp[ih];
   }
   RFSsplit /= RFSsplit_denominator;

   RNScalar RFSmerge = 0.0;
   RNScalar RFSmerge_denominator = 0.0;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         RFSmerge += p[ip][ih] * p[ip][ih];
      }
      RFSmerge_denominator += sp[ip] * sp[ip];
   }
   RFSmerge /= RFSmerge_denominator;

   RNScalar RFS = 0.0;
   RNScalar RFS_denominator = 0.0;
   RNScalar alpha = 0.5;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         RFS += p[ip][ih] * p[ip][ih];
      }
      RFS_denominator += alpha * sp[ip] * sp[ip];
   }
   for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
      RFS_denominator += (1.0 - alpha) * tp[ih] * tp[ih];
   }
   RFS /= RFS_denominator;
   
   if (print_verbose) {
      printf("Rand F Score Full: %lf\n", RFS);
      printf("Rand F Score Merge: %lf\n", RFSmerge);
      printf("Rand F Score Split: %lf\n", RFSsplit);
   }
   metrics[3] = RFS;
   metrics[4] = RFSmerge;
   metrics[5] = RFSsplit;

   RNScalar HS = 0.0;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      if (sp[ip] == 0.0) continue;
      HS -= sp[ip] * log(sp[ip]);
   }
   RNScalar HT = 0.0;
   for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
      if (tp[ih] == 0.0) continue;
      HT -= tp[ih] * log(tp[ih]);
   }
   RNScalar HST = 0.0;
   RNScalar HTS = 0.0;
   RNScalar IST = 0.0;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         if (p[ip][ih] == 0.0) continue;
         HST -= p[ip][ih] * log(p[ip][ih] / tp[ih]);
         HTS -= p[ip][ih] * log(p[ip][ih] / sp[ip]);
         IST += p[ip][ih] * log(p[ip][ih] / (sp[ip] * tp[ih]));
      }
   }

   RNScalar VIFS = IST / (alpha * HT + (1.0 - alpha) * HS);
   RNScalar VIFSmerge = IST / HT;
   RNScalar VIFSsplit = IST / HS;

   if (print_verbose) {
      printf("Variation F Score Full: %lf\n", VIFS);
      printf("Variation F Score Merge: %lf\n", VIFSmerge);
      printf("Variation F Score Split: %lf\n", VIFSsplit);
   }
   metrics[6] = VIFS;
   metrics[7] = VIFSmerge;
   metrics[8] = VIFSsplit;

   RNScalar VI = 0.0;
   RNScalar VIsplit = HST;
   RNScalar VImerge = HTS;
   for (int ip = 0; ip < nproposals + 1; ++ip) {
      for (int ih = 1; ih < nhumanlabels + 1; ++ih) {
         if (p[ip][ih] == 0.0) continue;
         VI -= p[ip][ih] * (log(p[ip][ih] / sp[ip]) + log(p[ip][ih] / tp[ih]));
      }
   }

   if (print_verbose) {
      printf("Variation of Information Full: %lf\n", VI);
      printf("Variation of Information Merge: %lf\n", VImerge);
      printf("Variation of Information Split: %lf\n", VIsplit);
   }
   metrics[9] = VI;
   metrics[10] = VImerge;
   metrics[11] = VIsplit;

   // free memory
   delete[] voxel_truths;
   delete[] s;
   delete[] t;
   for (int ip = 0; ip < nproposals + 1; ++ip)
      delete[] c[ip];
   delete[] c;
   delete[] sp;
   for (int ip = 0; ip < nproposals + 1; ++ip)
      delete[] p[ip];
   delete[] p;

   // return success
   return 1;
}



int NeuronData::
SegmentationMetric(int *cellular_proposals, RNScalar *metrics) const
{
   // just checking
   rn_assertion(cellular_proposals != NULL);
   rn_assertion(metrics != NULL);

   // get segmentation results
   int true_positives = 0;
   int false_positives = 0;
   int false_negatives = 0;
   int true_negatives = 0;

   // go through all pairs of cellulars
   for (int ic1 = 0; ic1 < NCellulars(); ++ic1) {
      NeuronCellular *cellular_one = Cellular(ic1);
      for (int ic2 = ic1 + 1; ic2 < NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = Cellular(ic2);

         RNBoolean same_proposal = (cellular_proposals[ic1] == cellular_proposals[ic2]);
         RNBoolean same_human_label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());

         // increment the correct counter
         if (same_proposal && same_human_label) true_positives++;
         else if (same_proposal && !same_human_label) false_positives++;
         else if (!same_proposal && same_human_label) false_negatives++;
         else if (!same_proposal && !same_human_label) true_negatives++;
         else rn_assertion(FALSE);
      }
   }

   // set the metrics
   metrics[0] = true_positives;
   metrics[1] = false_positives;
   metrics[2] = false_negatives;
   metrics[3] = true_negatives;
   // merge precision
   metrics[4] = true_positives / (RNScalar)(true_positives + false_positives);
   // merge recall
   metrics[5] = true_positives / (RNScalar)(true_positives + false_negatives);
   // split precision
   metrics[6] = true_negatives / (RNScalar)(true_negatives + false_negatives);
   // split recall
   metrics[7] = true_negatives / (RNScalar)(true_negatives + false_positives);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void NeuronData::
UpdateBBox(void)
{
   for (int is = 0; is < NSupervoxels(); ++is) {
      Supervoxel(is)->ReadVoxels();
      Supervoxel(is)->UpdateBBox();
      Supervoxel(is)->ReleaseVoxels();
   }

   for (int ih = 0; ih < NHumanLabels(); ++ih) {
      HumanLabel(ih)->ReadVoxels();
      HumanLabel(ih)->UpdateBBox();
      HumanLabel(ih)->ReleaseVoxels();
   }

   for (int ip = 0; ip < NPredictions(); ++ip) {
      Prediction(ip)->ReadVoxels();
      Prediction(ip)->UpdateBBox();
      Prediction(ip)->ReleaseVoxels();
   }

   bbox = R3null_box;
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (!Supervoxel(is)->BBox().IsEmpty())
         bbox.Union(Supervoxel(is)->BBox());
   }
}



void NeuronData::
InvalidateBBox(void)
{
   // invalidata bounding box
   if (bbox.XMin() != FLT_MAX) {
      bbox = R3Box(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
   }
}



RNBoolean NeuronData::
DoesBBoxNeedUpdate(void) const
{
   // return whether bounding box needs update
   if (bbox.XMin() == FLT_MAX) return TRUE;
   else return FALSE;
}



void NeuronData::
SetBBox(const R3Box& bbox)
{
   // set bounding box
   this->bbox = bbox;
}



////////////////////////////////////////////////////////////////////////
// Read/Write functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
ReadFile(const char *input_filename, RNBoolean read_voxels, RNBoolean read_affinities)
{
   // parse input filename extension
   const char *extension = strchr(input_filename, '.');
   if (!extension) { fprintf(stderr, "Filename %s has no extension (e.g., .neuron, .txt)\n", input_filename); return 0; }

   if (!strncmp(extension, ".txt", 4)) return ReadASCIIFile(input_filename, read_voxels, read_affinities);
   else if (!strncmp(extension, ".neuron", 7)) return ReadNeuronFile(input_filename, read_voxels, read_affinities);
   else { fprintf(stderr, "Unrecognized file extension %s in %s\n", extension, input_filename); return 0; }
}



int NeuronData::
WriteFile(const char *output_filename)
{
   // parse output filename extension
   const char *extension = strchr(output_filename, '.');
   if (!extension) { fprintf(stderr, "Filename %s has no extension (e.g., .neuron, .txt)\n", output_filename); return 0; }

   if (!strncmp(extension, ".txt", 4)) return WriteASCIIFile(output_filename);
   else if (!strncmp(extension, ".neuron", 7)) return WriteNeuronFile(output_filename);
   else { fprintf(stderr, "Unrecognized file extension %s in %s\n", extension, output_filename); return 0; }
}



////////////////////////////////////////////////////////////////////////
// ASCII to Neuron helper functions
////////////////////////////////////////////////////////////////////////

void NeuronData::
LabelExtracellulars(R3Grid *machine_labels, int ncellulars) const
{
   // go through every voxel in machine labels
   int nextracellular_voxels = 0;
   int *extracellular_mapping = new int[machine_labels->NEntries()];

   for (int ix = 0; ix < machine_labels->XResolution(); ++ix) {
      for (int iy = 0; iy < machine_labels->YResolution(); ++iy) {
         for (int iz = 0; iz < machine_labels->ZResolution(); ++iz) {
            int voxel_index = IndicesToIndex(ix, iy, iz);

            // skip all cellulars that are labeled
            int label = (int)(machine_labels->GridValue(ix, iy, iz) + 0.5);
            if (label > 0) {
               extracellular_mapping[voxel_index] = -1;
            }
            else {
               extracellular_mapping[voxel_index] = nextracellular_voxels;
               // increment extracellular counters
               nextracellular_voxels++;
            }
         }
      }
   }
   // set helper variables
   int nlabeled = 0;
   int nextracellulars = 0;
   int left_to_label = nextracellular_voxels;

   // create a stack of current voxels
   RNQueue<int *> voxel_queue = RNQueue<int *>();

   // continually run depth first search with new lead voxel
   int voxel_index = 0;
   if (NEURON_DEBUG) { printf("Creating extracellular supervoxels from machine labels..."); fflush(stdout); }
   while (left_to_label) {
      // find the lowest voxel yet to belong to extracellular
      while (extracellular_mapping[voxel_index] == -1) voxel_index++;

      voxel_queue.Push(new int(voxel_index));
      extracellular_mapping[voxel_index] = -1;

      // run stack based depth first search from voxel
      while (!voxel_queue.IsEmpty()) {
         int *popped_index = voxel_queue.Pop();
         int index = *popped_index;
         delete popped_index;

         rn_assertion(machine_labels->GridValue(index) < 1);
         // add the one so extracellulars not added to existing cellular
         machine_labels->SetGridValue(index, ncellulars + nextracellulars + 1);

         left_to_label--;
         nlabeled++;

         int ix, iy, iz;
         IndexToIndices(index, ix, iy, iz);
         // go through x coordinates
         for (int dx = -1; dx <= 1; ++dx) {
            // check x dimensions
            int neighbor_x = ix + dx;
            if (neighbor_x < 0) continue;
            if (neighbor_x > resolution[RN_X] - 1) continue;

            // go through y coordinates
            for (int dy = -1; dy <= 1; ++dy) {
               // check y dimensions
               int neighbor_y = iy + dy;
               if (neighbor_y < 0) continue;
               if (neighbor_y > resolution[RN_Y] - 1) continue;

               // go through z coordinates
               for (int dz = -1; dz <= 1; ++dz) {
                  // check z dimensions
                  int neighbor_z = iz + dz;
                  if (neighbor_z < 0) continue;
                  if (neighbor_z > resolution[RN_Z] - 1) continue;

                  int neighbor_index = IndicesToIndex(neighbor_x, neighbor_y, neighbor_z);
                  if (extracellular_mapping[neighbor_index] == -1) continue;

                  // reset extracellular mappiing so this neighbor is not revisited
                  extracellular_mapping[neighbor_index] = -1;

                  // push the voxel onto the queue
                  voxel_queue.Push(new int(neighbor_index));
               }
            }
         }
      }
      // increment the number of extracellulars
      nextracellulars++;
   }
   if (NEURON_DEBUG) printf("done.\n");

   delete[] extracellular_mapping;
}



void NeuronData::
CreateSupervoxels(R3Grid *machine_labels)
{
   int ncellulars = (int)(machine_labels->Maximum() + 0.5);

   // update the extracellular in machine labels to be non zero
   LabelExtracellulars(machine_labels, ncellulars);

   int nsupervoxels = (int)(machine_labels->Maximum() + 0.5);
   rn_assertion((int)(machine_labels->Minimum() + 0.5) > 0);

   // create array of number of voxels for each supervoxel
   int *supervoxel_nvoxels = new int[nsupervoxels];
   int *nvoxels_seen = new int[nsupervoxels];
   for (int is = 0; is < nsupervoxels; ++is) {
      supervoxel_nvoxels[is] = 0;
      nvoxels_seen[is] = 0;
   }

   // go through every index of the machine labels grid
   for (int ix = 0; ix < machine_labels->XResolution(); ++ix) {
      for (int iy = 0; iy < machine_labels->YResolution(); ++iy) {
         for (int iz = 0; iz < machine_labels->ZResolution(); ++iz) {
            int label = (int)(machine_labels->GridValue(ix, iy, iz) + 0.5) - 1;
            rn_assertion((0 <= label) && (label < nsupervoxels));
            supervoxel_nvoxels[label]++;
         }
      }
   }

   // create voxels for every supervoxel
   for (int is = 0; is < nsupervoxels; ++is) {
      if (is < ncellulars) {
         NeuronCellular *cellular = new NeuronCellular();
         InsertCellular(cellular);
         cellular->CreateVoxels(supervoxel_nvoxels[is]);
      }
      else {
         NeuronExtracellular *extracellular = new NeuronExtracellular();
         InsertExtracellular(extracellular);
         extracellular->CreateVoxels(supervoxel_nvoxels[is]);
      }
   }

   // go through grids to initialize voxel mapping
   if (NEURON_DEBUG) printf("Setting voxel mapping for supervoxels...\n  ");
   for (int iv = 0; iv < machine_labels->NEntries(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, machine_labels->NEntries());
      int supervoxel_index = (int)(machine_labels->GridValue(iv) + 0.5) - 1;
      voxel_mapping[iv] = supervoxel_index;

      // set voxel mapping within supervoxel
      NeuronSupervoxel *supervoxel = Supervoxel(supervoxel_index);
      supervoxel->SetVoxelMapping(nvoxels_seen[supervoxel_index], iv);
      nvoxels_seen[supervoxel_index]++;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // go through grid to set preliminary values
   if (NEURON_DEBUG) printf("Setting preliminary voxel values...\n  ");
   for (int ix = 0; ix < machine_labels->XResolution(); ++ix) {
      RNProgressBar(ix, machine_labels->XResolution());
      for (int iy = 0; iy < machine_labels->YResolution(); ++iy) {
         for (int iz = 0; iz < machine_labels->ZResolution(); ++iz) {
            int voxel_index = IndicesToIndex(ix, iy, iz);
            // get voxel and set data
            NeuronVoxel *voxel = Voxel(ix, iy, iz);
            voxel->data = this;

            // get supervoxel
            NeuronSupervoxel *supervoxel = Supervoxel(voxel_mapping[voxel_index]);
            voxel->supervoxel = supervoxel;
            voxel->supervoxel_index = supervoxel->LocalIndex(voxel_index);
            rn_assertion(voxel->supervoxel_index != -1);

            // set voxel coordinates
            voxel->coordinates[RN_X] = ix;
            voxel->coordinates[RN_Y] = iy;
            voxel->coordinates[RN_Z] = iz;

            // set boundary conditions for supervoxel
            if (voxel->IsOnBoundary()) supervoxel->SetBoundaryFlag(TRUE);
         }
      }
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // update bboxes for all supervoxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      Supervoxel(is)->UpdateBBox();
   }

   // free memory
   delete[] supervoxel_nvoxels;
   delete[] nvoxels_seen;
}



void NeuronData::
SetAffinities(R3Grid *affinities[3])
{
   // go through each grid to update affinity
   for (int dim = 0; dim <= 2; ++dim) {
      if (NEURON_DEBUG) printf("Updating voxel affinity values for dimension %d...\n  ", dim);
      for (int ix = 0; ix < affinities[dim]->XResolution(); ++ix) {
         if (NEURON_DEBUG) RNProgressBar(ix, affinities[dim]->XResolution());
         for (int iy = 0; iy < affinities[dim]->YResolution(); ++iy) {
            for (int iz = 0; iz < affinities[dim]->ZResolution(); ++iz) {
               NeuronVoxel *voxel = Voxel(ix, iy, iz);
               voxel->affinities[dim] = affinities[dim]->GridValue(ix, iy, iz);
            }
         }
      }
      if (NEURON_DEBUG) printf("\ndone.\n");
   }
}



void NeuronData::
SetImageValues(R3Grid *image)
{
   if (NEURON_DEBUG) printf("Updating image values...\n  ");
   for (int ix = 0; ix < image->XResolution(); ++ix) {
      if (NEURON_DEBUG) RNProgressBar(ix, image->XResolution());
      for (int iy = 0; iy < image->YResolution(); ++iy) {
         for (int iz = 0; iz < image->ZResolution(); ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);
            voxel->image = image->GridValue(ix, iy, iz);
         }
      }
   }
   if (NEURON_DEBUG) printf("\ndone.\n");
}



void NeuronData::
CreateHumanLabels(R3Grid *human_labels)
{
   int nhuman_labels = (int)(human_labels->Maximum() + 0.5);

   // create all of the human labels
   for (int ih = 0; ih < nhuman_labels; ++ih) {
      NeuronHumanLabel *human_label = new NeuronHumanLabel();
      InsertHumanLabel(human_label);
   }

   int *human_labels_nvoxels = new int[nhuman_labels];
   for (int ih = 0; ih < nhuman_labels; ++ih)
      human_labels_nvoxels[ih] = 0;

   // go through every grid voxel
   if (NEURON_DEBUG) printf("Creating human labels...\n  ");
   for (int iv = 0; iv < human_labels->NEntries(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, human_labels->NEntries());
      int label = (int)(human_labels->GridValue(iv) + 0.5) - 1;
      if (label < 0) continue;
      human_labels_nvoxels[label]++;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // create voxel mappings
   for (int ih = 0; ih < nhuman_labels; ++ih) {
      NeuronHumanLabel *human_label = HumanLabel(ih);
      human_label->CreateVoxelMapping(human_labels_nvoxels[ih]);
   }

   // set voxel human label mapping
   int *nvoxels_seen = new int[nhuman_labels];
   for (int ih = 0; ih < nhuman_labels; ++ih)
      nvoxels_seen[ih] = 0;

   // set voxel mapping
   if (NEURON_DEBUG) printf("Setting voxel mapping for human labels...\n  ");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      int human_label_index = (int)(human_labels->GridValue(iv) + 0.5) - 1;
      if (human_label_index < 0) continue;

      // set voxel mapping within human labels
      NeuronHumanLabel *human_label = HumanLabel(human_label_index);
      human_label->SetVoxelMapping(nvoxels_seen[human_label_index], iv);
      nvoxels_seen[human_label_index]++;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // go through grid to set preliminary values
   if (NEURON_DEBUG) printf("Setting human label voxel values...\n  ");
   for (int ix = 0; ix < human_labels->XResolution(); ++ix) {
      if (NEURON_DEBUG) RNProgressBar(ix, human_labels->XResolution());
      for (int iy = 0; iy < human_labels->YResolution(); ++iy) {
         for (int iz = 0; iz < human_labels->ZResolution(); ++iz) {
            int voxel_index = IndicesToIndex(ix, iy, iz);
            int human_label_index = (int)(human_labels->GridValue(ix, iy, iz) + 0.5) - 1;
            if (human_label_index < 0) continue;

            // get voxel and set data
            NeuronVoxel *voxel = Voxel(ix, iy, iz);

            // get human_label
            NeuronHumanLabel *human_label = HumanLabel(human_label_index);
            voxel->human_label = human_label;
            voxel->human_label_index = human_label->LocalIndex(voxel_index);
            rn_assertion(voxel->human_label_index != -1);

            // set voxel coordinates
            voxel->coordinates[RN_X] = ix;
            voxel->coordinates[RN_Y] = iy;
            voxel->coordinates[RN_Z] = iz;

            // set boundary conditions for human label
            if (voxel->IsOnBoundary()) human_label->SetBoundaryFlag(TRUE);
         }
      }
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) printf("Adding supervoxels to reconstruction objects and vice versa...\n  ");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      NeuronVoxel *voxel = Voxel(iv);
      NeuronSupervoxel *supervoxel = voxel->supervoxel;
      NeuronHumanLabel *human_label = voxel->human_label;
      if (!human_label) continue;

      human_label->InsertSupervoxel(supervoxel);
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   delete[] human_labels_nvoxels;
   delete[] nvoxels_seen;
}



// create arrays of coordinates
struct Coordinate {
   // constructors
   Coordinate(int voxel_index, int coordinate) {
      this->voxel_index = voxel_index;
      this->coordinate = coordinate;
   }

   // instance variables
   int voxel_index;
   int coordinate;
};



int CoordinateSort(Coordinate a, Coordinate b) {
   return (a.coordinate < b.coordinate);
}



void NeuronData::
CreatePredictions(R3Grid *predictions)
{
   int npredictions = (int)(predictions->Maximum() + 0.5);

   // create all of the predictions
   for (int ip = 0; ip < npredictions; ++ip) {
      NeuronPrediction *prediction = new NeuronPrediction();
      InsertPrediction(prediction);
      prediction->read_voxel_count++;
   }

   int *prediction_nvoxels = new int[npredictions];
   for (int ip = 0; ip < npredictions; ++ip)
      prediction_nvoxels[ip] = 0;

   // go through every grid voxel
   if (NEURON_DEBUG) printf("Creating predictions...\n  ");
   for (int iv = 0; iv < predictions->NEntries(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, predictions->NEntries());
      int label = (int)(predictions->GridValue(iv) + 0.5) - 1;
      if (label < 0) continue;
      prediction_nvoxels[label]++;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // create voxel mappings
   for (int ip = 0; ip < npredictions; ++ip) {
      NeuronPrediction *prediction = Prediction(ip);
      prediction->CreateVoxelMapping(prediction_nvoxels[ip]);
   }

   // set voxel prediction mapping
   int *nvoxels_seen = new int[npredictions];
   for (int ip = 0; ip < npredictions; ++ip)
      nvoxels_seen[ip] = 0;

   // set voxel mapping
   if (NEURON_DEBUG) printf("Setting voxel mapping for predictions...\n  ");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      int prediction_index = (int)(predictions->GridValue(iv) + 0.5) - 1;
      if (prediction_index < 0) continue;

      // set voxel mapping within predictions
      NeuronPrediction *prediction = Prediction(prediction_index);
      prediction->SetVoxelMapping(nvoxels_seen[prediction_index], iv);
      nvoxels_seen[prediction_index]++;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // go through grid to set preliminary values
   if (NEURON_DEBUG) printf("Setting prediction voxel values...\n  ");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      int prediction_index = (int)(predictions->GridValue(iv) + 0.5) - 1;
      if (prediction_index < 0) continue;

      // get voxel and set data
      NeuronVoxel *voxel = Voxel(iv);

      // get prediction
      NeuronPrediction *prediction = Prediction(prediction_index);
      voxel->prediction = prediction;
      voxel->prediction_index = prediction->LocalIndex(iv);
      rn_assertion(voxel->prediction_index != -1);

      // set boundary conditions for prediction
      if (voxel->IsOnBoundary()) prediction->SetBoundaryFlag(TRUE);
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) printf("Adding supervoxels to prediction objects and vice versa...\n  ");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      NeuronVoxel *voxel = Voxel(iv);
      NeuronSupervoxel *supervoxel = voxel->supervoxel;
      NeuronPrediction *prediction = voxel->prediction;
      if (!prediction) continue;

      prediction->InsertSupervoxel(supervoxel);
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) { printf("Mapping predictions to voxels...\n  "); fflush(stdout); }
   // for every supervoxel, find the best voxel
   for (int ip = 0; ip < NPredictions(); ++ip) {
      if (NEURON_DEBUG) RNProgressBar(ip, NPredictions());
      NeuronPrediction *prediction = Prediction(ip);

      // make sure the voxels are in memory
      rn_assertion(prediction->AreVoxelsResident());

      // get number of voxels
      int nvoxels = prediction->NVoxels();

      std::vector<Coordinate> xcoordinates = std::vector<Coordinate>();
      std::vector<Coordinate> ycoordinates = std::vector<Coordinate>();
      std::vector<Coordinate> zcoordinates = std::vector<Coordinate>();

      unsigned long long *total_cost = new unsigned long long[nvoxels];

      for (int iv = 0; iv < nvoxels; ++iv) {
         NeuronVoxel *voxel = prediction->LocalVoxel(iv);
         Coordinate xcoord = Coordinate(iv, voxel->XCoordinate());
         Coordinate ycoord = Coordinate(iv, voxel->YCoordinate());
         Coordinate zcoord = Coordinate(iv, voxel->ZCoordinate());

         // add to the coordinate vectors
         xcoordinates.push_back(xcoord);
         ycoordinates.push_back(ycoord);
         zcoordinates.push_back(zcoord);
         total_cost[iv] = 0;
      }

      // sort the coordinates on the axes
      sort(xcoordinates.begin(), xcoordinates.end(), &CoordinateSort);
      sort(ycoordinates.begin(), ycoordinates.end(), &CoordinateSort);
      sort(zcoordinates.begin(), zcoordinates.end(), &CoordinateSort);

      // calculate the distance between this point and every point to the left and right
      for (int dim = 0; dim <= 2; ++dim) {
         std::vector<Coordinate> *coordinates;
         if (dim == RN_X) coordinates = &xcoordinates;
         else if (dim == RN_Y) coordinates = &ycoordinates;
         else coordinates = &zcoordinates;

         unsigned long long *left_cost = new unsigned long long[nvoxels];
         unsigned long long *right_cost = new unsigned long long[nvoxels];

         left_cost[0] = 0;
         for (int iv = 1; iv < nvoxels; ++iv) {
            // get the difference between the last coordinate and this
            int diff = coordinates->at(iv).coordinate - coordinates->at(iv - 1).coordinate;
            rn_assertion(diff >= 0);

            left_cost[iv] = left_cost[iv - 1] + diff * iv;
         }

         right_cost[nvoxels - 1] = 0;
         for (int iv = nvoxels - 2; iv >= 0; --iv) {
            // get the difference between the last coordinate and this
            int diff = coordinates->at(iv + 1).coordinate - coordinates->at(iv).coordinate;
            rn_assertion(diff >= 0);

            right_cost[iv] = right_cost[iv + 1] + diff * ((nvoxels - 1) - iv);
         }

         // update the total cost
         for (int iv = 0; iv < nvoxels; ++iv) {
            // get index for this coordinate
            int voxel_index = coordinates->at(iv).voxel_index;

            // make sure the costs are positive (sanity check)
            rn_assertion(left_cost[iv] >= 0);
            rn_assertion(right_cost[iv] >= 0);

            total_cost[voxel_index] += left_cost[iv] + right_cost[iv];
         }

         // free memory
         delete[] left_cost;
         delete[] right_cost;
      }

      // get the best voxel from total cost array
      unsigned long long manhattan_minimum = LLONG_MAX;
      int manhattan_index = -1;
      for (int iv = 0; iv < nvoxels; ++iv) {
         if (total_cost[iv] < manhattan_minimum) {
            manhattan_minimum = total_cost[iv];
            manhattan_index = iv;
         }
      }

      // update center voxel and center human label
      prediction->center_voxel_index = manhattan_index;
      if (prediction->CenterVoxel()->HumanLabel()) prediction->center_human_label_index = prediction->CenterVoxel()->HumanLabel()->DataIndex();
      else prediction->center_human_label_index = -1;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // free memory
   delete[] prediction_nvoxels;
   delete[] nvoxels_seen;

   // go through every supervoxel
   if (NEURON_DEBUG) { printf("Creating a mapping from predictions to human labels...\n  "); fflush(stdout); }
   for (int ip = 0; ip < NPredictions(); ++ip) {
      if (NEURON_DEBUG) { RNProgressBar(ip, NPredictions()); }
      NeuronPrediction *prediction = Prediction(ip);

      // create counter for all human label occurrences within cellular
      int *human_label_count = new int[NHumanLabels()];
      for (int ih = 0; ih < NHumanLabels(); ++ih) {
         human_label_count[ih] = 0;
      }

      // go through all supervoxel voxels
      for (int iv = 0; iv < prediction->NVoxels(); ++iv) {
         NeuronVoxel *voxel = prediction->LocalVoxel(iv);
         if (voxel->HumanLabel()) human_label_count[voxel->HumanLabel()->DataIndex()]++;
      }

      // find the most prominent human label
      int best_human_label_value = 0;
      int best_human_label_index = -1;

      for (int ih = 0; ih < NHumanLabels(); ++ih) {
         if (human_label_count[ih] > best_human_label_value) {
            best_human_label_value = human_label_count[ih];
            best_human_label_index = ih;
         }
      }

      if (best_human_label_index != -1) prediction->majority_human_label_index = best_human_label_index;
      else prediction->majority_human_label_index = -1;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   // create prediction boundaries
   CreatePredictionBoundaries();
}



void NeuronData::
CreatePredictionBoundaries(void)
{
   std::vector<std::vector<RNScalar> > prediction_boundary_affinities = std::vector<std::vector<RNScalar> >();

   if (NEURON_DEBUG) { printf("Determining boundaries between predictions...\n  "); fflush(stdout); }
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      NeuronVoxel *voxel = Voxel(iv);
      NeuronPrediction *prediction = voxel->Prediction();
      if (!prediction) continue;

      for (int in = 0; in < voxel->NNeighbors() / 2; ++in) {
         NeuronVoxel *neighbor = voxel->NextVoxel(in);
         if (!neighbor) continue;

         NeuronPrediction *neighbor_prediction = neighbor->Prediction();
         if (!neighbor_prediction) continue;
         if (prediction == neighbor_prediction) continue;

         // create boundary if one does not exist
         NeuronPredictionBoundary *prediction_boundary = NULL;
         for (int ib = 0; ib < prediction->NBoundaries(); ++ib) {
            NeuronPredictionBoundary *candidate = prediction->Boundary(ib);
            NeuronPrediction *other = candidate->OtherPrediction(prediction);
            
            // check for equality
            if (other == neighbor_prediction) { prediction_boundary = candidate; break; }
         }

         // create new prediction boundary if it is null
         if (!prediction_boundary) {
            prediction_boundary = new NeuronPredictionBoundary();
            prediction_boundary_affinities.push_back(std::vector<RNScalar>());

            prediction->InsertBoundary(prediction_boundary);
            neighbor_prediction->InsertBoundary(prediction_boundary);

            InsertPredictionBoundary(prediction_boundary);
         }

         // add the voxel affinity
         RNScalar affinity = voxel->AffinityToNeighbor(neighbor);
         prediction_boundary_affinities[prediction_boundary->DataIndex()].push_back(affinity);
      }
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) { printf("Creating all prediction boundary affinities...\n  "); fflush(stdout); }
   for (int ipb = 0; ipb < NPredictionBoundaries(); ++ipb) {
      if (NEURON_DEBUG) RNProgressBar(ipb, NPredictionBoundaries());
      PredictionBoundary(ipb)->CreateNAffinities(prediction_boundary_affinities[ipb].size());
      for (unsigned int ia = 0; ia < prediction_boundary_affinities[ipb].size(); ++ia) {
         PredictionBoundary(ipb)->UpdateVoxelAffinity(prediction_boundary_affinities[ipb][ia], ia);
      }
      prediction_boundary_affinities[ipb].clear();
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) { printf("Updating prediction boundary statistics...\n  "); fflush(stdout); }
   for (int ipb = 0; ipb < NPredictionBoundaries(); ++ipb) {
      if (NEURON_DEBUG) RNProgressBar(ipb, NPredictionBoundaries());
      PredictionBoundary(ipb)->UpdateStatistics();
   }
   if (NEURON_DEBUG) printf("\ndone.\n");
}



void NeuronData::
CreateBoundaries(void)
{
   std::vector<std::vector<RNScalar> > boundary_affinities = std::vector<std::vector<RNScalar> >();

   // go through every supervoxel
   if (NEURON_DEBUG) { printf("Determining boundaries between supervoxels...\n  "); fflush(stdout); }
   for (int iv = 0; iv < NVoxels(); ++iv) {
      if (NEURON_DEBUG) RNProgressBar(iv, NVoxels());
      NeuronVoxel *voxel = Voxel(iv);
      NeuronSupervoxel *supervoxel = voxel->Supervoxel();

      for (int in = 0; in < voxel->NNeighbors() / 2; ++in) {
         NeuronVoxel *neighbor = voxel->NextVoxel(in);
         if (!neighbor) continue;

         NeuronSupervoxel *neighbor_supervoxel = neighbor->Supervoxel();
         if (supervoxel == neighbor_supervoxel) continue;

         // create boundary if one does not exist
         NeuronBoundary *boundary = NULL;
         for (int ib = 0; ib < supervoxel->NBoundaries(); ++ib) {
            NeuronBoundary *candidate = supervoxel->Boundary(ib);
            NeuronSupervoxel *other = candidate->OtherSupervoxel(supervoxel);

            // check for equality
            if (other == neighbor_supervoxel) { boundary = candidate; break; }
         }

         // create new boundary if it is null
         if (!boundary) {
            boundary = new NeuronBoundary();
            boundary_affinities.push_back(std::vector<RNScalar>());

            supervoxel->InsertBoundary(boundary);
            neighbor_supervoxel->InsertBoundary(boundary);

            InsertBoundary(boundary);
         }

         // add the voxel affinity
         RNScalar affinity = voxel->AffinityToNeighbor(neighbor);
         boundary_affinities[boundary->DataIndex()].push_back(affinity);
      }
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) { printf("Creating all boundary affinities...\n  "); fflush(stdout); }
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      if (NEURON_DEBUG) RNProgressBar(ib, NBoundaries());
      Boundary(ib)->CreateNAffinities(boundary_affinities[ib].size());
      for (unsigned int ia = 0; ia < boundary_affinities[ib].size(); ++ia) {
         Boundary(ib)->UpdateVoxelAffinity(boundary_affinities[ib][ia], ia);
      }
      boundary_affinities[ib].clear();
   }
   if (NEURON_DEBUG) printf("\ndone.\n");

   if (NEURON_DEBUG) { printf("Updating boundary statistics...\n  "); fflush(stdout); }
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      if (NEURON_DEBUG) RNProgressBar(ib, NBoundaries());
      Boundary(ib)->UpdateStatistics();
   }
   if (NEURON_DEBUG) printf("\ndone.\n");
}



void NeuronData::
CreateHumanLabelMapping(void)
{
   // go through every supervoxel
   if (NEURON_DEBUG) { printf("Creating a mapping from supervoxels to human labels...\n  "); fflush(stdout); }
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (NEURON_DEBUG) { RNProgressBar(is, NSupervoxels()); }
      NeuronSupervoxel *supervoxel = Supervoxel(is);

      // create counter for all human label occurrences within cellular
      int *human_label_count = new int[NHumanLabels()];
      for (int ih = 0; ih < NHumanLabels(); ++ih) {
         human_label_count[ih] = 0;
      }

      // go through all supervoxel voxels
      for (int iv = 0; iv < supervoxel->NVoxels(); ++iv) {
         NeuronVoxel *voxel = supervoxel->LocalVoxel(iv);
         if (voxel->HumanLabel()) human_label_count[voxel->HumanLabel()->DataIndex()]++;
      }

      // find the most prominent human label
      int best_human_label_value = 0;
      int best_human_label_index = -1;

      for (int ih = 0; ih < NHumanLabels(); ++ih) {
         if (human_label_count[ih] > best_human_label_value) {
            best_human_label_value = human_label_count[ih];
            best_human_label_index = ih;
         }
      }

      if (best_human_label_index != -1) supervoxel->majority_human_label_index = best_human_label_index;
      else supervoxel->majority_human_label_index = -1;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");
}


void NeuronData::
CreateSupervoxelMapping(void)
{
   if (NEURON_DEBUG) { printf("Mapping supervoxels to voxels...\n  "); fflush(stdout); }
   // for every supervoxel, find the best voxel
   for (int is = 0; is < NSupervoxels(); ++is) {
      if (NEURON_DEBUG) RNProgressBar(is, NSupervoxels());
      NeuronSupervoxel *supervoxel = Supervoxel(is);

      // make sure the voxels are in memory
      rn_assertion(supervoxel->AreVoxelsResident());

      // get number of voxels
      int nvoxels = supervoxel->NVoxels();



      std::vector<Coordinate> xcoordinates = std::vector<Coordinate>();
      std::vector<Coordinate> ycoordinates = std::vector<Coordinate>();
      std::vector<Coordinate> zcoordinates = std::vector<Coordinate>();

      unsigned long long *total_cost = new unsigned long long[nvoxels];

      for (int iv = 0; iv < nvoxels; ++iv) {
         NeuronVoxel *voxel = supervoxel->LocalVoxel(iv);
         Coordinate xcoord = Coordinate(iv, voxel->XCoordinate());
         Coordinate ycoord = Coordinate(iv, voxel->YCoordinate());
         Coordinate zcoord = Coordinate(iv, voxel->ZCoordinate());

         // add to the coordinate vectors
         xcoordinates.push_back(xcoord);
         ycoordinates.push_back(ycoord);
         zcoordinates.push_back(zcoord);
         total_cost[iv] = 0;
      }

      // sort the coordinates on the axes
      sort(xcoordinates.begin(), xcoordinates.end(), &CoordinateSort);
      sort(ycoordinates.begin(), ycoordinates.end(), &CoordinateSort);
      sort(zcoordinates.begin(), zcoordinates.end(), &CoordinateSort);

      // calculate the distance between this point and every point to the left and right
      for (int dim = 0; dim <= 2; ++dim) {
         std::vector<Coordinate> *coordinates;
         if (dim == RN_X) coordinates = &xcoordinates;
         else if (dim == RN_Y) coordinates = &ycoordinates;
         else coordinates = &zcoordinates;

         unsigned long long *left_cost = new unsigned long long[nvoxels];
         unsigned long long *right_cost = new unsigned long long[nvoxels];

         left_cost[0] = 0;
         for (int iv = 1; iv < nvoxels; ++iv) {
            // get the difference between the last coordinate and this
            int diff = coordinates->at(iv).coordinate - coordinates->at(iv - 1).coordinate;
            rn_assertion(diff >= 0);

            left_cost[iv] = left_cost[iv - 1] + diff * iv;
         }

         right_cost[nvoxels - 1] = 0;
         for (int iv = nvoxels - 2; iv >= 0; --iv) {
            // get the difference between the last coordinate and this
            int diff = coordinates->at(iv + 1).coordinate - coordinates->at(iv).coordinate;
            rn_assertion(diff >= 0);

            right_cost[iv] = right_cost[iv + 1] + diff * ((nvoxels - 1) - iv);
         }

         // update the total cost
         for (int iv = 0; iv < nvoxels; ++iv) {
            // get index for this coordinate
            int voxel_index = coordinates->at(iv).voxel_index;

            // make sure the costs are positive (sanity check)
            rn_assertion(left_cost[iv] >= 0);
            rn_assertion(right_cost[iv] >= 0);

            total_cost[voxel_index] += left_cost[iv] + right_cost[iv];
         }

         // free memory
         delete[] left_cost;
         delete[] right_cost;
      }

      // get the best voxel from total cost array
      unsigned long long manhattan_minimum = LLONG_MAX;
      int manhattan_index = -1;
      for (int iv = 0; iv < nvoxels; ++iv) {
         if (total_cost[iv] < manhattan_minimum) {
            manhattan_minimum = total_cost[iv];
            manhattan_index = iv;
         }
      }

      // update center voxel and center human label
      supervoxel->center_voxel_index = manhattan_index;
      if (supervoxel->CenterVoxel()->HumanLabel()) supervoxel->center_human_label_index = supervoxel->CenterVoxel()->HumanLabel()->DataIndex();
      else supervoxel->center_human_label_index = -1;
   }
   if (NEURON_DEBUG) printf("\ndone.\n");
}



////////////////////////////////////////////////////////////////////////
// ASCII Read/Write functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
ReadASCIIFile(const char *ascii_filename, RNBoolean read_voxels, RNBoolean read_affinities)
{
   // open file
   FILE *fp = fopen(ascii_filename, "r");
   if (!fp) { fprintf(stderr, "Failed to open file: %s\n", ascii_filename); return 0; }

   // get file path and filename
   char path[4096];
   strncpy(path, ascii_filename, 4096);

   char *file = strchr(path, '/');
   *file = '\0'; file++;

   strncpy(file_path, path, 4096);
   strncpy(filename, file, 4096);

   // get affinities filenames
   fgets(affinities_filename, 4096, fp);

   char *affinities_extp = strchr(affinities_filename, '.');
   *affinities_extp = '\0';

   // get human labels filenames
   fgets(human_labels_filename, 4096, fp);

   // remove new line character
   human_labels_filename[strlen(human_labels_filename) - 1] = '\0';
   if (human_labels_filename[strlen(human_labels_filename) - 1] == '\r')
      human_labels_filename[strlen(human_labels_filename) - 1] = '\0';

   char *human_labels_extp = strchr(human_labels_filename, '.');

   // get image filenames
   fgets(image_filename, 4096, fp);

   char *image_extp = strchr(image_filename, '.');
   *image_extp = '\0';

   // get machine labels filename
   fgets(machine_labels_filename, 4096, fp);

   char *machine_labels_extp = strchr(machine_labels_filename, '.');
   *machine_labels_extp = '\0';

   // get scaling factor
   fscanf(fp, "Scale(%lf,%lf,%lf)", &(scaling[RN_X]), &(scaling[RN_Y]), &(scaling[RN_Z]));

   // close file
   fclose(fp);

   char affinities_root_filename[4096];
   sprintf(affinities_root_filename, "%s/%s", file_path, affinities_filename);
   char image_root_filename[4096];
   sprintf(image_root_filename, "%s/%s", file_path, image_filename);
   char machine_labels_root_filename[4096];
   sprintf(machine_labels_root_filename, "%s/%s", file_path, machine_labels_filename);

   if (NEURON_DEBUG) { printf("Reading machine labels file..."); fflush(stdout); }
   RNTime machine_labels_time;
   machine_labels_time.Read();
   R3Grid *machine_labels = RNReadNeuronMetaRawFile(machine_labels_root_filename);
   if (!machine_labels) return 0;
   if (NEURON_DEBUG) printf("done in %0.2f seconds.\n", machine_labels_time.Elapsed());

   // update resolution instance variables
   resolution[RN_X] = machine_labels->XResolution();
   resolution[RN_Y] = machine_labels->YResolution();
   resolution[RN_Z] = machine_labels->ZResolution();

   // update various size variables
   data_size = resolution[RN_X] * resolution[RN_Y] * resolution[RN_Z];
   data_sheet_size = resolution[RN_X] * resolution[RN_Y];
   data_row_size = resolution[RN_X];

   // create voxel_mapping
   voxel_mapping = new int[data_size];
   for (int iv = 0; iv < data_size; ++iv) {
      voxel_mapping[iv] = -1;
   }

   // set up all supervoxels
   CreateSupervoxels(machine_labels);

   // free memory
   delete machine_labels;

   // read in affinities, human_labels, image, and machine_labels files
   if (NEURON_DEBUG) { printf("Reading affinities files..."); fflush(stdout); }
   RNTime affinity_time;
   affinity_time.Read();
   R3Grid **affinities = RNReadNeuronMetaRawFile(affinities_root_filename, TRUE);
   if (!affinities) return 0;
   if (NEURON_DEBUG) printf("done in %0.2f seconds.\n", affinity_time.Elapsed());

   for (int dim = 0; dim <= 2; ++dim) {
      rn_assertion(affinities[dim]->XResolution() == resolution[RN_X]);
      rn_assertion(affinities[dim]->YResolution() == resolution[RN_Y]);
      rn_assertion(affinities[dim]->ZResolution() == resolution[RN_Z]);
   }

   // set affinities
   SetAffinities(affinities);

   // free memory
   for (int dim = 0; dim <= 2; ++dim)
      delete affinities[dim];
   delete[] affinities;

   if (NEURON_DEBUG) { printf("Reading image file..."); fflush(stdout); }
   RNTime image_time;
   image_time.Read();
   R3Grid *image = RNReadNeuronMetaRawFile(image_root_filename);
   if (!image) return 0;
   if (NEURON_DEBUG) printf("done in %0.2f seconds.\n", image_time.Elapsed());

   rn_assertion(image->XResolution() == resolution[RN_X]);
   rn_assertion(image->YResolution() == resolution[RN_Y]);
   rn_assertion(image->ZResolution() == resolution[RN_Z]);

   // set image values
   SetImageValues(image);

   // free memory
   delete image;

   // make sure voxel mapping set
   for (int iv = 0; iv < data_size; ++iv)
      rn_assertion(voxel_mapping[iv] != -1);

   // see if truth is in .txt format or .meta format
   RNBoolean complete_truth = FALSE;

   if (!strcmp(human_labels_extp, ".txt")) complete_truth = FALSE;
   else if (!strcmp(human_labels_extp, ".meta")) complete_truth = TRUE;
   else { fprintf(stderr, "Failed to recognize file extension %s.\n", human_labels_extp); return 0; }

   *human_labels_extp = '\0';
   char human_labels_root_filename[4096];
   sprintf(human_labels_root_filename, "%s/%s", file_path, human_labels_filename);

   if (NEURON_DEBUG) { printf("Reading human labels file..."); fflush(stdout); }
   R3Grid *human_labels = NULL;
   if (complete_truth) {
      human_labels = RNReadNeuronMetaRawFile(human_labels_root_filename);
   }
   else {
      // create the root filename before removing the extension
      char human_labels_root_filename[4096];
      sprintf(human_labels_root_filename, "%s/%s.txt", file_path, human_labels_filename);

      human_labels = ReadTruthTxtFile(human_labels_root_filename);

      // write RNMetaRaw file
      human_labels_extp = strchr(human_labels_root_filename, '.');
      *human_labels_extp = '\0';

      RNMeta meta = RNMeta("Uint16", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 1);
      RNWriteNeuronMetaRawFile(human_labels_root_filename, meta, human_labels);
   }
   if (!human_labels) return 0;
   if (NEURON_DEBUG) printf("done.\n");

   // set human labels
   CreateHumanLabels(human_labels);

   rn_assertion(human_labels->XResolution() == resolution[RN_X]);
   rn_assertion(human_labels->YResolution() == resolution[RN_Y]);
   rn_assertion(human_labels->ZResolution() == resolution[RN_Z]);

   // mapping cellulars to human labels
   CreateHumanLabelMapping();

   // mapping supervoxels to voxels
   CreateSupervoxelMapping();

   // free all temporary memory
   delete human_labels;

   // create supervoxel boundaries
   CreateBoundaries();

   // return OK status
   return 1;
}



int NeuronData::
WriteASCIIFile(const char *ascii_filename)
{
   // open file
   FILE *fp = fopen(ascii_filename, "w");
   if (!fp) { fprintf(stderr, "Failed to open file: %s\n", ascii_filename); return 0; }

   // get file path and filename
   char path[4096];
   strncpy(path, ascii_filename, 4096);

   char *file = strchr(path, '/');
   *file = '\0'; file++;

   // remove extension
   char *extp = strchr(file, '.');
   *extp = '\0';

   strncpy(file_path, path, 4096);
   strncpy(filename, file, 4096);

   // save all filenames
   R3Grid *affinities[3];
   for (int dim = 0; dim <= 2; ++dim) {
      affinities[dim] = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
      if (!affinities[dim]) return 0;
   }
   R3Grid *human_labels = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
   R3Grid *image = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
   R3Grid *machine_labels = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);

   // read all voxels
   ReadVoxels();

   // go through every voxel
   for (int is = 0; is < NSupervoxels(); ++is) {
      NeuronSupervoxel *supervoxel = Supervoxel(is);
      for (int iv = 0; iv < supervoxel->NVoxels(); ++iv) {
         NeuronVoxel *voxel = supervoxel->LocalVoxel(iv);
         int ix = voxel->XCoordinate();
         int iy = voxel->YCoordinate();
         int iz = voxel->ZCoordinate();

         // set affinities
         for (int dim = 0; dim <= 2; ++dim)
            affinities[dim]->SetGridValue(ix, iy, iz, voxel->affinities[dim]);

         // set human labels
         if (voxel->HumanLabel())
            human_labels->SetGridValue(ix, iy, iz, voxel->HumanLabel()->DataIndex() + 1);
         else
            human_labels->SetGridValue(ix, iy, iz, 0);

         // set image value
         image->SetGridValue(ix, iy, iz, voxel->Image());

         // set machine labels
         if (voxel->IsCellular())
            machine_labels->SetGridValue(ix, iy, iz, voxel->Supervoxel()->DataIndex() + 1);
         else
            machine_labels->SetGridValue(ix, iy, iz, 0);
      }
   }

   // save affinity file
   char affinities_filename[4096];
   sprintf(affinities_filename, "%s/affinities/%s_affinities", path, file);

   // save human label file
   char human_labels_filename[4096];
   sprintf(human_labels_filename, "%s/human_labels/%s_human_labels", path, file);

   // save image
   char image_filename[4096];
   sprintf(image_filename, "%s/images/%s_image", path, file);

   // save machine labels
   char machine_labels_filename[4096];
   sprintf(machine_labels_filename, "%s/machine_labels/%s_machine_labels", path, file);

   RNMeta affinities_meta = RNMeta("Float32", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 3);
   RNMeta human_labels_meta = RNMeta("Uint16", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 1);
   RNMeta image_meta = RNMeta("Float32", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 1);
   RNMeta machine_labels_meta = RNMeta("Int32", resolution[RN_X], resolution[RN_Y], resolution[RN_Z], 1);

   // write all neuron file
   if (NEURON_DEBUG) { printf("Writing affinities file %s...", affinities_filename); fflush(stdout); }
   if (!RNWriteNeuronMetaRawFile(affinities_filename, affinities_meta, affinities)) return 0;
   if (NEURON_DEBUG) printf("done.\n");

   if (NEURON_DEBUG) { printf("Writing human labels file %s...", human_labels_filename); fflush(stdout); }
   if (!RNWriteNeuronMetaRawFile(human_labels_filename, human_labels_meta, human_labels)) return 0;
   if (NEURON_DEBUG) printf("done.\n");

   if (NEURON_DEBUG) { printf("Writing image file %s...", image_filename); fflush(stdout); }
   if (!RNWriteNeuronMetaRawFile(image_filename, image_meta, image)) return 0;
   if (NEURON_DEBUG) printf("done.\n");

   if (NEURON_DEBUG) { printf("Writing machine labels file %s...", machine_labels_filename); fflush(stdout); }
   if (!RNWriteNeuronMetaRawFile(machine_labels_filename, machine_labels_meta, machine_labels)) return 0;
   if (NEURON_DEBUG) printf("done.\n");

   char relative_affinities[4096];
   sprintf(relative_affinities, "affinities/%s_affinities", file);
   char relative_human_labels[4096];
   sprintf(relative_human_labels, "human_labels/%s_human_labels", file);
   char relative_image[4096];
   sprintf(relative_image, "images/%s_image", file);
   char relative_machine_labels[4096];
   sprintf(relative_machine_labels, "machine_labels/%s_machine_labels", file);

   fprintf(fp, "%s.meta\n", relative_affinities);
   fprintf(fp, "%s.meta\n", relative_human_labels);
   fprintf(fp, "%s.meta\n", relative_image);
   fprintf(fp, "%s.meta\n", relative_machine_labels);
   fprintf(fp, "Scale(%f,%f,%f)", scaling[RN_X], scaling[RN_Y], scaling[RN_Z]);

   // release all voxels
   ReleaseVoxels();

   // free memory
   for (int dim = 0; dim <= 2; ++dim)
      delete affinities[dim];
   delete human_labels;
   delete image;
   delete machine_labels;

   // return OK status
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Neuron Read/Write functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
ReadNeuronFile(const char *neuron_filename, RNBoolean read_voxels, RNBoolean read_affinities)
{
   // open file
   FILE *fp = fopen(neuron_filename, "r");
   if (!fp) { fprintf(stderr, "Failed to open file: %s\n", neuron_filename); return 0; }

   // get file path and filename
   char path[4096];
   strncpy(path, neuron_filename, 4096);

   char *file = strchr(path, '/');
   *file = '\0'; file++;

   strncpy(file_path, path, 4096);
   strncpy(filename, file, 4096);

   // convenient variables
   int dummy = 0;

   // read magic keyword
   char magic[16] = { '\0' };
   if (fread(magic, sizeof(char), 16, fp) != (unsigned int)16) {
      fprintf(stderr, "Unable to read first bytes from Neuron data file\n");
      fclose(fp);
      return 0;
   }

   if (strcmp(magic, "Neuron")) {
      fprintf(stderr, "Corrupted file %s\n", neuron_filename);
      fclose(fp);
      return 0;
   }

   if (NEURON_DEBUG) { printf("Reading header..."); fflush(stdout); }
   // read header
   int endian_test;
   int major_version;
   int minor_version;
   unsigned long long data_file_offset = 0;
   fread(&dummy, sizeof(int), 1, fp);
   fread(&dummy, sizeof(int), 1, fp);
   fread(&data_file_offset, sizeof(unsigned long long), 1, fp);
   fread(&endian_test, sizeof(int), 1, fp);
   fread(&major_version, sizeof(int), 1, fp);
   fread(&minor_version, sizeof(int), 1, fp);
   if ((major_version != 0) || (minor_version != 1)) {
      for (int j = 0; j < HEADER_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("done.\n");

   // seek to beginning of data stuff
   rn_assertion(((major_version != 0) || (minor_version != 1)) || (data_file_offset == 0));
   rn_assertion(((major_version == 0) && (minor_version == 1)) || (data_file_offset > 0));
   unsigned long long voxel_file_offset = RNFileTell(fp);
   if (data_file_offset > 0) {
      if (!RNFileSeek(fp, data_file_offset, RN_FILE_SEEK_SET)) {
         fprintf(stderr, "Unable to seek to data file offset\n");
         fclose(fp);
         return 0;
      }
   }

   ///////////////////
   //// Read data ////
   ///////////////////

   // start statistics
   RNTime data_time;
   data_time.Read();

   if (NEURON_DEBUG) { printf("Reading data..."); fflush(stdout); }
   // read data stuff
   int data_nboundaries;
   int data_ncellulars;
   int data_nextracellulars;
   int data_nhuman_labels;
   int data_npredictions;
   int resolution[3];
   RNScalar scaling[3];
   R3Box data_bbox;
   char affinities_meta_filename[4096];
   char human_labels_meta_filename[4096];
   char image_meta_filename[4096];
   char machine_labels_meta_filename[4096];
   fread(&data_nboundaries, sizeof(int), 1, fp);
   fread(&data_ncellulars, sizeof(int), 1, fp);
   fread(&data_nextracellulars, sizeof(int), 1, fp);
   fread(&data_nhuman_labels, sizeof(int), 1, fp);
   fread(&data_npredictions, sizeof(int), 1, fp);
   fread(&(resolution[0]), sizeof(int), 3, fp);
   fread(&(scaling[0]), sizeof(RNScalar), 3, fp);
   fread(&(data_bbox[0][0]), sizeof(RNCoord), 6, fp);
   fread(&(affinities_meta_filename[0]), sizeof(char), 4096, fp);
   fread(&(human_labels_meta_filename[0]), sizeof(char), 4096, fp);
   fread(&(image_meta_filename[0]), sizeof(char), 4096, fp);
   fread(&(machine_labels_meta_filename[0]), sizeof(char), 4096, fp);

   // update resolutions and scalings
   this->resolution[RN_X] = resolution[RN_X];
   this->resolution[RN_Y] = resolution[RN_Y];
   this->resolution[RN_Z] = resolution[RN_Z];
   this->scaling[RN_X] = scaling[RN_X];
   this->scaling[RN_Y] = scaling[RN_Y];
   this->scaling[RN_Z] = scaling[RN_Z];
   this->transformation = R3Affine(R4Matrix(scaling[RN_X], 0, 0, 0, 0, scaling[RN_Y], 0, 0, 0, 0, scaling[RN_Z], 0, 0, 0, 0, 1));
   this->data_size = resolution[RN_X] * resolution[RN_Y] * resolution[RN_Z];
   this->data_sheet_size = resolution[RN_X] * resolution[RN_Y];
   this->data_row_size = resolution[RN_X];
   strncpy(affinities_filename, affinities_meta_filename, 4096);
   strncpy(human_labels_filename, human_labels_meta_filename, 4096);
   strncpy(image_filename, image_meta_filename, 4096);
   strncpy(machine_labels_filename, machine_labels_meta_filename, 4096);

   // read in the voxel mapping
   this->voxel_mapping = new int[NVoxels()];
   if (!voxel_mapping) {
      fprintf(stderr, "Failed to allcoate memory for voxel_mapping\n");
      fclose(fp);
      return 0;
   }
   fread(&voxel_mapping[0], sizeof(int), NVoxels(), fp);
   for (int j = 0; j < DATA_FILLER; ++j) {
      fread(&dummy, sizeof(int), 1, fp);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", data_time.Elapsed());


   /////////////////////////
   //// Read boundaries ////
   /////////////////////////

   // start statistics
   RNTime boundary_time;
   boundary_time.Read();

   if (NEURON_DEBUG) { printf("Reading %d boundaries...", data_nboundaries); fflush(stdout); }
   // reading boundaries
   for (int ib = 0; ib < data_nboundaries; ++ib) {
      NeuronBoundary *boundary = new NeuronBoundary();
      int boundary_naffinities;
      int boundary_data_index;
      int boundary_supervoxel_one_index;
      int boundary_supervoxel_two_index;
      unsigned long long boundary_file_offset;
      RNScalar boundary_maximum;
      RNScalar boundary_mean;
      RNScalar boundary_median;
      RNScalar boundary_minimum;
      RNScalar boundary_stddev;
      RNScalar boundary_skew;
      fread(&boundary_naffinities, sizeof(int), 1, fp);
      fread(&boundary_data_index, sizeof(int), 1, fp);
      fread(&boundary_supervoxel_one_index, sizeof(int), 1, fp);
      fread(&boundary_supervoxel_two_index, sizeof(int), 1, fp);
      fread(&boundary_file_offset, sizeof(unsigned long long), 1, fp);
      fread(&boundary_maximum, sizeof(RNScalar), 1, fp);
      fread(&boundary_mean, sizeof(RNScalar), 1, fp);
      fread(&boundary_median, sizeof(RNScalar), 1, fp);
      fread(&boundary_minimum, sizeof(RNScalar), 1, fp);
      fread(&boundary_stddev, sizeof(RNScalar), 1, fp);
      fread(&boundary_skew, sizeof(RNScalar), 1, fp);
      // read spare room
      for (int j = 0; j < BOUNDARY_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }

      // set current variables
      boundary->supervoxel_one_index = boundary_supervoxel_one_index;
      boundary->supervoxel_two_index = boundary_supervoxel_two_index;
      boundary->SetNAffinities(boundary_naffinities);
      boundary->SetFileOffset(boundary_file_offset);
      boundary->maximum = boundary_maximum;
      boundary->mean = boundary_mean;
      boundary->median = boundary_median;
      boundary->minimum = boundary_minimum;
      boundary->stddev = boundary_stddev;
      boundary->skew = boundary_skew;

      this->InsertBoundary(boundary);

      // make sure the file is not corrupted
      rn_assertion(boundary_data_index == boundary->data_index);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", boundary_time.Elapsed());


   ////////////////////////
   //// Read cellulars ////
   ////////////////////////

   // start statistics
   RNTime cellular_time;
   cellular_time.Read();

   if (NEURON_DEBUG) { printf("Reading %d cellulars...", data_ncellulars); fflush(stdout); }
   // read cellulars
   for (int is = 0; is < data_ncellulars; ++is) {
      NeuronCellular *cellular = new NeuronCellular();
      int cellular_nvoxels;
      int cellular_nboundaries;
      int cellular_data_index;
      RNBoolean cellular_boundary_flag;
      R3Box cellular_bbox;
      unsigned long long cellular_file_offset;
      int cellular_majority_human_label_index;
      int cellular_center_human_label_index;
      int cellular_center_voxel_index;
      fread(&cellular_nvoxels, sizeof(int), 1, fp);
      fread(&cellular_nboundaries, sizeof(int), 1, fp);
      fread(&cellular_data_index, sizeof(int), 1, fp);
      fread(&cellular_boundary_flag, sizeof(RNBoolean), 1, fp);
      fread(&(cellular_bbox[0][0]), sizeof(RNCoord), 6, fp);
      fread(&cellular_file_offset, sizeof(unsigned long long), 1, fp);
      fread(&cellular_majority_human_label_index, sizeof(int), 1, fp);
      fread(&cellular_center_human_label_index, sizeof(int), 1, fp);
      fread(&cellular_center_voxel_index, sizeof(int), 1, fp);
      cellular->SetNVoxels(cellular_nvoxels);
      fread(&(cellular->voxel_mapping[0]), sizeof(int), cellular_nvoxels, fp);
      for (int j = 0; j < CELLULAR_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }

      // create cellular
      cellular->SetBoundaryFlag(cellular_boundary_flag);
      cellular->SetBBox(cellular_bbox);
      cellular->SetFileOffset(cellular_file_offset);
      cellular->majority_human_label_index = cellular_majority_human_label_index;
      cellular->center_human_label_index = cellular_center_human_label_index;
      cellular->center_voxel_index = cellular_center_voxel_index;
      this->InsertCellular(cellular);

      // add boundaries
      for (int ib = 0; ib < cellular_nboundaries; ++ib) {
         int boundary_index;
         fread(&boundary_index, sizeof(int), 1, fp);
         cellular->InsertBoundary(Boundary(boundary_index));
      }
      // make sure the file is not corrupted
      rn_assertion(cellular_data_index == cellular->data_index);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", cellular_time.Elapsed());


   /////////////////////////////
   //// Read extracellulars ////
   /////////////////////////////

   // start statistics
   RNTime extracellular_time;
   extracellular_time.Read();

   if (NEURON_DEBUG) { printf("Reading %d extracellulars...", data_nextracellulars); fflush(stdout); }
   // read extracellulars
   for (int ie = 0; ie < data_nextracellulars; ++ie) {
      NeuronExtracellular *extracellular = new NeuronExtracellular();
      int extracellular_nvoxels;
      int extracellular_nboundaries;
      int extracellular_data_index;
      RNBoolean extracellular_boundary_flag;
      R3Box extracellular_bbox;
      unsigned long long extracellular_file_offset;
      int extracellular_majority_human_label_index;
      int extracellular_center_human_label_index;
      int extracellular_center_voxel_index;
      fread(&extracellular_nvoxels, sizeof(int), 1, fp);
      fread(&extracellular_nboundaries, sizeof(int), 1, fp);
      fread(&extracellular_data_index, sizeof(int), 1, fp);
      fread(&extracellular_boundary_flag, sizeof(RNBoolean), 1, fp);
      fread(&(extracellular_bbox[0][0]), sizeof(RNCoord), 6, fp);
      fread(&extracellular_file_offset, sizeof(unsigned long long), 1, fp);
      fread(&extracellular_majority_human_label_index, sizeof(int), 1, fp);
      fread(&extracellular_center_human_label_index, sizeof(int), 1, fp);
      fread(&extracellular_center_voxel_index, sizeof(int), 1, fp);
      extracellular->SetNVoxels(extracellular_nvoxels);
      fread(&(extracellular->voxel_mapping[0]), sizeof(int), extracellular_nvoxels, fp);
      for (int j = 0; j < EXTRACELLULAR_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }
      // create extracellular
      extracellular->SetBoundaryFlag(extracellular_boundary_flag);
      extracellular->SetBBox(extracellular_bbox);
      extracellular->SetFileOffset(extracellular_file_offset);
      extracellular->majority_human_label_index = extracellular_majority_human_label_index;
      extracellular->center_human_label_index = extracellular_center_human_label_index;
      extracellular->center_voxel_index = extracellular_center_voxel_index;
      InsertExtracellular(extracellular);
      // add boundaries
      for (int ib = 0; ib < extracellular_nboundaries; ++ib) {
         int boundary_index;
         fread(&boundary_index, sizeof(int), 1, fp);
         extracellular->InsertBoundary(Boundary(boundary_index));
      }

      // make sure the file is not corrupted 
      rn_assertion(extracellular_data_index == extracellular->data_index);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", extracellular_time.Elapsed());


   ///////////////////////////
   //// Read human labels ////
   ///////////////////////////

   // start statistics
   RNTime human_label_time;
   human_label_time.Read();

   if (NEURON_DEBUG) { printf("Reading %d human labels...", data_nhuman_labels); fflush(stdout); }
   // read truths
   for (int it = 0; it < data_nhuman_labels; ++it) {
      // create NeuronHumanLabel
      NeuronHumanLabel *human_label = new NeuronHumanLabel();
      int human_label_nvoxels;
      int human_label_nsupervoxels;
      int human_label_data_index;
      R3Box human_label_bbox;
      fread(&human_label_nvoxels, sizeof(int), 1, fp);
      fread(&human_label_nsupervoxels, sizeof(int), 1, fp);
      fread(&human_label_data_index, sizeof(int), 1, fp);
      fread(&(human_label_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int iv = 0; iv < human_label_nsupervoxels; ++iv) {
         int supervoxel_index;
         fread(&supervoxel_index, sizeof(int), 1, fp);
         // add supervoxel
         human_label->InsertSupervoxel(Supervoxel(supervoxel_index));
      }
      human_label->SetNVoxels(human_label_nvoxels);
      fread(&(human_label->voxel_mapping[0]), sizeof(int), human_label_nvoxels, fp);
      for (int j = 0; j < HUMAN_LABEL_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }
      human_label->SetBBox(human_label_bbox);
      InsertHumanLabel(human_label);

      // make sure the file is not corrupted
      rn_assertion(human_label_data_index == human_label->data_index);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", human_label_time.Elapsed());


   //////////////////////////
   //// Read predictions ////
   //////////////////////////

   // start statistics
   RNTime prediction_time;
   prediction_time.Read();

   if (NEURON_DEBUG) { printf("Reading %d predictions...", data_npredictions); fflush(stdout); }
   // read truths
   for (int it = 0; it < data_npredictions; ++it) {
      // create NeuronHumanLabel
      NeuronPrediction *prediction = new NeuronPrediction();
      int prediction_nvoxels;
      int prediction_nsupervoxels;
      int prediction_data_index;
      R3Box prediction_bbox;
      fread(&prediction_nvoxels, sizeof(int), 1, fp);
      fread(&prediction_nsupervoxels, sizeof(int), 1, fp);
      fread(&prediction_data_index, sizeof(int), 1, fp);
      fread(&(prediction_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int iv = 0; iv < prediction_nsupervoxels; ++iv) {
         int supervoxel_index;
         fread(&supervoxel_index, sizeof(int), 1, fp);
         // add supervoxel
         prediction->InsertSupervoxel(Supervoxel(supervoxel_index));
      }
      prediction->SetNVoxels(prediction_nvoxels);
      fread(&(prediction->voxel_mapping[0]), sizeof(int), prediction_nvoxels, fp);
      for (int j = 0; j < HUMAN_LABEL_FILLER; ++j) {
         fread(&dummy, sizeof(int), 1, fp);
      }
      prediction->SetBBox(prediction_bbox);
      InsertPrediction(prediction);

      // make sure the file is not corrupted
      rn_assertion(prediction_data_index == prediction->data_index);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", prediction_time.Elapsed());


   ////////////////////////////////////
   //// Read voxels and affinities ////
   ////////////////////////////////////

   // read voxels
   if (read_voxels) {
      // start statistics
      RNTime voxel_time;
      voxel_time.Read();

      if (NEURON_DEBUG) { printf("Reading voxels..."); fflush(stdout); }
      if ((major_version != 0) || (minor_version != 1)) {
         if (!RNFileSeek(fp, voxel_file_offset, RN_FILE_SEEK_SET)) {
            fprintf(stderr, "Unable to seek to voxel file offset\n");
            fclose(fp);
            return 0;
         }
         for (int is = 0; is < NSupervoxels(); ++is) {
            NeuronSupervoxel *supervoxel = Supervoxel(is);
            if (!supervoxel->ReadVoxels(fp, TRUE)) return 0;
         }
      }

      ++read_voxel_count;

      // print statistics
      if (NEURON_DEBUG) printf("in %.2f seconds.\n", voxel_time.Elapsed());
   }

   // read affinities
   if (read_affinities) {
      // start statistics
      RNTime affinity_time;
      affinity_time.Read();

      if (NEURON_DEBUG) { printf("Reading affinities..."); fflush(stdout); }
      for (int ib = 0; ib < NBoundaries(); ++ib) {
         NeuronBoundary *boundary = Boundary(ib);
         if (!boundary->ReadAffinities(fp, TRUE)) return 0;
      }

      ++read_affinities;

      // print statistics
      if (NEURON_DEBUG) printf("in %.2f seconds.\n\n", affinity_time.Elapsed());
   }

   // set data stuff
   SetBBox(data_bbox);

   // close file
   fclose(fp);

   // return OK status
   return 1;
}



int NeuronData::
WriteNeuronFile(const char *neuron_filename)
{
   if (NEURON_DEBUG) printf("Updating bounding box...");
   // update all bounding boxes
   UpdateBBox();
   if (NEURON_DEBUG) printf("...done.\n");

   // convenient variables
   int dummy = 0;

   // make sure everything is uptodate
   R3Box foobar = GridBox();
   if (foobar[0][0] == FLT_MAX) printf("Bounding Box not defined for this data");

   // open file
   FILE *fp = fopen(neuron_filename, "wb");
   if (!fp) {
      fprintf(stderr, "Unable to open Neuron data file: %s\n", filename);
      return 0;
   }

   // write magic keyword
   char magic[16] = { '\0' };
   strncpy(magic, "Neuron", 16);
   if (fwrite(magic, sizeof(char), 16, fp) != (unsigned int)16) {
      fprintf(stderr, "Unable to write Neuron data file: %s\n", filename);
      return 0;
   }


   ///////////////////////
   //// Write headers ////
   ///////////////////////

   if (NEURON_DEBUG) { printf("Writing header..."); fflush(stdout); }
   // write header
   int endian_test = 1;
   int major_version = 0;
   int minor_version = 2;
   unsigned long long data_file_offset = 0;
   fwrite(&dummy, sizeof(int), 1, fp);
   fwrite(&dummy, sizeof(int), 1, fp);
   fwrite(&data_file_offset, sizeof(unsigned long long), 1, fp);
   fwrite(&endian_test, sizeof(int), 1, fp);
   fwrite(&major_version, sizeof(int), 1, fp);
   fwrite(&minor_version, sizeof(int), 1, fp);
   for (int j = 0; j < HEADER_FILLER; ++j) {
      fwrite(&dummy, sizeof(int), 1, fp);
   }
   // print statistics
   if (NEURON_DEBUG) printf("done.\n");


   //////////////////////
   //// Write voxels ////
   //////////////////////

   // start statistics
   RNTime voxel_time;
   voxel_time.Read();

   if (NEURON_DEBUG) { printf("Writing voxels..."); fflush(stdout); }
   // write voxels (through supervoxels)
   for (int is = 0; is < NSupervoxels(); ++is) {
      NeuronSupervoxel *supervoxel = Supervoxel(is);
      if (!supervoxel->ReadVoxels()) return 0;
      if (!supervoxel->WriteVoxels(fp, FALSE)) return 0;
      if (!supervoxel->ReleaseVoxels()) return 0;
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", voxel_time.Elapsed());


   //////////////////////////
   //// Write affinities ////
   //////////////////////////

   // start statistics
   RNTime affinity_time;
   affinity_time.Read();

   if (NEURON_DEBUG) { printf("Writing affinities..."); fflush(stdout); }
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      NeuronBoundary *boundary = Boundary(ib);
      if (!boundary->ReadAffinities()) return 0;
      if (!boundary->WriteAffinities(fp, FALSE)) return 0;
      if (!boundary->ReleaseAffinities()) return 0;
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", affinity_time.Elapsed());

   // get data file offset
   data_file_offset = RNFileTell(fp);


   ////////////////////
   //// Write data ////
   ////////////////////

   // start statistics
   RNTime data_time;
   data_time.Read();

   if (NEURON_DEBUG) { printf("Writing data..."); fflush(stdout); }
   // write data
   int data_nboundaries = NBoundaries();
   int data_ncellulars = NCellulars();
   int data_nextracellulars = NExtracellulars();
   int data_nhuman_labels = NHumanLabels();
   int data_npredictions = NPredictions();
   R3Box data_bbox = GridBox();
   fwrite(&data_nboundaries, sizeof(int), 1, fp);
   fwrite(&data_ncellulars, sizeof(int), 1, fp);
   fwrite(&data_nextracellulars, sizeof(int), 1, fp);
   fwrite(&data_nhuman_labels, sizeof(int), 1, fp);
   fwrite(&data_npredictions, sizeof(int), 1, fp);
   fwrite(&(resolution[0]), sizeof(int), 3, fp);
   fwrite(&(scaling[0]), sizeof(RNScalar), 3, fp);
   fwrite(&(data_bbox[0][0]), sizeof(RNCoord), 6, fp);
   fwrite(&(affinities_filename[0]), sizeof(char), 4096, fp);
   fwrite(&(human_labels_filename[0]), sizeof(char), 4096, fp);
   fwrite(&(image_filename[0]), sizeof(char), 4096, fp);
   fwrite(&(machine_labels_filename[0]), sizeof(char), 4096, fp);
   fwrite(&(voxel_mapping[0]), sizeof(int), NVoxels(), fp);
   for (int j = 0; j < DATA_FILLER; ++j) {
      fwrite(&dummy, sizeof(int), 1, fp);
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", data_time.Elapsed());

   /////////////////////////
   //// Write bondaries ////
   /////////////////////////

   // start statistics
   RNTime boundary_time;
   boundary_time.Read();

   if (NEURON_DEBUG) { printf("Writing %d boundaries...", data_nboundaries); fflush(stdout); }
   for (int ib = 0; ib < NBoundaries(); ++ib) {
      NeuronBoundary *boundary = Boundary(ib);
      int boundary_naffinities = boundary->NAffinities();
      int boundary_data_index = boundary->DataIndex();
      int boundary_supervoxel_one_index = boundary->SupervoxelOneIndex();
      int boundary_supervoxel_two_index = boundary->SupervoxelTwoIndex();
      unsigned long long boundary_file_offset = boundary->FileOffset();
      RNScalar boundary_maximum = boundary->Maximum();
      RNScalar boundary_mean = boundary->Mean();
      RNScalar boundary_median = boundary->Median();
      RNScalar boundary_minimum = boundary->Minimum();
      RNScalar boundary_stddev = boundary->StdDev();
      RNScalar boundary_skew = boundary->Skew();
      fwrite(&boundary_naffinities, sizeof(int), 1, fp);
      fwrite(&boundary_data_index, sizeof(int), 1, fp);
      fwrite(&boundary_supervoxel_one_index, sizeof(int), 1, fp);
      fwrite(&boundary_supervoxel_two_index, sizeof(int), 1, fp);
      fwrite(&boundary_file_offset, sizeof(unsigned long long), 1, fp);
      fwrite(&boundary_maximum, sizeof(RNScalar), 1, fp);
      fwrite(&boundary_mean, sizeof(RNScalar), 1, fp);
      fwrite(&boundary_median, sizeof(RNScalar), 1, fp);
      fwrite(&boundary_minimum, sizeof(RNScalar), 1, fp);
      fwrite(&boundary_stddev, sizeof(RNScalar), 1, fp);
      fwrite(&boundary_skew, sizeof(RNScalar), 1, fp);
      // write spare room
      for (int j = 0; j < BOUNDARY_FILLER; ++j) {
         fwrite(&dummy, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", boundary_time.Elapsed());


   /////////////////////////
   //// Write cellulars ////
   /////////////////////////

   // start statistics
   RNTime cellular_time;
   cellular_time.Read();

   if (NEURON_DEBUG) { printf("Writing %d cellulars...", data_ncellulars); fflush(stdout); }
   // write supervoxels
   for (int is = 0; is < NCellulars(); ++is) {
      NeuronCellular *cellular = Cellular(is);
      int cellular_nvoxels = cellular->NVoxels();
      int cellular_nboundaries = cellular->NBoundaries();
      int cellular_data_index = cellular->DataIndex();
      RNBoolean cellular_boundary_flag = cellular->IsOnBoundary();
      R3Box cellular_bbox = cellular->BBox();
      unsigned long long cellular_file_offset = cellular->FileOffset();
      int cellular_majority_human_label_index = cellular->majority_human_label_index;
      int cellular_center_human_label_index = cellular->center_human_label_index;
      int cellular_center_voxel_index = cellular->center_voxel_index;
      fwrite(&cellular_nvoxels, sizeof(int), 1, fp);
      fwrite(&cellular_nboundaries, sizeof(int), 1, fp);
      fwrite(&cellular_data_index, sizeof(int), 1, fp);
      fwrite(&cellular_boundary_flag, sizeof(RNBoolean), 1, fp);
      fwrite(&(cellular_bbox[0][0]), sizeof(RNCoord), 6, fp);
      fwrite(&cellular_file_offset, sizeof(unsigned long long), 1, fp);
      fwrite(&cellular_majority_human_label_index, sizeof(int), 1, fp);
      fwrite(&cellular_center_human_label_index, sizeof(int), 1, fp);
      fwrite(&cellular_center_voxel_index, sizeof(int), 1, fp);
      fwrite(&(cellular->voxel_mapping[0]), sizeof(int), cellular_nvoxels, fp);
      for (int j = 0; j < CELLULAR_FILLER; ++j) {
         fwrite(&dummy, sizeof(int), 1, fp);
      }
      for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = cellular->Boundary(ib);
         int boundary_index = boundary->data_index;
         fwrite(&boundary_index, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", cellular_time.Elapsed());


   //////////////////////////////
   //// Write extracellulars ////
   //////////////////////////////

   // start statistics
   RNTime extracellular_time;
   extracellular_time.Read();


   if (NEURON_DEBUG) { printf("Writing %d extracellulars...", data_nextracellulars); fflush(stdout); }
   // write extracellulars
   for (int ie = 0; ie < NExtracellulars(); ++ie) {
      NeuronExtracellular *extracellular = Extracellular(ie);
      int extracellular_nvoxels = extracellular->NVoxels();
      int extracellular_nboundaries = extracellular->NBoundaries();
      int extracellular_data_index = extracellular->DataIndex();
      RNBoolean extracellular_boundary_flag = extracellular->IsOnBoundary();
      R3Box extracellular_bbox = extracellular->BBox();
      unsigned long long extracellular_file_offset = extracellular->FileOffset();
      int extracellular_majority_human_label_index = extracellular->majority_human_label_index;
      int extracellular_center_human_label_index = extracellular->center_human_label_index;
      int extracellular_center_voxel_index = extracellular->center_voxel_index;
      fwrite(&extracellular_nvoxels, sizeof(int), 1, fp);
      fwrite(&extracellular_nboundaries, sizeof(int), 1, fp);
      fwrite(&extracellular_data_index, sizeof(int), 1, fp);
      fwrite(&extracellular_boundary_flag, sizeof(RNBoolean), 1, fp);
      fwrite(&(extracellular_bbox[0][0]), sizeof(RNCoord), 6, fp);
      fwrite(&extracellular_file_offset, sizeof(unsigned long long), 1, fp);
      fwrite(&extracellular_majority_human_label_index, sizeof(int), 1, fp);
      fwrite(&extracellular_center_human_label_index, sizeof(int), 1, fp);
      fwrite(&extracellular_center_voxel_index, sizeof(int), 1, fp);
      fwrite(&(extracellular->voxel_mapping[0]), sizeof(int), extracellular_nvoxels, fp);
      for (int j = 0; j < EXTRACELLULAR_FILLER; ++j) {
         fwrite(&dummy, sizeof(int), 1, fp);
      }
      for (int ib = 0; ib < extracellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = extracellular->Boundary(ib);
         int boundary_index = boundary->data_index;
         fwrite(&boundary_index, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", extracellular_time.Elapsed());


   ////////////////////////////
   //// Write human labels ////
   ////////////////////////////

   // start statistics
   RNTime human_labels_time;
   human_labels_time.Read();

   if (NEURON_DEBUG) { printf("Writing %d human labels...", data_nhuman_labels); fflush(stdout); }
   // write truths
   for (int it = 0; it < NHumanLabels(); ++it) {
      NeuronHumanLabel *human_label = HumanLabel(it);
      int human_label_nvoxels = human_label->NVoxels();
      int human_label_nsupervoxels = human_label->NSupervoxels();
      int human_label_data_index = human_label->DataIndex();
      R3Box human_label_bbox = human_label->BBox();
      fwrite(&human_label_nvoxels, sizeof(int), 1, fp);
      fwrite(&human_label_nsupervoxels, sizeof(int), 1, fp);
      fwrite(&human_label_data_index, sizeof(int), 1, fp);
      fwrite(&(human_label_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int iv = 0; iv < human_label_nsupervoxels; ++iv) {
         int supervoxel_index = human_label->Supervoxel(iv)->SupervoxelIndex();
         fwrite(&supervoxel_index, sizeof(int), 1, fp);
      }
      fwrite(&(human_label->voxel_mapping[0]), sizeof(int), human_label_nvoxels, fp);
      for (int j = 0; j < HUMAN_LABEL_FILLER; ++j) {
         fwrite(&dummy, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", human_labels_time.Elapsed());


   ///////////////////////////
   //// Write predictions ////
   ///////////////////////////
   /* TODO allow saving of boundaries of predictions/human labels */
   // start statistics
   RNTime prediction_time;
   prediction_time.Read();

   if (NEURON_DEBUG) { printf("Writing %d predictions...", data_npredictions); fflush(stdout); }
   // write segments
   for (int is = 0; is < NPredictions(); ++is) {
      NeuronPrediction *prediction = Prediction(is);
      int prediction_nvoxels = prediction->NVoxels();
      int prediction_nsupervoxels = prediction->NSupervoxels();
      int prediction_data_index = prediction->DataIndex();
      R3Box prediction_bbox = prediction->BBox();
      fwrite(&prediction_nvoxels, sizeof(int), 1, fp);
      fwrite(&prediction_nsupervoxels, sizeof(int), 1, fp);
      fwrite(&prediction_data_index, sizeof(int), 1, fp);
      fwrite(&(prediction_bbox[0][0]), sizeof(RNCoord), 6, fp);
      for (int iv = 0; iv < prediction_nsupervoxels; ++iv) {
         int supervoxel_index = prediction->Supervoxel(iv)->SupervoxelIndex();
         fwrite(&supervoxel_index, sizeof(int), 1, fp);
      }
      fwrite(&(prediction->voxel_mapping[0]), sizeof(int), prediction_nvoxels, fp);
      for (int j = 0; j < PREDICTION_FILLER; ++j) {
         fwrite(&dummy, sizeof(int), 1, fp);
      }
   }
   // print statistics
   if (NEURON_DEBUG) printf("in %.2f seconds.\n", prediction_time.Elapsed());


   ////////////////////////
   //// Rewrite header ////
   ////////////////////////

   if (NEURON_DEBUG) { printf("Rewriting header..."); fflush(stdout); }
   // re-write header with updated data file offset
   RNFileSeek(fp, 0, RN_FILE_SEEK_SET);
   fwrite(magic, sizeof(char), 16, fp);
   fwrite(&dummy, sizeof(int), 1, fp);
   fwrite(&dummy, sizeof(int), 1, fp);
   fwrite(&data_file_offset, sizeof(unsigned long long), 1, fp);
   fwrite(&endian_test, sizeof(int), 1, fp);
   fwrite(&major_version, sizeof(int), 1, fp);
   fwrite(&minor_version, sizeof(int), 1, fp);
   if (NEURON_DEBUG) printf("done.\n");

   // close file
   fclose(fp);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int NeuronData::
ReadAffinities(R3Grid *affinities_grids[3]) const
{
   for (int dim = 0; dim < 3; ++dim)
      rn_assertion(affinities_grids[dim] != NULL);

   // make sure voxels are resident
   if (!AreVoxelsResident()) {
      char input_filename[4096];
      sprintf(input_filename, "%s/%s", file_path, affinities_filename);
      R3Grid **affinities = RNReadNeuronMetaRawFile(input_filename, TRUE);
      for (int dim = 0; dim <= 2; ++dim) {
         for (int ix = 0; ix < resolution[RN_X]; ++ix) {
            for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
               for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
                  affinities_grids[dim]->SetGridValue(ix, iy, iz, affinities[dim]->GridValue(ix, iy, iz));
               }
            }
         }
      }
      // delete affinities
      for (int dim = 0; dim <= 2; ++dim)
         delete affinities[dim];
      return 1;
   }

   // go through every voxel
   for (int dim = 0; dim <= 2; ++dim) {
      for (int ix = 0; ix < XResolution(); ++ix) {
         for (int iy = 0; iy < YResolution(); ++iy) {
            for (int iz = 0; iz < ZResolution(); ++iz) {
               NeuronVoxel *voxel = Voxel(ix, iy, iz);

               // set affinity
               affinities_grids[dim]->SetGridValue(ix, iy, iz, voxel->affinities[dim]);
            }
         }
      }
   }

   // return success
   return 1;
}



R3Grid *NeuronData::
ReadTruthTxtFile(const char human_labels_txt_filename[4096]) const
{
   // read the txt file
   FILE *truth_fp = fopen(human_labels_txt_filename, "r");
   if (!truth_fp) { fprintf(stderr, "Failed to read %s\n", human_labels_txt_filename); return NULL; }

   // create an R3Grid
   R3Grid *human_labels_grid = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);

   // initialize values to 0
   for (int ix = 0; ix < resolution[RN_X]; ++ix) {
      for (int iy = 0; iy < resolution[RN_Y]; ++iy) {
         for (int iz = 0; iz < resolution[RN_Z]; ++iz) {
            human_labels_grid->SetGridValue(ix, iy, iz, 0);
         }
      }
   }


   // read every line of the file for each supervoxel
   char input_line[4096];
   while (fgets(input_line, 4096, truth_fp) != NULL) {
      int cellular_index;
      RNScalar agreement;
      sscanf(input_line, "%d %lf", &cellular_index, &agreement);

      // every voxel in this cellular gets labeled as truth
      NeuronCellular *cellular = Cellular(cellular_index);
      for (int iv = 0; iv < cellular->NVoxels(); ++iv) {
         NeuronVoxel *voxel = cellular->LocalVoxel(iv);
         human_labels_grid->SetGridValue(voxel->XCoordinate(), voxel->YCoordinate(), voxel->ZCoordinate(), 1);
      }
   }

   // close file
   fclose(truth_fp);

   // return grid
   return human_labels_grid;
}



int NeuronData::
ReadHumanLabels(R3Grid *human_label_grid) const
{
   // make sure voxels are resident
   rn_assertion(human_label_grid != NULL);

   // make sure voxels are resident
   if (!AreVoxelsResident()) {
      char input_filename[4096];
      sprintf(input_filename, "%s/%s", file_path, human_labels_filename);
      *human_label_grid = *(RNReadNeuronMetaRawFile(input_filename));
      return 1;
   }

   // go through every voxel
   for (int ix = 0; ix < XResolution(); ++ix) {
      for (int iy = 0; iy < YResolution(); ++iy) {
         for (int iz = 0; iz < ZResolution(); ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);

            // set human label grid
            if (voxel->HumanLabel())
               human_label_grid->SetGridValue(ix, iy, iz, voxel->HumanLabel()->DataIndex() + 1);
            else
               human_label_grid->SetGridValue(ix, iy, iz, 0);
         }
      }
   }

   // return success
   return 1;
}



int NeuronData::
ReadImages(R3Grid *image_grid) const
{
   // make sure voxels are resident
   rn_assertion(image_grid != NULL);

   // make sure voxels are resident
   if (!AreVoxelsResident()) {
      char input_filename[4096];
      sprintf(input_filename, "%s/%s", file_path, image_filename);
      *image_grid = *(RNReadNeuronMetaRawFile(input_filename));
      return 1;
   }

   // go through every voxel
   for (int ix = 0; ix < XResolution(); ++ix) {
      for (int iy = 0; iy < YResolution(); ++iy) {
         for (int iz = 0; iz < ZResolution(); ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);

            // set image grid
            image_grid->SetGridValue(ix, iy, iz, voxel->Image());
         }
      }
   }

   // return success
   return 1;
}



int NeuronData::
ReadMachineLabels(R3Grid *machine_labels_grid) const
{
   // make sure voxels are resident
   rn_assertion(machine_labels_grid != NULL);

   // make sure voxels are resident
   if (!AreVoxelsResident()) {
      char input_filename[4096];
      sprintf(input_filename, "%s/%s", file_path, machine_labels_filename);
      *machine_labels_grid = *(RNReadNeuronMetaRawFile(input_filename));
      return 1;
   }

   // go through every voxel
   for (int ix = 0; ix < XResolution(); ++ix) {
      for (int iy = 0; iy < YResolution(); ++iy) {
         for (int iz = 0; iz < ZResolution(); ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);

            // set machine label grid value
            if (voxel->IsCellular())
               machine_labels_grid->SetGridValue(ix, iy, iz, voxel->Supervoxel()->DataIndex() + 1);
            else
               machine_labels_grid->SetGridValue(ix, iy, iz, 0);
         }
      }
   }

   // return success
   return 1;
}



int NeuronData::
ReadPredictions(R3Grid *predictions_grid) const
{
   // make sure voxels are resident
   rn_assertion(predictions_grid != NULL);

   // go through every voxel
   for (int ix = 0; ix < XResolution(); ++ix) {
      for (int iy = 0; iy < YResolution(); ++iy) {
         for (int iz = 0; iz < ZResolution(); ++iz) {
            NeuronVoxel *voxel = Voxel(ix, iy, iz);

            // set machine label grid value
            if (voxel->Prediction())
               predictions_grid->SetGridValue(ix, iy, iz, voxel->Prediction()->DataIndex() + 1);
            else
               predictions_grid->SetGridValue(ix, iy, iz, 0);
         }
      }
   }

   // return success
   return 1;
}



int NeuronData::
ReadPredictions(const char meta_filename[4096])
{
   // get root filename
   char root_filename[4096];
   strncpy(root_filename, meta_filename, 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // get the grid for the predictions and convert into int array
   R3Grid *prediction_grid = RNReadNeuronMetaRawFile(root_filename);
   int *predictions_array = new int[NVoxels()];
   printf("Here\n");
   for (int iv = 0; iv < NVoxels(); ++iv) {
      int proposal_value = (int)(prediction_grid->GridValue(iv) + 0.5);
      predictions_array[iv] = proposal_value;
   }

   // free memory
   delete prediction_grid;

   // deflate the proposals
   RNDeflateIntegerArray(predictions_array, NVoxels());

   R3Grid *predictions = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
   for (int iv = 0; iv < NVoxels(); ++iv) {
      predictions->SetGridValue(iv, predictions_array[iv]);
   }

   // call the actual function to change data structure
   CreatePredictions(predictions);

   // free memory
   delete[] predictions_array;

   // return success
   return 1;
}



int NeuronData::
ReadVoxels(void)
{
   // check read counter
   if (++read_voxel_count > 1) return 1;

   // read all voxels
   for (int is = 0; is < NSupervoxels(); ++is) {
      Supervoxel(is)->ReadVoxels();
   }

   // return success
   return 1;
}



int NeuronData::
ReleaseVoxels(void)
{
   // check read counter
   if (--read_voxel_count > 0) return 1;

   // release all voxels
   for (int is = 0; is < NSupervoxels(); ++is)
      Supervoxel(is)->ReleaseVoxels();

   // return success
   return 1;
}
