// Source file for the path statistics algorithm



// include files 


#include <algorithm>
#include <vector>
#include "Neuron/Neuron.h"

#define _DEBUG TRUE
#define min(a,b) (((a) < (b)) ? (a) : (b))

// useful constants

static int nbins = 100;



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static RNScalar scaling[3];
static RNScalar max_affinity = 0.80;
static int cdf = 0;
static int pdf = 0;
static int trace = 0;



// global variables

static NeuronData *nd = NULL;
static char *input_filename = NULL;
static char root_filename[4096];
static char root_input_filename[4096];
static char root_output_filename[4096];



// path variables

static RNBoolean *voxel_last_correct = NULL;
static int *voxel_ncorrect = NULL;
static int *voxel_path_length = NULL;



// dijkstra input

static RNScalar *dijkstra_distances = NULL;
static int *dijkstra_prev_node = NULL;


////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate neuron data\n");
      return 0;
   }

   nd->ReadFile(input_filename, TRUE);

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }
   if (print_debug) {
      printf("  Bounding Box: (%0.2f, %0.2f, %0.2f) to (%0.2f, %0.2f, %0.2f)\n", nd->WorldBox().XMin(), nd->WorldBox().YMin(), nd->WorldBox().ZMin(), nd->WorldBox().XMax(), nd->WorldBox().YMax(), nd->WorldBox().ZMax());
      printf("  Voxels: %d\n", nd->NVoxels());
      printf("  Supervoxels: %d\n", nd->NSupervoxels());
      printf("  Extracellulars: %d\n", nd->NExtracellulars());
      printf("  Boundaries: %d\n", nd->NBoundaries());
      printf("  Human Labels: %d\n", nd->NHumanLabels());
      printf("  Predictions: %d\n", nd->NPredictions());
      fflush(stdout);
   }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Path Update functions
////////////////////////////////////////////////////////////////////////

static void
UpdatePath(int voxel_index)
{
   // get truth for this voxel
   NeuronVoxel *voxel = nd->Voxel(voxel_index);

   // skip extracellulars
   if (voxel->HumanLabel() == NULL) return;

   // get the entire path
   std::vector<int> voxel_path = std::vector<int>();
   int previous = voxel_index;
   int index;
   do {
      index = previous;
      previous = dijkstra_prev_node[index];
      voxel_path.push_back(index);
   } while (previous != index);

   // get path length for this voxel
   voxel_path_length[voxel_index] = (int)voxel_path.size();

   // find the number of matching along the path
   for (int iv = 0; iv < voxel_path_length[voxel_index]; ++iv) {
      int path_index = voxel_path[iv];
      NeuronVoxel *path_voxel = nd->Voxel(path_index);

      // increment for correct truth
      if (voxel->HumanLabel() == path_voxel->HumanLabel())
         voxel_ncorrect[voxel_index]++;
   }

   // see if voxel last correct is correct
   NeuronHumanLabel *human_label = voxel->HumanLabel();
   int last_path_index = voxel_path[voxel_path.size() - 1];
   NeuronHumanLabel *last_human_label = nd->Voxel(last_path_index)->HumanLabel();
   if (human_label == last_human_label) {
      voxel_last_correct[voxel_index] = TRUE;
   }

   // just checking...
   NeuronVoxel *last_path = nd->Voxel(last_path_index);
   rn_assertion(last_path->IsOnBoundary());
}



static void
UpdatePathVariables(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // go through each voxel and update the path
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel() == NULL) continue;
      UpdatePath(iv);
   }

   // print statistics
   if (print_verbose) printf("Read all paths in %.2f seconds\n", start_time.Elapsed());
}



///////////////////////////////////////////////////////////////////////
// Struct definitions
///////////////////////////////////////////////////////////////////////

struct dijkstra_distance_metric {
   RNScalar distance;
   int index;
};



struct dijkstra_path_length_metric {
   RNScalar path_length;
   int index;
};



struct voxel_distance_metric {
   int voxel_distance;
   int index;
};



///////////////////////////////////////////////////////////////////////
// Struct comparison functions
///////////////////////////////////////////////////////////////////////

bool
distance_compare(dijkstra_distance_metric i, dijkstra_distance_metric j) {
   return (i.distance < j.distance);
}



bool
path_length_compare(dijkstra_path_length_metric i, dijkstra_path_length_metric j) {
   return (i.path_length < j.path_length);
}



bool
voxel_distance_compare(voxel_distance_metric i, voxel_distance_metric j) {
   return (i.voxel_distance < j.voxel_distance);
}



///////////////////////////////////////////////////////////////////////
// Distribution Output functions
///////////////////////////////////////////////////////////////////////

static int
DijkstraDistanceCDF(RNScalar min_distance, RNScalar max_distance, std::vector<dijkstra_distance_metric>& distances)
{
   printf("Starting dijkstra distance cdf..."); fflush(stdout);

   // get distance cdf filename
   char distance_cdf_filename[4096];
   sprintf(distance_cdf_filename, "%s_dijkstra_distance_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(distance_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", distance_cdf_filename); return 0; }

   // get distance range and bin size
   RNScalar distance_range = max_distance - min_distance;
   // add a small factor for double precision errors
   RNScalar distance_range_bin_size = distance_range / nbins + 10e-6;

   RNScalar next_distance = distance_range_bin_size + min_distance;
   RNScalar prev_distance = min_distance;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;

   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].distance < next_distance) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_distance, ((RNScalar)ncorrect) / nnodes);

      prev_distance = next_distance;
      next_distance = next_distance + distance_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraDistanceCDFBoundary(RNScalar min_distance, RNScalar max_distance, std::vector<dijkstra_distance_metric>& distances)
{
   printf("Starting dijkstra distance cdf for boundaries..."); fflush(stdout);

   // get distance cdf filename
   char distance_cdf_filename[4096];
   sprintf(distance_cdf_filename, "%s_dijkstra_distance_boundary_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(distance_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", distance_cdf_filename); return 0; }

   // get distance range and bin size
   RNScalar distance_range = max_distance - min_distance;
   // add a small factor for double precision errors
   RNScalar distance_range_bin_size = distance_range / nbins + 10e-6;

   RNScalar next_distance = distance_range_bin_size + min_distance;
   RNScalar prev_distance = min_distance;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;

   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].distance < next_distance) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index])
            ncorrect++;
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_distance, ((RNScalar)ncorrect) / nnodes);

      prev_distance = next_distance;
      next_distance = next_distance + distance_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraDistancePDF(RNScalar min_distance, RNScalar max_distance, std::vector<dijkstra_distance_metric>& distances)
{
   printf("Starting dijkstra distance pdf..."); fflush(stdout);

   // get distance pdf filename
   char distance_pdf_filename[4096];
   sprintf(distance_pdf_filename, "%s_dijkstra_distance_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(distance_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", distance_pdf_filename); return 0; }

   // get distance range and bin size
   RNScalar distance_range = max_distance - min_distance;
   // add a small factor for double precision errors
   RNScalar distance_range_bin_size = distance_range / nbins + 10e-6;

   RNScalar next_distance = distance_range_bin_size + min_distance;
   RNScalar prev_distance = min_distance;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].distance < next_distance) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_distance, ((RNScalar)ncorrect) / nnodes);

      prev_distance = next_distance;
      next_distance = next_distance + distance_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraDistancePDFBoundary(RNScalar min_distance, RNScalar max_distance, std::vector<dijkstra_distance_metric>& distances)
{
   printf("Starting dijkstra distance pdf for boundaries..."); fflush(stdout);

   // get distance pdf filename
   char distance_pdf_filename[4096];
   sprintf(distance_pdf_filename, "%s_dijkstra_distance_boundary_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(distance_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", distance_pdf_filename); return 0; }

   // get distance range and bin size
   RNScalar distance_range = max_distance - min_distance;
   // add a small factor for double precision errors
   RNScalar distance_range_bin_size = distance_range / nbins + 10e-6;

   RNScalar next_distance = distance_range_bin_size + min_distance;
   RNScalar prev_distance = min_distance;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].distance < next_distance) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index])
            ncorrect++;
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_distance, ((RNScalar)ncorrect) / nnodes);

      prev_distance = next_distance;
      next_distance = next_distance + distance_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraPathLengthCDF(int min_path_length, int max_path_length, std::vector<dijkstra_path_length_metric>& path_lengths)
{
   printf("Starting dijkstra path length cdf..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_dijkstra_path_length_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < path_lengths.size() && path_lengths[iv].path_length < next_path_length) {
         int index = path_lengths[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraPathLengthCDFBoundary(int min_path_length, int max_path_length, std::vector<dijkstra_path_length_metric>& path_lengths)
{
   printf("Starting dijkstra path length cdf for boundaries..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_dijkstra_path_length_boundary_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < path_lengths.size() && path_lengths[iv].path_length < next_path_length) {
         int index = path_lengths[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index])
            ncorrect++;
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraPathLengthPDF(int min_path_length, int max_path_length, std::vector<dijkstra_path_length_metric>& path_lengths)
{
   printf("Starting dijkstra path length pdf..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_dijkstra_path_length_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < path_lengths.size() && path_lengths[iv].path_length < next_path_length) {
         int index = path_lengths[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
DijkstraPathLengthPDFBoundary(int min_path_length, int max_path_length, std::vector<dijkstra_path_length_metric>& path_lengths)
{
   printf("Starting dijkstra path length pdf for boundaries..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_dijkstra_path_length_boundary_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < path_lengths.size() && path_lengths[iv].path_length < next_path_length) {
         int index = path_lengths[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index])
            ncorrect++;
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelDistanceCDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel distance cdf..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_distance_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelDistanceCDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel distance cdf for boundaries..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_distance_boundary_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // increment number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelDistancePDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel distance pdf..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_distance_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelDistancePDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel distance pdf for boundaries..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_distance_boundary_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // incerement number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelScaledDistanceCDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel scaled distance cdf..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_scaled_distance_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelScaledDistanceCDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel scaled distance cdf for boundaries..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_scaled_distance_boundary_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // increment number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelScaledDistancePDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel scaled distance pdf..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_scaled_distance_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelScaledDistancePDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel scaled distance pdf for boundaries..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_scaled_distance_boundary_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // incerement number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelNormalizedDistanceCDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel normalized distance cdf..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_normalized_distance_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelNormalizedDistanceCDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel normalized distance cdf for boundaries..."); fflush(stdout);

   // get path length cdf filename
   char path_length_cdf_filename[4096];
   sprintf(path_length_cdf_filename, "%s_voxel_normalized_distance_boundary_cdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_cdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_cdf_filename); return 0; }

   // get path length range and bin size
   RNScalar path_range = (max_path_length - min_path_length);
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;
   // make sure bin is not too small
   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned long long ncorrect = 0;
   unsigned long long nnodes = 0;
   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // increment number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelNormalizedDistancePDF(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel normalized distance pdf..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_normalized_distance_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         ncorrect += voxel_ncorrect[index];
         nnodes += voxel_path_length[index];

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
VoxelNormalizedDistancePDFBoundary(int min_path_length, int max_path_length, std::vector<voxel_distance_metric>& distances)
{
   printf("Starting voxel normalized distance pdf for boundaries..."); fflush(stdout);

   // get path length pdf filename
   char path_length_pdf_filename[4096];
   sprintf(path_length_pdf_filename, "%s_voxel_normalized_distance_boundary_pdf.csv", root_output_filename);

   // open file
   FILE *fp = fopen(path_length_pdf_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", path_length_pdf_filename); return 0; }

   // get path length range and bin sizes
   RNScalar path_range = max_path_length - min_path_length;
   // add a small factor for double precision errors
   RNScalar path_range_bin_size = path_range / nbins + 10e-6;

   RNScalar next_path_length = path_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   unsigned int iv = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      unsigned long long ncorrect = 0;
      unsigned long long nnodes = 0;

      while (iv < distances.size() && distances[iv].voxel_distance < next_path_length) {
         int index = distances[iv].index;

         // update previous and next distances
         if (voxel_last_correct[index]) ncorrect++;

         // incerement number of nodes
         nnodes++;

         // increment iv
         iv++;
      }
      fprintf(fp, "%f, %f\n", prev_path_length, ((RNScalar)ncorrect) / nnodes);

      prev_path_length = next_path_length;
      next_path_length = next_path_length + path_range_bin_size;
   }

   // close file
   fclose(fp);

   printf("done.\n");

   // return success
   return 1;
}



static int
OutputDistributionFiles(void)
{
   // create vectors
   std::vector<dijkstra_distance_metric> distances = std::vector<dijkstra_distance_metric>();
   std::vector<dijkstra_path_length_metric> path_lengths = std::vector<dijkstra_path_length_metric>();
   std::vector<voxel_distance_metric> voxel_distances = std::vector<voxel_distance_metric>();
   std::vector<voxel_distance_metric> voxel_scaled_distances = std::vector<voxel_distance_metric>();
   std::vector<voxel_distance_metric> voxel_normalized_distances = std::vector<voxel_distance_metric>();

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      // skip extracellulars
      if (nd->Voxel(iv)->HumanLabel() == NULL) continue;

      // create metric objects
      dijkstra_distance_metric distance = dijkstra_distance_metric();
      dijkstra_path_length_metric path_length = dijkstra_path_length_metric();

      // update attributes
      distance.distance = dijkstra_distances[iv];
      distance.index = iv;
      path_length.path_length = voxel_path_length[iv];
      path_length.index = iv;

      // add to vectors
      distances.push_back(distance);
      path_lengths.push_back(path_length);


      // get the other distances
      int ix, iy, iz;
      nd->IndexToIndices(iv, ix, iy, iz);

      // find x, y, and z closest distances
      int xdistance = ix;
      int ydistance = iy;
      int zdistance = iz;
      if (ix >= nd->XResolution() / 2) xdistance = nd->XResolution() - 1 - ix;
      if (iy >= nd->YResolution() / 2) ydistance = nd->YResolution() - 1 - iy;
      if (iz >= nd->ZResolution() / 2) zdistance = nd->ZResolution() - 1 - iz;

      // find x, y, and z scale distances
      int xscaled_distance = scaling[RN_X] * xdistance;
      int yscaled_distance = scaling[RN_Y] * ydistance;
      int zscaled_distance = scaling[RN_Z] * zdistance;

      // find x, y, and z normalized distances
      int xnormalized_distance = xdistance * nd->XScale();
      int ynormalized_distance = ydistance * nd->YScale();
      int znormalized_distance = zdistance * nd->ZScale();

      int closest_distance = min(min(xdistance, ydistance), zdistance);
      int closest_scaled_distance = min(min(xscaled_distance, yscaled_distance), zscaled_distance);
      int closest_normalized_distance = min(min(xnormalized_distance, ynormalized_distance), znormalized_distance);

      voxel_distance_metric voxel_distance = voxel_distance_metric();
      voxel_distance.voxel_distance = closest_distance;
      voxel_distance.index = iv;

      voxel_distance_metric voxel_scaled_distance = voxel_distance_metric();
      voxel_scaled_distance.voxel_distance = closest_scaled_distance;
      voxel_scaled_distance.index = iv;

      voxel_distance_metric voxel_normalized_distance = voxel_distance_metric();
      voxel_normalized_distance.voxel_distance = closest_normalized_distance;
      voxel_normalized_distance.index = iv;

      // add to vectors
      voxel_distances.push_back(voxel_distance);
      voxel_scaled_distances.push_back(voxel_scaled_distance);
      voxel_normalized_distances.push_back(voxel_normalized_distance);
   }

   // sort the function
   std::sort(distances.begin(), distances.end(), distance_compare);
   std::sort(path_lengths.begin(), path_lengths.end(), path_length_compare);
   std::sort(voxel_distances.begin(), voxel_distances.end(), voxel_distance_compare);
   std::sort(voxel_scaled_distances.begin(), voxel_scaled_distances.end(), voxel_distance_compare);
   std::sort(voxel_normalized_distances.begin(), voxel_normalized_distances.end(), voxel_distance_compare);

   // get minimum and maximum distance
   RNScalar min_distance = distances[0].distance;
   RNScalar max_distance = distances[distances.size() - 1].distance;

   // get minimum and maximum path length
   int min_path_length = path_lengths[0].path_length;
   int max_path_length = path_lengths[path_lengths.size() - 1].path_length;

   // get minimum and maximum voxel distance
   int min_voxel_distance = voxel_distances[0].voxel_distance;
   int max_voxel_distance = voxel_distances[voxel_distances.size() - 1].voxel_distance;

   // get minimum and maximum voxel scaled distance
   int min_voxel_scaled_distance = voxel_scaled_distances[0].voxel_distance;
   int max_voxel_scaled_distance = voxel_scaled_distances[voxel_scaled_distances.size() - 1].voxel_distance;

   // get minimum and maximum voxel normalized distance
   int min_voxel_normalized_distance = voxel_normalized_distances[0].voxel_distance;
   int max_voxel_normalized_distance = voxel_normalized_distances[voxel_normalized_distances.size() - 1].voxel_distance;


   ////////////////////////////
   //// DIJKSTRA DISTANCES ////
   ////////////////////////////

   if (cdf) {
      if (!DijkstraDistanceCDF(min_distance, max_distance, distances)) return 0;
      if (!DijkstraDistanceCDFBoundary(min_distance, max_distance, distances)) return 0;
   }
   if (pdf) {
      if (!DijkstraDistancePDF(min_distance, max_distance, distances)) return 0;
      if (!DijkstraDistancePDFBoundary(min_distance, max_distance, distances)) return 0;
   }


   ///////////////////////////////
   //// DIJKSTRA PATH LENGTHS ////
   ///////////////////////////////

   if (cdf) {
      if (!DijkstraPathLengthCDF(min_path_length, max_path_length, path_lengths)) return 0;
      if (!DijkstraPathLengthCDFBoundary(min_path_length, max_path_length, path_lengths)) return 0;
   }
   if (pdf) {
      if (!DijkstraPathLengthPDF(min_path_length, max_path_length, path_lengths)) return 0;
      if (!DijkstraPathLengthPDFBoundary(min_path_length, max_path_length, path_lengths)) return 0;
   }


   /////////////////////////
   //// VOXEL DISTANCES ////
   /////////////////////////

   if (cdf) {
      if (!VoxelDistanceCDF(min_voxel_distance, max_voxel_distance, voxel_distances)) return 0;
      if (!VoxelDistanceCDFBoundary(min_voxel_distance, max_voxel_distance, voxel_distances)) return 0;
   }
   if (pdf) {
      if (!VoxelDistancePDF(min_voxel_distance, max_voxel_distance, voxel_distances)) return 0;
      if (!VoxelDistancePDFBoundary(min_voxel_distance, max_voxel_distance, voxel_distances)) return 0;
   }


   ////////////////////////////////
   //// VOXEL SCALED DISTANCES ////
   ////////////////////////////////

   if (cdf) {
      if (!VoxelScaledDistanceCDF(min_voxel_scaled_distance, max_voxel_scaled_distance, voxel_scaled_distances)) return 0;
      if (!VoxelScaledDistanceCDFBoundary(min_voxel_scaled_distance, max_voxel_scaled_distance, voxel_scaled_distances)) return 0;
   }
   if (pdf) {
      if (!VoxelScaledDistancePDF(min_voxel_scaled_distance, max_voxel_scaled_distance, voxel_scaled_distances)) return 0;
      if (!VoxelScaledDistancePDFBoundary(min_voxel_scaled_distance, max_voxel_scaled_distance, voxel_scaled_distances)) return 0;
   }


   ////////////////////////////////////
   //// VOXEL NORMALIZED DISTANCES ////
   ////////////////////////////////////

   if (cdf) {
      if (!VoxelNormalizedDistanceCDF(min_voxel_normalized_distance, max_voxel_normalized_distance, voxel_normalized_distances)) return 0;
      if (!VoxelNormalizedDistanceCDFBoundary(min_voxel_normalized_distance, max_voxel_normalized_distance, voxel_normalized_distances)) return 0;
   }
   if (pdf) {
      if (!VoxelNormalizedDistancePDF(min_voxel_normalized_distance, max_voxel_normalized_distance, voxel_normalized_distances)) return 0;
      if (!VoxelNormalizedDistancePDFBoundary(min_voxel_normalized_distance, max_voxel_normalized_distance, voxel_normalized_distances)) return 0;
   }

   // return success
   return 1;
}



///////////////////////////////////////////////////////////////////////
// Trace Output functions
///////////////////////////////////////////////////////////////////////

static int 
DijkstraDistanceTrace(void) {
   // start statistics
   RNTime start_time;
   start_time.Read();

   printf("Starting dijkstra distance trace metric...\n"); fflush(stdout);

   // get distance cdf filename
   char distance_trace_filename[4096];
   sprintf(distance_trace_filename, "%s_dijkstra_distance_trace.csv", root_output_filename);

   // get the maximum and minimum distance
   RNScalar min_distance = FLT_MAX;
   RNScalar max_distance = FLT_MIN;
   int min_path_length = INT_MAX;
   int max_path_length = INT_MIN;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel() == NULL) continue;

      if (dijkstra_distances[iv] < min_distance) {
         min_distance = dijkstra_distances[iv];
      }
      if (dijkstra_distances[iv] > max_distance) {
         max_distance = dijkstra_distances[iv];
      }
      if (voxel_path_length[iv] < min_path_length) {
         min_path_length = voxel_path_length[iv];
      }
      if (voxel_path_length[iv] > max_path_length) {
         max_path_length = voxel_path_length[iv];
      }
   }

   // set the bin size
   nbins = 500;

   // get distance range and bin size
   RNScalar distance_range = max_distance - min_distance;
   RNScalar distance_range_bin_size = distance_range / nbins + 10e-6;

   // get path length range and bin size
   RNScalar path_length_range = max_path_length - min_path_length;
   RNScalar path_length_range_bin_size = path_length_range / nbins + 10e-6;

   // open file
   FILE *fp = fopen(distance_trace_filename, "w");
   if (!fp) { fprintf(stderr, "Unable to write to %s\n", distance_trace_filename); return 0; }

   // create bins for distance
   unsigned long long *distance_ncorrect = new unsigned long long[nbins];
   unsigned long long *distance_nnodes = new unsigned long long[nbins];
   // create bins for path length
   unsigned long long *path_length_ncorrect = new unsigned long long[nbins];
   unsigned long long *path_length_nnodes = new unsigned long long[nbins];

   for (int ib = 0; ib < nbins; ++ib) {
      distance_ncorrect[ib] = 0;
      distance_nnodes[ib] = 0;
      path_length_ncorrect[ib] = 0;
      path_length_nnodes[ib] = 0;
   }

   int nvoxels = 0;
   unsigned long long nsupervoxel_changes = 0;
   // traverse through every dijkstra path
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (iv % 100000 == 0) printf("   %.2f%% done; %.2f seconds remaining\n", 100 * ((RNScalar)iv) / nd->NVoxels(), (nd->NVoxels() - iv) / ((RNScalar)iv) * start_time.Elapsed());
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel() == NULL) continue;
      nvoxels++;
      // get the entire path with distances
      std::vector<int> voxel_path = std::vector<int>();
      std::vector<RNScalar> distances = std::vector<RNScalar>();
      int previous = iv;
      int index;
      RNScalar cumulative_distance = 0.0;
      while (true) {
         index = previous;
         previous = dijkstra_prev_node[index];

         voxel_path.push_back(index);
         distances.push_back(cumulative_distance);

         if (index == previous) break;

         // get distance between previous and int
         NeuronVoxel *previous_voxel = nd->Voxel(previous);
         NeuronVoxel *index_voxel = nd->Voxel(index);

         if (cumulative_distance < 0.076643)
            if (previous_voxel->Supervoxel() != index_voxel->Supervoxel()) nsupervoxel_changes++;

         // get indices
         int ix, iy, iz;
         nd->IndexToIndices(index, ix, iy, iz);
         int ii, ij, ik;
         nd->IndexToIndices(previous, ii, ij, ik);

         RNScalar affinity = previous_voxel->AffinityToNeighbor(index_voxel);

         // get dijkstra distance
         RNScalar dijkstra_distance;

         if (ix - ii)
            dijkstra_distance = scaling[RN_X] * (max_affinity - affinity);
         else if (iy - ij)
            dijkstra_distance = scaling[RN_Y] * (max_affinity - affinity);
         else if (iz - ik)
            dijkstra_distance = scaling[RN_Z] * (max_affinity - affinity);
         else {
            dijkstra_distance = FLT_MAX;
            rn_assertion(FALSE);
         }

         // distance must be non-negative
         if (dijkstra_distance < 0.0) dijkstra_distance = 0.0;

         cumulative_distance += dijkstra_distance;
      }

      // fill in the bins
      for (unsigned int ip = 0; ip < voxel_path.size(); ++ip) {
         RNScalar dijkstra_distance = distances[ip];

         int distance_bin = (int)((dijkstra_distance - min_distance) / distance_range_bin_size);
         int path_length_bin = (int)((ip + 1 - min_path_length) / path_length_range_bin_size);

         rn_assertion((0 <= distance_bin) && (distance_bin < nbins));
         rn_assertion((0 <= path_length_bin) && (path_length_bin < nbins));
         // see if the truth index matches
         NeuronHumanLabel *human_label = nd->Voxel(iv)->HumanLabel();
         NeuronHumanLabel *path_human_label = nd->Voxel(voxel_path[ip])->HumanLabel();
         if (human_label == path_human_label) {
            distance_ncorrect[distance_bin]++;
            path_length_ncorrect[path_length_bin]++;
         }
         distance_nnodes[distance_bin]++;
         path_length_nnodes[path_length_bin]++;
      }
   }

   // create next and prev metrics
   RNScalar next_distance = distance_range_bin_size + min_distance;
   RNScalar prev_distance = min_distance;
   RNScalar next_path_length = path_length_range_bin_size + min_path_length;
   RNScalar prev_path_length = min_path_length;

   /* rn_assertion variables */
   unsigned long long distance_total_correct = 0;
   unsigned long long distance_total_nodes = 0;
   unsigned long long path_length_total_correct = 0;
   unsigned long long path_length_total_nodes = 0;
   for (int ib = 0; ib < nbins; ++ib) {
      fprintf(fp, "%f,%f,%f,%f\n", prev_distance, ((RNScalar)distance_ncorrect[ib]) / distance_nnodes[ib], prev_path_length, ((RNScalar)path_length_ncorrect[ib]) / path_length_nnodes[ib]);
      distance_total_correct += distance_ncorrect[ib];
      distance_total_nodes += distance_nnodes[ib];
      path_length_total_correct += path_length_ncorrect[ib];
      path_length_total_nodes += path_length_nnodes[ib];

      // update prev and next distances
      prev_distance = next_distance;
      next_distance = distance_range_bin_size + prev_distance;
      // update prev and next path lengths
      prev_path_length = next_path_length;
      next_path_length = path_length_range_bin_size + prev_path_length;
   }
   rn_assertion(distance_total_correct == path_length_total_correct);
   rn_assertion(distance_total_nodes == path_length_total_nodes);

   printf("%llu/%d = %f\n", nsupervoxel_changes, nvoxels, (RNScalar)nsupervoxel_changes / nvoxels);

   // close the file
   fclose(fp);

   // get total number of correct in both
   unsigned long long predicted_correct = 0;
   unsigned long long predicted_path_length = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (nd->Voxel(iv)->HumanLabel() == NULL) continue;
      predicted_correct += voxel_ncorrect[iv];
      predicted_path_length += voxel_path_length[iv];
   }

   // free memory
   delete[] distance_ncorrect;
   delete[] distance_nnodes;
   delete[] path_length_ncorrect;
   delete[] path_length_nnodes;


   printf("done.\n");

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
   // Parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) { print_debug = 1; print_verbose = 1; }
         else if (!strcmp(*argv, "-max_affinity")) { argv++; argc--; max_affinity = ((RNScalar)atoi(*argv)) / 100.0; }
         else if (!strcmp(*argv, "-all")) { cdf = 1; pdf = 1; trace = 1; }
         else if (!strcmp(*argv, "-cdf")) cdf = 1;
         else if (!strcmp(*argv, "-pdf")) pdf = 1;
         else if (!strcmp(*argv, "-trace")) trace = 1;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) { fprintf(stderr, "Need to supply input filename.\n"); return 0; }

   // make sure at least one metric selected
   if (!cdf && !pdf && !trace) { fprintf(stderr, "Must select at least one statistic metric (e.g. cdf, pdf, trace)\n");  return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{

   // Parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in the nueron data
   if (!ReadData()) exit(-1);

   // update scaling
   for (int dim = 0; dim <= 2; ++dim) {
      scaling[dim] = nd->Resolution(RN_X) / nd->Resolution(dim);
   }

   // create input filename root
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strchr(root_filename, '.');
   *extp = '\0';

   sprintf(root_input_filename, "results/dijkstra/%s_voxel_%03d", root_filename, (int) (100 * max_affinity + 0.5));
   sprintf(root_output_filename, "statistics/dijkstra/%s_voxel_%03d", root_filename, (int)(100 * max_affinity + 0.5));

   // read previous node file
   char input_byte_name[4096];
   sprintf(input_byte_name, "%s_prev_node.bin", root_input_filename);
   dijkstra_prev_node = new int[nd->NVoxels()];
   if (!dijkstra_prev_node) { fprintf(stderr, "Failed to allocate memory for previous nodes\n"); exit(-1); }
   FILE *fp = fopen(input_byte_name, "rb");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", input_byte_name); exit(-1); }
   fread(&(dijkstra_prev_node[0]), sizeof(int), nd->NVoxels(), fp);
   fclose(fp);

   // read distance file
   sprintf(input_byte_name, "%s_distances.bin", root_input_filename);
   dijkstra_distances = new RNScalar[nd->NVoxels()];
   if (!dijkstra_distances) { fprintf(stderr, "Failed to allocate memory for distances\n"); exit(-1); }
   fp = fopen(input_byte_name, "rb");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", input_byte_name); exit(-1); }
   fread(&(dijkstra_distances[0]), sizeof(RNScalar), nd->NVoxels(), fp);
   fclose(fp);


   // create array for path length
   voxel_path_length = new int[nd->NVoxels()];
   if (!voxel_path_length) { fprintf(stderr, "Failed to allocate memory for path length\n"); exit(-1); }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      voxel_path_length[iv] = 0;
   }

   // create array for voxel ncorrect
   voxel_ncorrect = new int[nd->NVoxels()];
   if (!voxel_ncorrect) { fprintf(stderr, "Failed to allocate memory for ncollect\n"); exit(-1); }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      voxel_ncorrect[iv] = 0;
   }

   // create array for last path variable
   voxel_last_correct = new RNBoolean[nd->NVoxels()];
   if (!voxel_last_correct) { fprintf(stderr, "Failed to allocate memory for last on path\n"); exit(-1); }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      voxel_last_correct[iv] = FALSE;
   }

   // update various path variables
   UpdatePathVariables();

   // only perform the following steps if cdf or pdf flag is present
   if (cdf || pdf) {
      // output csv files
      if (!OutputDistributionFiles()) exit(-1);
   }

   // create trace files if need
   if (trace) {
      if (!DijkstraDistanceTrace()) exit(-1);
   }

   // free memory
   delete[] voxel_last_correct;
   delete[] voxel_ncorrect;
   delete[] voxel_path_length;

   delete[] dijkstra_distances;
   delete[] dijkstra_prev_node;
   delete nd;

   // Return success
   return 0;
}
