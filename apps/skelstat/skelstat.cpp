// Source file for getting statistics about human label/prediction skeletons



// include files 

#include "Neuron/Neuron.h"
#include <limits.h>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;



// global variables

static R3Grid *predictions = NULL;
static std::vector<int> prediction_points = std::vector<int>();
static R3Grid *human_labels = NULL;
static std::vector<int> human_label_points = std::vector<int>();
static NeuronData *nd = NULL;



// directory structure

static const char *skeleton_directory = "neuron_data/skeletons";



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(const char *filename)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate memory for neuron data\n");
      return 0;
   }

   // read in the file
   if (!nd->ReadFile(filename, TRUE)) {
      fprintf(stderr, "Failed to read file\n");
      return 0;
   }

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }
   if (print_debug) {
      printf("  Bounding Box: (%0.2f, %0.2f, %0.2f) to (%0.2f, %0.2f, %0.2f)\n", nd->GridBox().XMin(), nd->GridBox().YMin(), nd->GridBox().ZMin(), nd->GridBox().XMax(), nd->GridBox().YMax(), nd->GridBox().ZMax());
      printf("  Voxels: %d\n", nd->NVoxels());
      printf("  Cellulars: %d\n", nd->NCellulars());
      printf("  Extracellulars: %d\n", nd->NExtracellulars());
      printf("  Boundaries: %d\n", nd->NBoundaries());
      printf("  Human Labels: %d\n", nd->NHumanLabels());
      printf("  Predictions: %d\n", nd->NPredictions());
   }

   // return success
   return 1;
}




////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) {
      fprintf(stderr, "Need to supply input filename.\n");
      return 0;
   }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read neuron file
   if (!ReadData(input_filename)) exit(-1);

   // get input root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   int vessel_number = root_filename[strlen(root_filename) - 1] - '0';

   char prediction_filename[4096];
   sprintf(prediction_filename, "%s/skeleton_predictions%d", skeleton_directory, vessel_number);
   char human_label_filename[4096];
   sprintf(human_label_filename, "%s/skeleton_human_labels%d", skeleton_directory, vessel_number);

   predictions = RNReadNeuronMetaRawFile(prediction_filename);
   if (!predictions) { fprintf(stderr, "Failed to read %s\n", prediction_filename); exit(-1); }
   human_labels = RNReadNeuronMetaRawFile(human_label_filename);
   if (!human_labels) { fprintf(stderr, "Failed to read %s\n", human_label_filename); exit(-1); }

   // get list of points on the skeleton
   for (int iv = 0; iv < predictions->NEntries(); ++iv) {
      int prediction = (int)(predictions->GridValue(iv) + 0.5);
      int human_label = (int)(human_labels->GridValue(iv) + 0.5);

      // do these points belong to a skeleton
      if (prediction != 0) prediction_points.push_back(iv);
      if (human_label != 0) human_label_points.push_back(iv);
   }

   char distance_filename[4096];
   sprintf(distance_filename, "%s_skeleton_stats.csv", root_filename);

   FILE *distance_fp = fopen(distance_filename, "w");
   if (!distance_fp) { fprintf(stderr, "Failed to write to %s\n", distance_filename); exit(-1); }
   RNScalar sum_distance = 0.0;
   // for every point in the prediction, find the closest point on the human label
   for (unsigned int ip = 0; ip < prediction_points.size(); ++ip) {
      int prediction = prediction_points[ip];
      int ix, iy, iz;
      predictions->IndexToIndices(prediction, ix, iy, iz);
      unsigned int min_distance = UINT_MAX;
      for (unsigned int ih = 0; ih < human_label_points.size(); ++ih) {
         int human_label = human_label_points[ih];
         int ii, ij, ik;
         human_labels->IndexToIndices(human_label, ii, ij, ik);

         unsigned int distance_squared = (ix - ii) * (ix - ii) + (iy - ij) * (iy - ij) + (iz - ik) * (iz - ik);
         if (distance_squared < min_distance) {
            min_distance = distance_squared;
         }
      }

      // distance in microns
      RNScalar distance = sqrt(min_distance);
      sum_distance += distance;
      fprintf(distance_fp, "%lf\n", distance);
   }

   // close file
   fclose(distance_fp);

   // find the percent of my predictions that are in the ground truth
   int ncorrect_predictions = 0;
   for (unsigned int ip = 0; ip < prediction_points.size(); ++ip) {
      int prediction = prediction_points[ip];
      NeuronHumanLabel *human_label = nd->Voxel(prediction)->HumanLabel();
      if (human_label) ncorrect_predictions++;
   }

   printf("Number of voxels in prediction skeleton: %lu\n", prediction_points.size());
   printf("Number of voxels in human label skeleton: %lu\n", human_label_points.size());

   printf("Percent of predictions in ground truth: %lf%%\n", 100 * ncorrect_predictions / (RNScalar)prediction_points.size());

   int ncorrect_human_labels = 0;
   for (unsigned int ih = 0; ih < human_label_points.size(); ++ih) {
      int human_label = human_label_points[ih];
      NeuronSupervoxel *supervoxel = nd->Voxel(human_label)->Supervoxel();
      if (supervoxel->IsCellular()) ncorrect_human_labels++;
   }

   printf("Percent of ground truth in predictions: %lf%%\n", 100 * ncorrect_human_labels / (RNScalar)human_label_points.size());
   printf("Average distance: %lf\n", sum_distance / prediction_points.size());
   // free memory
   delete predictions;
   delete human_labels;

   // return success
   return 0;
}
