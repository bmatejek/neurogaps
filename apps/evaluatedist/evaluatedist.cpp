// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *distances_directory = "distances";
static const char *results_directory = "results/distances";
static const char *truth_directory = "truth";

static const int NDISTANCE_METRICS = 3;
static const char *metrics[NDISTANCE_METRICS] = { "average_commute", "first_pass", "dijkstra" };



////////////////////////////////////////////////////////////////////////
// Useful structs
////////////////////////////////////////////////////////////////////////

struct distance_pair {
   distance_pair(RNBoolean same_neuron, RNScalar distance, int index_one, int index_two) {
      this->same_neuron = same_neuron;
      this->distance = distance;
      this->index_one = index_one;
      this->index_two = index_two;
   }

   RNBoolean same_neuron;
   RNScalar distance;
   int index_one;
   int index_two;
};



int distance_compare(distance_pair a, distance_pair b) {
   return a.distance < b.distance;
}



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
   if (!nd->ReadFile(filename)) {
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

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // get truth filename
   char truth_filename[4096];
   sprintf(truth_filename, "%s/%s_cellular_center.truth", truth_directory, root_filename);

   // open truth filename
   FILE *truth_fp = fopen(truth_filename, "rb");
   if (!truth_fp) { fprintf(stderr, "Failed to read %s\n", truth_filename); return 0; }

   int truth_ncellulars;
   fread(&truth_ncellulars, sizeof(int), 1, truth_fp);

   rn_assertion(truth_ncellulars == nd->NCellulars());

   RNBoolean **truth = new RNBoolean *[truth_ncellulars];
   for (int i = 0; i < truth_ncellulars; ++i) {
      truth[i] = new RNBoolean[truth_ncellulars];
      for (int j = 0; j < truth_ncellulars; ++j) {
         fread(&(truth[i][j]), sizeof(RNBoolean), 1, truth_fp);
      }
   }

   // close truth filename
   fclose(truth_fp);

   for (int id = 0; id < NDISTANCE_METRICS; ++id) {
      const char *distance_metric = metrics[id];

      printf("Calculating results for distance metric: %s\n", distance_metric);

      // get distance filename
      char distances_filename[4096];
      sprintf(distances_filename, "%s/%s_%s.distance", distances_directory, root_filename, distance_metric);

      // open distance filename
      FILE *distances_fp = fopen(distances_filename, "rb");
      if (!distances_fp) { fprintf(stderr, "Failed to read %s\n", distances_filename); return 0; }

      int distances_ncellulars;
      fread(&distances_ncellulars, sizeof(int), 1, distances_fp);
      rn_assertion(distances_ncellulars == truth_ncellulars);

      // get distance double array
      RNScalar **distances = new RNScalar *[distances_ncellulars];
      for (int ic = 0; ic < distances_ncellulars; ++ic) {
         distances[ic] = new RNScalar[distances_ncellulars];
         fread(distances[ic], sizeof(RNScalar), nd->NCellulars(), distances_fp);
      }

      // what percent of closest boundary distances are correct
      std::vector<std::vector<distance_pair> > all_distances = std::vector<std::vector<distance_pair> >();
      std::vector<std::vector<distance_pair> > boundary_distances = std::vector<std::vector<distance_pair> >();

      int first_boundary_correct_counter = 0;
      int first_distance_correct_counter = 0;

      for (int i = 0; i < distances_ncellulars; ++i) {
         all_distances.push_back(std::vector<distance_pair>());
         boundary_distances.push_back(std::vector<distance_pair>());
         for (int j = 0; j < distances_ncellulars; ++j) {
            // make sure distance is reasonable (very much hard coded)
            rn_assertion((distances[i][j] > -10e-6) && (distances[i][j] < 10e12));

            NeuronCellular *cellular_two = nd->Cellular(j);
            distance_pair pair = distance_pair(truth[i][j], distances[i][j], i, j);
            if (i != j) all_distances[i].push_back(pair);
            if (cellular_two->IsOnBoundary()) boundary_distances[i].push_back(pair);
         }

         // sort the vectors in increasing distance
         sort(all_distances[i].begin(), all_distances[i].end(), &distance_compare);
         sort(boundary_distances[i].begin(), boundary_distances[i].end(), &distance_compare);

         if (boundary_distances[i][0].same_neuron) first_boundary_correct_counter++;
         if (all_distances[i][0].same_neuron) first_distance_correct_counter++;
      }

      printf("   Proportion of boundary first correct: %f\n", (RNScalar)first_boundary_correct_counter / boundary_distances.size());
      printf("   Proportion of all first correct: %f\n", (RNScalar)first_distance_correct_counter / all_distances.size());

      // close distance filename
      fclose(distances_fp);

      // free memory
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         delete[] distances[ic];
      }
      delete[] distances;
   }

   // free memory
   for (int ic = 0; ic < truth_ncellulars; ++ic) {
      delete[] truth[ic];
   }
   delete[] truth;

   // free memory
   delete nd;

   // return success
   return 0;
}
