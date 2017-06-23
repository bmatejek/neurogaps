// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int cellular_truth = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *truth_directory = "truth";



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
// Truth generation functions
////////////////////////////////////////////////////////////////////////

static int GenerateCellularTruth(char root_filename[4096])
{
   // write truth file
   char majority_truth_filename[4096];
   sprintf(majority_truth_filename, "%s/%s_cellular_majority.truth", truth_directory, root_filename);
   char center_truth_filename[4096];
   sprintf(center_truth_filename, "%s/%s_cellular_center.truth", truth_directory, root_filename);

   // open files
   FILE *majority_fp = fopen(majority_truth_filename, "wb");
   if (!majority_fp) { fprintf(stderr, "Failed to write to %s\n", majority_truth_filename); return 0; }
   FILE *center_fp = fopen(center_truth_filename, "wb");
   if (!center_fp) { fclose(majority_fp); fprintf(stderr, "Failed to write to %s\n", center_truth_filename); return 0; }

   // write truth between all supervoxels (do they belong to the same human label)
   int ncellulars = nd->NCellulars();
   fwrite(&ncellulars, sizeof(int), 1, majority_fp);
   fwrite(&ncellulars, sizeof(int), 1, center_fp);

   for (int i = 0; i < ncellulars; ++i) {
      NeuronCellular *cellular_one = nd->Cellular(i);
      // get the human labels
      int majority_human_label_one = cellular_one->MajorityHumanLabelIndex();

      for (int j = 0; j < ncellulars; ++j) {
         NeuronCellular *cellular_two = nd->Cellular(j);
         // get neighbor human labels
         int majority_human_label_two = cellular_two->MajorityHumanLabelIndex();

         RNBoolean majority_match = FALSE;
         RNBoolean center_match = FALSE;
         if (majority_human_label_one == majority_human_label_two) majority_match = TRUE;

         // write match results
         fwrite(&majority_match, sizeof(RNBoolean), 1, majority_fp);
         fwrite(&center_match, sizeof(RNBoolean), 1, center_fp);
      }
   }

   // close files
   fclose(majority_fp);
   fclose(center_fp);

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
         else if (!strcmp(*argv, "-cellular")) cellular_truth = 1;
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

   // truth creation functions
   if (cellular_truth && !GenerateCellularTruth(root_filename)) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}
