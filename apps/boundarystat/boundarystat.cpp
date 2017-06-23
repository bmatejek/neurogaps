// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static float **data_terms = NULL;



// directory structure

static const char *graphcut_directory = "algs_data/graphcut";



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
// Result functions
////////////////////////////////////////////////////////////////////////

static int InteriorStats(void)
{
   // counter variables
   int true_positives = 0;
   int true_negatives = 0;
   int false_positives = 0;
   int false_negatives = 0;

   // go through all pairs of boundary cellulars
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      NeuronCellular *cellular_one = nd->Cellular(ic1);
      if (!cellular_one->IsOnBoundary()) continue;
      for (int ic2 = ic1 + 1; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);
         if (!cellular_two->IsOnBoundary()) continue;

         int cellular_one_index = cellular_one->DataIndex();
         int cellular_two_index = cellular_two->DataIndex();

         int prediction = 1 - (int)(data_terms[cellular_one_index][cellular_two_index] + 0.5);
         int human_label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());

         if (prediction && human_label) true_positives++;
         else if (prediction && !human_label) false_positives++;
         else if (!prediction && human_label) false_negatives++;
         else if (!prediction && !human_label) true_negatives++;
         else rn_assertion(FALSE);

      }
   }

   // print results
   printf("Interior Statistics: \n");
   printf("TP: %8d\n", true_positives);
   printf("TN: %8d\n", true_negatives);
   printf("FP: %8d\n", false_positives);
   printf("FN: %8d\n", false_negatives);
   printf("------------\n");
   printf("Total Proportion Correct: %lf\n", (true_positives + true_negatives) / (RNScalar)(true_positives + true_negatives + false_positives + false_negatives));
   printf("Total Proportion of Positive Guesses: %lf\n", true_positives / (RNScalar)(true_positives + false_positives));
   printf("Total Proportion of Negative Guesses: %lf\n", true_negatives / (RNScalar)(true_negatives + false_negatives));

   // return success
   return 1;
}



static int BoundaryStats(void)
{
   // counter variables
   int true_positives = 0;
   int true_negatives = 0;
   int false_positives = 0;
   int false_negatives = 0;

   // go through all pairs of boundary supervoxels
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      NeuronCellular *cellular_one = nd->Cellular(ic1);
      if (!cellular_one->IsOnBoundary()) continue;
      for (int ic2 = ic1 + 1; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);
         if (!cellular_two->IsOnBoundary()) continue;

         int prediction = 1 - (int)(data_terms[ic1][ic2] + 0.5);
         int human_label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());

         if (prediction && human_label) true_positives++;
         else if (prediction && !human_label) false_positives++;
         else if (!prediction && human_label) false_negatives++;
         else if (!prediction && !human_label) true_negatives++;
         else rn_assertion(FALSE);
      }
   }

   // print results
   printf("Boundary Statistics: \n");
   printf("TP: %8d\n", true_positives);
   printf("TN: %8d\n", true_negatives);
   printf("FP: %8d\n", false_positives);
   printf("FN: %8d\n", false_negatives);
   printf("------------\n");
   printf("Total Proportion Correct: %lf\n", (true_positives + true_negatives) / (RNScalar)(true_positives + true_negatives + false_positives + false_negatives));
   printf("Total Proportion of Positive Guesses: %lf\n", true_positives / (RNScalar)(true_positives + false_positives));
   printf("Total Proportion of Negative Guesses: %lf\n", true_negatives / (RNScalar)(true_negatives + false_negatives));

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
   strncpy(root_filename, nd->Filename(), '.');
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // get data term filename
   char data_filename[4096];
   sprintf(data_filename, "%s/%s_random_forest.unary", graphcut_directory, root_filename);

   // open file
   FILE *fp = fopen(data_filename, "rb");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", data_filename); exit(-1); }

   // read in data terms
   data_terms = new float *[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic)
      data_terms[ic] = new float[nd->NCellulars()];

   // check
   int ncellulars;
   fread(&ncellulars, sizeof(int), 1, fp);
   rn_assertion(ncellulars == nd->NCellulars());

   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      fread(data_terms[ic], sizeof(float), nd->NCellulars(), fp);
   }

   // close file
   fclose(fp);

   // run statistic functions
   //if (!InteriorStats()) exit(-1);
   if (!BoundaryStats()) exit(-1);

   // free up memory
   delete nd;
   for (int ic = 0; ic < nd->NCellulars(); ++ic)
      delete[] data_terms[ic];
   delete[] data_terms;

   // return success
   return 0;
}
