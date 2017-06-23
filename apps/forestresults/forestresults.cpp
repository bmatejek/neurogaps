// Source file for the simple file conversion algorithm



// include files

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int smoothing_terms = 0;
static int data_terms = 0;
static int boundary_terms = 0;
// feature options
static int boundary_features = 0;
static int shape_features = 0;
static int global_features = 0;
// string arguments
static const char *extension = NULL;
static char *testing_input_filename = NULL;



// global data sets and random forest

static NeuronData *testing_nd = NULL;



// directory structure

static const char *graphcut_directory = "algs_data/graphcut";



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static NeuronData *ReadData(const char *filename)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   NeuronData *nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate memory for neuron data\n");
      return NULL;
   }

   // read in the file
   if (!nd->ReadFile(filename)) {
      fprintf(stderr, "Failed to read file\n");
      return NULL;
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

   // return the neuron structure
   return nd;
}



static int GetResults(NeuronData *nd)
{
   // get root filenames
   char root_filename[4096];
   strcpy(root_filename, nd->Filename());
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // number of labeling classes
   const int nclasses = 2;

   if (smoothing_terms) {
      // allocate memory for array
      float *smoothing_terms = new float[nd->NBoundaries()];

      // get output filename
      char input_filename[4096];
      sprintf(input_filename, "%s/%s_random_forest_%s.binary", graphcut_directory, root_filename, extension);

      FILE *binary_fp = fopen(input_filename, "rb");
      if (!binary_fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

      // write boundary terms
      int nboundaries = nd->NBoundaries();
      fread(&nboundaries, sizeof(int), 1, binary_fp);
      fread(smoothing_terms, sizeof(float), nd->NBoundaries(), binary_fp);

      fclose(binary_fp);

      // print out the results
      unsigned int **confusion_matrix = new unsigned int *[nclasses];
      for (unsigned int i = 0; i < nclasses; ++i) {
         confusion_matrix[i] = new unsigned int[nclasses];
         for (unsigned int j = 0; j < nclasses; ++j) {
            confusion_matrix[i][j] = 0;
         }
      }

      int nentries = 0;
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);
         if (boundary->SupervoxelOne()->IsExtracellular()) continue;
         if (boundary->SupervoxelTwo()->IsExtracellular()) continue;

         float result = smoothing_terms[ib];
         int categorical_result = (int)(0.5 + result);
         int label = (boundary->SupervoxelOne()->MajorityHumanLabel() == boundary->SupervoxelTwo()->MajorityHumanLabel());
         confusion_matrix[label][categorical_result]++;
         nentries++;
      }
      printf("------------------------------------------------------\n");
      printf("Confusion Matrix:\n");
      printf("                  Predicted Class\n");
      printf("Actual Class");
      for (unsigned int i = 0; i < nclasses; ++i) {
         printf("%10d", i);
      }
      printf("\n");
      for (unsigned int i = 0; i < nclasses; ++i) {
         printf("%12d", i);
         for (unsigned int j = 0; j < nclasses; ++j) {
            printf("%10u", confusion_matrix[i][j]);
         }
         printf("\n");
      }

      unsigned int total_correct = 0;
      unsigned int *total_per_class = new unsigned int[nclasses];
      for (unsigned int i = 0; i < nclasses; ++i) {
         total_correct += confusion_matrix[i][i];
         total_per_class[i] = 0;
         for (unsigned int j = 0; j < nclasses; ++j) {
            total_per_class[i] += confusion_matrix[i][j];
         }
      }
      printf("------------------------------------------------------\n");
      printf("Correctly Labeled: %lf%%\n", 100 * (RNScalar)total_correct / nentries);
      for (unsigned int i = 0; i < nclasses; ++i)
         printf("Class %d Correctly Labeled: %lf%%\n", i, 100 * (RNScalar)confusion_matrix[i][i] / total_per_class[i]);
      printf("------------------------------------------------------\n");

      delete[] total_per_class;
      for (int i = 0; i < nclasses; ++i)
         delete[] confusion_matrix[i];
      delete[] confusion_matrix;

      // free memory
      delete[] smoothing_terms;
   }
   else if (data_terms) {
      // allocate data term array
      int ncellulars = nd->NCellulars();
      float **data_terms = new float *[ncellulars];
      for (int ic1 = 0; ic1 < ncellulars; ++ic1) {
         data_terms[ic1] = new float[ncellulars];
      }

      // get output filenames
      char input_filename[4096];
      sprintf(input_filename, "%s/%s_random_forest.unary", graphcut_directory, root_filename);

      // open files
      FILE *data_fp = fopen(input_filename, "rb");
      if (!data_fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

      fread(&ncellulars, sizeof(int), 1, data_fp);
      for (int ic = 0; ic < ncellulars; ++ic) {
         fread(data_terms[ic], sizeof(float), ncellulars, data_fp);
      }

      // close files
      fclose(data_fp);


      // print out the results
      unsigned int **confusion_matrix = new unsigned int *[nclasses];
      for (unsigned int i = 0; i < nclasses; ++i) {
         confusion_matrix[i] = new unsigned int[nclasses];
         for (unsigned int j = 0; j < nclasses; ++j) {
            confusion_matrix[i][j] = 0;
         }
      }

      int nentries = 0;
      for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
         NeuronCellular *cellular_one = nd->Cellular(ic1);
         if (!cellular_one->IsOnBoundary()) continue;
         for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
            NeuronCellular *cellular_two = nd->Cellular(ic2);

            float result = data_terms[ic1][ic2];
            int categorical_result = 1 - (int)(0.5 + result);
            int label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());
            confusion_matrix[label][categorical_result]++;
            nentries++;
         }
      }
      printf("------------------------------------------------------\n");
      printf("Confusion Matrix:\n");
      printf("                  Predicted Class\n");
      printf("Actual Class");
      for (unsigned int i = 0; i < nclasses; ++i) {
         printf("%10d", i);
      }
      printf("\n");
      for (unsigned int i = 0; i < nclasses; ++i) {
         printf("%12d", i);
         for (unsigned int j = 0; j < nclasses; ++j) {
            printf("%10u", confusion_matrix[i][j]);
         }
         printf("\n");
      }

      unsigned int total_correct = 0;
      unsigned int *total_per_class = new unsigned int[nclasses];
      for (unsigned int i = 0; i < nclasses; ++i) {
         total_correct += confusion_matrix[i][i];
         total_per_class[i] = 0;
         for (unsigned int j = 0; j < nclasses; ++j) {
            total_per_class[i] += confusion_matrix[i][j];
         }
      }
      printf("------------------------------------------------------\n");
      printf("Correctly Labeled: %lf%%\n", 100 * (RNScalar)total_correct / nentries);
      for (unsigned int i = 0; i < nclasses; ++i)
         printf("Class %d Correctly Labeled: %lf%%\n", i, 100 * (RNScalar)confusion_matrix[i][i] / total_per_class[i]);
      printf("------------------------------------------------------\n");

      delete[] total_per_class;
      for (int i = 0; i < nclasses; ++i)
         delete[] confusion_matrix[i];
      delete[] confusion_matrix;


      // free memory
      for (int ic = 0; ic < ncellulars; ++ic)
         delete[] data_terms[ic];
      delete[] data_terms;
   }
   else { fprintf(stderr, "Neither data or smoothing term selected.\n"); return 0; }

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
         else if (!strcmp(*argv, "-smoothing_term")) { smoothing_terms = 1; }
         else if (!strcmp(*argv, "-data_term")) { data_terms = 1; }
         else if (!strcmp(*argv, "-boundary_term")) { boundary_terms = 1; }
         // feature options
         else if (!strcmp(*argv, "-boundary_features")) { boundary_features = 1; }
         else if (!strcmp(*argv, "-shape_features")) { shape_features = 1; }
         else if (!strcmp(*argv, "-global_features")) { global_features = 1; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!testing_input_filename) testing_input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!testing_input_filename) {
      fprintf(stderr, "Need to supply testing input filename.\n");
      return 0;
   }

   // make sure there are consistent labeling options
   if (data_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      extension = "data";
   }
   else if (boundary_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      data_terms = 1;
      extension = "data_boundary";
   }
   else if (smoothing_terms) {
      if (global_features) extension = "smoothing_global";
      else extension = "smoothing";

      // turn on all feature options besides global
      boundary_features = 1;
      shape_features = 1;
   }
   else if (boundary_features || shape_features || global_features) {
      // make sure only one feature is selected
      if (boundary_features + shape_features + global_features > 1) {
         fprintf(stderr, "Can only select one among boundary, shape, distance, and global features.\n");
         return 0;
      }

      // get the correct extension
      if (boundary_features) extension = "boundary";
      if (shape_features) extension = "shape";
      if (global_features) extension = "global";

      // turn on smoothing terms to go to the correct functions
      smoothing_terms = 1;
   }
   else {
      fprintf(stderr, "Must select either data or smoothing terms.\n");
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

   // read in neuron files
   testing_nd = ReadData(testing_input_filename);

   if (!GetResults(testing_nd)) exit(-1);

   // free up memory
   delete testing_nd;

   // return success
   return 0;
}
