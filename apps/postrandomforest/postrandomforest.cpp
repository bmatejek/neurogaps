// Source file for the simple file conversion algorithm



// include files

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
// string arguments
static const char *extension = "post";
static char *random_forest_filename = NULL;
static char *testing_input_filename = NULL;



// global data sets and random forest

static NeuronData *testing_nd = NULL;
static RNDataset *testing_dataset = NULL;
static RNRandomForest *random_forest = NULL;



// directory structure
static const char *dataset_directory = "algs_data/dataset";
static const char *graphcut_directory = "algs_data/graphcut";



// input directories where random forests are stored
static const char *random_forest_directory = "algs_data/random_forest";



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
   if (!nd->ReadFile(filename, TRUE)) {
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



static int WriteGraphcutFiles(NeuronData *nd, RNDataset *dataset)
{
   // get the vector results
   std::vector<float> results = random_forest->PredictBinaryProbabilities(*dataset);

   // get root filenames
   char root_filename[4096];
   strcpy(root_filename, nd->Filename());
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // allocate memory for array
   float *smoothing_terms = new float[nd->NPredictionBoundaries()];
   for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb)
      // default is no smoothing penalty
      smoothing_terms[ipb] = 0.0;

   // go through all dataset entries
   for (unsigned int in = 0; in < dataset->NEntries(); ++in) {
      // get boundary identification
      int identification = dataset->GetIdentification(in);

      // set the results
      smoothing_terms[identification] = results[in];
   }

   // get output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_random_forest_%s.binary", graphcut_directory, root_filename, extension);

   FILE *binary_fp = fopen(output_filename, "wb");
   if (!binary_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // write boundary terms
   int nboundaries = nd->NPredictionBoundaries();
   fwrite(&nboundaries, sizeof(int), 1, binary_fp);
   fwrite(smoothing_terms, sizeof(float), nd->NPredictionBoundaries(), binary_fp);

   fclose(binary_fp);

   // free memory
   delete[] smoothing_terms;

   // return success
   return 1;
}



static int ReadDataset(void)
{
   // get data root filename
   char testing_root_filename[4096];
   strncpy(testing_root_filename, testing_nd->Filename(), 4096);
   char *testing_extp = strrchr(testing_root_filename, '.');
   *testing_extp = '\0';

   // get dataset filenames
   char testing_dataset_filename[4096];
   sprintf(testing_dataset_filename, "%s/%s_%s.mldb", dataset_directory, testing_root_filename, extension);

   // train random forest classifier
   char input_filename[4096];
   sprintf(input_filename, "%s/%s_%s.xml", random_forest_directory, random_forest_filename, extension);

   // load the random forest
   random_forest = new RNRandomForest();
   random_forest->LoadClassifier(input_filename);

   // read in the testing dataset
   testing_dataset = new RNDataset();
   testing_dataset->ReadFile(testing_dataset_filename);
   if (!WriteGraphcutFiles(testing_nd, testing_dataset)) exit(-1);

   // print statistics
   printf("Testing dataset results for %s terms:\n", extension);
   unsigned int **confusion_matrix = random_forest->PrintPerformance(*testing_dataset);

   // free memory
   for (unsigned int ic = 0; ic < testing_dataset->NClasses(); ++ic)
      delete[] confusion_matrix[ic];
   delete[] confusion_matrix;
   delete testing_dataset;
   delete random_forest;

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
         else if (!strcmp(*argv, "-forest")) { argv++; argc--; random_forest_filename = *argv; }
	 else if (!strcmp(*argv, "-all")) { extension = "post_all"; }
         else if (!strcmp(*argv, "-test")) { argv++; argc--; testing_input_filename = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0;
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!random_forest_filename) {
      fprintf(stderr, "Need to supply random forest filename.\n");
      return 0;
   }
   if (!testing_input_filename) {
      fprintf(stderr, "Need to supply testing input filename.\n");
      return 0;
   }

   // additional input verification
   char *extp = strrchr(random_forest_filename, '.');
   if (extp) { fprintf(stderr, "Root filename should not have an extension.\n"); return 0; }
   char *prep = strrchr(random_forest_filename, '/');
   if (prep) { fprintf(stderr, "Root filename should not include directory structure.\n"); return 0; }

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
   if (!testing_nd) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, testing_nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in prediction file
   char meta_root_filename[4096];
   sprintf(meta_root_filename, "output/ordermerge/%s", root_filename);

   R3Grid *prediction_grid = RNReadNeuronMetaRawFile(meta_root_filename);
   if (!prediction_grid) exit(-1);

   // create boundary predictions
   testing_nd->CreatePredictions(prediction_grid);

   if (!ReadDataset()) exit(-1);

   // free up memory
   delete testing_nd;

   // return success
   return 0;
}
