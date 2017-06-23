// Source file for the simple file conversion algorithm



// include files

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
// string arguments
static const char *extension = "post";
static char *output_root_filename = NULL;



// random forest parameters

static int max_depth = 7;
static int min_sample_count = 1;
static RNBoolean rebalance = FALSE;
static int nactive_vars = 0;
static int ntrees = 1000;
static float forest_accuracy = 0.01;



// global variables

static std::vector<const char *> training_filenames = std::vector<const char *>();



// global data sets and random forest

static std::vector<NeuronData *> training_nds = std::vector<NeuronData *>();
static RNRandomForest *random_forest = NULL;



// directory structure

static const char *dataset_directory = "algs_data/dataset";
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



////////////////////////////////////////////////////////////////////////
// Random forest creation functions
////////////////////////////////////////////////////////////////////////

static int ReadDataset(void)
{
   RNDataset *training_dataset = new RNDataset();

   for (unsigned int in = 0; in < training_nds.size(); ++in) {
      NeuronData *training_nd = training_nds[in];

      // get data root filenames
      char training_root_filename[4096];
      strncpy(training_root_filename, training_nd->Filename(), 4096);
      char *training_extp = strrchr(training_root_filename, '.');
      *training_extp = '\0';

      // get dataset filenames
      char training_dataset_filename[4096];
      sprintf(training_dataset_filename, "%s/%s_%s.mldb", dataset_directory, training_root_filename, extension);

      // read in the training dataset
      RNDataset *this_training_dataset = new RNDataset();
      this_training_dataset->ReadFile(training_dataset_filename);

      // set the names for this dataset
      if (in == 0) {
         std::vector<std::string> names = this_training_dataset->GetNames();
         training_dataset->SetNames(names);
      }

      // copy this training dataset to training dataset
      for (unsigned int id = 0; id < this_training_dataset->NEntries(); ++id) {
         const std::vector<float> &attributes = this_training_dataset->GetFeature(id);
         unsigned int label = this_training_dataset->GetLabel(id);
         unsigned int identification = this_training_dataset->GetIdentification(id);

         training_dataset->InsertDatapoint(attributes, label, identification);
      }
      // known memory leak (OK for now)
   }

   // train random forest classifier
   random_forest = new RNRandomForest(*training_dataset, max_depth, min_sample_count, rebalance, nactive_vars, ntrees, forest_accuracy);

   // save the random forest to an xml file
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_%s.xml", random_forest_directory, output_root_filename, extension);

   printf("Saving random forest classifier to %s...", output_filename); fflush(stdout);
   RNTime save_time;
   save_time.Read();
   random_forest->SaveClassifier(output_filename);
   printf("done in %0.2f seconds\n", save_time.Elapsed());

   // free memory
   delete training_dataset;
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
         else if (!strcmp(*argv, "-output_root_filename")) { argv++; argc--; output_root_filename = *argv; }
	 else if (!strcmp(*argv, "-all")) { extension = "post_all"; }
         // random forest parameters
         else if (!strcmp(*argv, "-max_depth")) { argv++; argc--; max_depth = atoi(*argv); }
         else if (!strcmp(*argv, "-min_sample_count")) { argv++; argc--; min_sample_count = atoi(*argv); }
         else if (!strcmp(*argv, "-rebalance")) rebalance = TRUE;
         else if (!strcmp(*argv, "-nactive_vars")) { argv++; argc--; nactive_vars = atoi(*argv); }
         else if (!strcmp(*argv, "-ntrees")) { argv++; argc--; ntrees = atoi(*argv); }
         else if (!strcmp(*argv, "-forest_accuracy")) { argv++; argc--; forest_accuracy = atof(*argv); }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         training_filenames.push_back(*argv);
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!training_filenames.size()) { fprintf(stderr, "Need to supply training input filename.\n"); return 0; }
   if (!output_root_filename) { fprintf(stderr, "Need to supply output root filename.\n"); return 0; }

   char *extp = strrchr(output_root_filename, '.');
   if (extp) { fprintf(stderr, "Root filename should not have an extension.\n"); return 0; }
   char *prep = strrchr(output_root_filename, '/');
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

   // write tree arguments
   if (print_verbose) {
      printf("Tree arguments:\n");
      printf("  Maximum Number of Trees: %d\n", ntrees);
      printf("  Number of Active Variables: %d\n", nactive_vars);
      printf("  Maximum Depth: %d\n", max_depth);
      if (rebalance) printf("  Rebalancing On\n");
      else printf("  Rebalancing Off\n");
   }

   // read in neuron files
   for (unsigned int in = 0; in < training_filenames.size(); ++in) {
      NeuronData *training_nd = ReadData(training_filenames[in]);
      if (!training_nd) exit(-1);
      training_nds.push_back(training_nd);
   }
   
   if (!ReadDataset()) exit(-1);
   
   // free up memory
   for (unsigned int in = 0; in < training_nds.size(); ++in)
      delete training_nds[in];

   // return success
   return 0;
}
