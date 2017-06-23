// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static const char *extension = NULL;



// global variables

static std::vector<const char *> filenames = std::vector<const char *>();
static std::vector<NeuronData *> nds = std::vector<NeuronData *>();



// directory structure

static const char *dataset_directory = "algs_data/dataset";



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
         else if (!strcmp(*argv, "-smoothing_term")) extension = "smoothing";
         else if (!strcmp(*argv, "-data_term")) extension = "data";
         else if (!strcmp(*argv, "-global_term")) extension = "smoothing_global";
	 else if (!strcmp(*argv, "-boundary_term")) extension = "data_boundary";
	 else if (!strcmp(*argv, "-boundary_features")) extension = "boundary";
	 else if (!strcmp(*argv, "-shape_features")) extension = "shape";
	 else if (!strcmp(*argv, "-postprocess")) extension = "postprocess";
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         filenames.push_back(*argv);
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!filenames.size()) {
      fprintf(stderr, "Need to supply input filename.\n");
      return 0;
   }

   // make sure there is an extension
   if (!extension) {
      fprintf(stderr, "Must have -smoothing_term, -data_term, or -global_term\n");
      return 0;
   }

   // return OK status 
   return 1;
}



static int ReadDataset(void)
{
   RNDataset *dataset = new RNDataset();

   for (unsigned int in = 0; in < nds.size(); ++in) {
      NeuronData *nd = nds[in];

      // get data root filenames
      char root_filename[4096];
      strncpy(root_filename, nd->Filename(), 4096);
      char *training_extp = strrchr(root_filename, '.');
      *training_extp = '\0';

      // get dataset filenames
      char dataset_filename[4096];
      sprintf(dataset_filename, "%s/%s_%s.mldb", dataset_directory, root_filename, extension);

      // read in the training dataset
      RNDataset *this_dataset = new RNDataset();
      this_dataset->ReadFile(dataset_filename);

      // set the names for this dataset
      if (in == 0) {
         std::vector<std::string> names = this_dataset->GetNames();
         dataset->SetNames(names);
      }

      // copy this training dataset to training dataset
      for (unsigned int id = 0; id < this_dataset->NEntries(); ++id) {
         const std::vector<float> &attributes = this_dataset->GetFeature(id);
         unsigned int label = this_dataset->GetLabel(id);
         unsigned int identification = this_dataset->GetIdentification(id);

         dataset->InsertDatapoint(attributes, label, identification);
      }
      // known memory leak (OK for now)
   }

   // calculate entropy
   dataset->CalculateEntropy();

   // free memory
   delete dataset;

   // return success
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
   for (unsigned int in = 0; in < filenames.size(); ++in) {
      NeuronData *nd = ReadData(filenames[in]);
      if (!nd) exit(-1);
      nds.push_back(nd);
   }

   if (!ReadDataset()) exit(-1);

   // free up memory
   for (unsigned int in = 0; in < nds.size(); ++in)
      delete nds[in];

   // return success
   return 0;
}
