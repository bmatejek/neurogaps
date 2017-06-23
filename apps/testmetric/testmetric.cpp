// Source file for the test metric algorithm



// include files 

#include "Neuron/Neuron.h"
#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static char *proposal_filename = NULL;



// global variables

static NeuronData *nd = NULL;



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

static int
ParseArgs(int argc, char **argv)
{
   // Parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-proposals")) { argc--; argv++; proposal_filename = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) { input_filename = *argv; }
         else if (!proposal_filename) proposal_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is a neuron file
   if (!input_filename) { fprintf(stderr, "Must supply an input filename\n"); return 0; }

   // make sure there is a proposed filename
   if (!proposal_filename) { fprintf(stderr, "Must supply a proposal filename\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
   // parse arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in neuron data
   if (!ReadData(input_filename)) exit(-1);

   // read in proposal grid
   char *proposal_extp = strrchr(proposal_filename, '.');
   if (proposal_extp) *proposal_extp = '\0';
   R3Grid *proposal_grid = RNReadNeuronMetaRawFile(proposal_filename);
   if (!proposal_grid) exit(-1);

   int *voxel_proposals = new int[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      voxel_proposals[iv] = (int)(proposal_grid->GridValue(iv) + 0.5);
   }
   int *cellular_proposals = new int[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      int voxel_index = nd->Cellular(ic)->CenterVoxel()->DataIndex();
      cellular_proposals[ic] = (int)(proposal_grid->GridValue(voxel_index) + 0.5);
   }

   // get the rand error for this proposal grid
   int ntest_metrics = nd->NTestMetrics();
   RNScalar *test_metrics = new RNScalar[ntest_metrics];
   if (!nd->TestMetric(voxel_proposals, test_metrics)) return 0;

   // get the segmentation metrics
   int nsegmentation_metrics = nd->NSegmentationMetrics();
   RNScalar *segmentation_metrics = new RNScalar[nsegmentation_metrics];
   if (!nd->SegmentationMetric(cellular_proposals, segmentation_metrics)) return 0;

   // free memory
   delete[] voxel_proposals;
   delete[] cellular_proposals;
   delete[] test_metrics;
   delete[] segmentation_metrics;
   delete proposal_grid;

   // return success
   return 0;
}
