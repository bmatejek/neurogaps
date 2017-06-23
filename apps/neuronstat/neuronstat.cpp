// Source file for the neuron statistics algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int watershed = 0;
static int boundaries = 0;
static int affinities = 0;
static int cellulars = 0;
static int voxels = 0;



// global variables

static NeuronData *nd = NULL;
static char *input_filename = NULL;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadData(const char *filename)
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
// General statistic functions
////////////////////////////////////////////////////////////////////////

static int AffinityStatistics(void)
{
   // read in voxels
   nd->ReadVoxels();

   RNScalar mean[4];
   RNScalar stddev[4];
   nd->AffinityStatistics(mean, stddev);

   printf("------------------------------------------------------\n");
   printf("Affinity Statistics\n");
   printf("------------------------------------------------------\n");
   printf("X Mean:   %0.6lf\n", mean[RN_X]);
   printf("X StdDev: %0.6lf\n", stddev[RN_X]);
   printf("Y Mean:   %0.6lf\n", mean[RN_Y]);
   printf("Y StdDev: %0.6lf\n", stddev[RN_Y]);
   printf("Z Mean:   %0.6lf\n", mean[RN_Z]);
   printf("Z StdDev: %0.6lf\n", stddev[RN_Z]);
   printf("Mean:     %0.6lf\n", mean[3]);
   printf("StdDev:   %0.6lf\n", stddev[3]);
   printf("------------------------------------------------------\n");

   // return success
   return 1;
}



static int WatershedStatistics(void)
{
   // read in voxels
   nd->ReadVoxels();

   // get cellular and extracellular results
   unsigned int true_positives = 0;
   unsigned int false_positives = 0;
   unsigned int false_negatives = 0;
   unsigned int true_negatives = 0;

   // iterate through all voxels
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel()) {
         if (voxel->IsCellular()) true_positives++;
         else false_negatives++;
      }
      else {
         if (voxel->IsCellular()) false_positives++;
         else true_negatives++;
      }
   }

   // print out confusion matrix
   printf("------------------------------------------------------\n");
   printf("Watershed Statistics\n");
   printf("------------------------------------------------------\n");
   printf("Confusion Matrix:\n");
   printf("                  Predicted Class\n");
   printf("Actual Class");
   printf("    Cellular    Extracellular\n");
   printf("     Cellular %10u       %10u\n", true_positives, false_negatives);
   printf("Extracellular %10u       %10u\n", false_positives, true_negatives);

   // print out sums
   printf("------------------------------------------------------\n");
   printf("Cellular correctly labeled: %lf\n", true_positives / (RNScalar)(true_positives + false_negatives));
   printf("Extracellular correctly labeled: %lf\n", true_negatives / (RNScalar)(false_positives + true_negatives));
   printf("Total correctly labeled: %lf\n", (true_negatives + true_positives) / (RNScalar)(true_positives + true_negatives + false_positives + false_negatives));
   printf("------------------------------------------------------\n");

   // return success
   return 1;
}



static int BoundaryStatistics(void)
{
   // print the number of boundaries that should oversegment
   int oversegmented_boundaries = 0;
   int correct_boundaries = 0;
   int extracellular_boundaries = 0;

   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      if (supervoxel_one->IsExtracellular() || supervoxel_two->IsExtracellular()) {
         extracellular_boundaries++;
         continue;
      }

      if (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel()) {
         oversegmented_boundaries++;
      }
      else {
         correct_boundaries++;
      }
   }

   // print out stats
   printf("------------------------------------------------------\n");
   printf("Boundary Statistics\n");
   printf("------------------------------------------------------\n");
   printf("Number Correct Boundaries: %d (%0.2f%%)\n", correct_boundaries, 100.0 * correct_boundaries / (RNScalar)nd->NBoundaries());
   printf("Number Oversegmented Boundaries: %d (%0.2f%%)\n", oversegmented_boundaries, 100.0 * oversegmented_boundaries / (RNScalar)nd->NBoundaries());
   printf("Number Cellular-Extracellular Boundaries: %d (%0.2f%%)\n", extracellular_boundaries, 100.0 * extracellular_boundaries / (RNScalar)nd->NBoundaries());
   printf("Total Number of Boundaries: %d\n", nd->NBoundaries());

   printf("------------------------------------------------------\n");
   printf("Without Cellular-Extracellular Boundaries\n");
   printf("------------------------------------------------------\n");
   printf("Number of Correct Boundaries: %d (%0.2f%%)\n", correct_boundaries, 100 * correct_boundaries / (RNScalar)(correct_boundaries + oversegmented_boundaries));
   printf("Number of Oversegmented Boundaries: %d (%0.2f%%)\n", oversegmented_boundaries, 100 * oversegmented_boundaries / (RNScalar)(correct_boundaries + oversegmented_boundaries));
   printf("Total Number of Boundaries: %d\n", correct_boundaries + oversegmented_boundaries);
   printf("------------------------------------------------------\n");

   // return success
   return 1;
}



static int CellularStatistics(void)
{
   int nboundary_cellulars = 0;
   int nhuman_labels = 0;
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      if (cellular->IsOnBoundary()) nboundary_cellulars++;
      nhuman_labels += cellular->NHumanLabels();
   }

   // print out stats
   printf("------------------------------------------------------\n");
   printf("Cellular Statistics\n");
   printf("------------------------------------------------------\n");
   printf("Number Boundary Cellulars: %d (%0.2f%%)\n", nboundary_cellulars, 100.0 * nboundary_cellulars / (RNScalar)nd->NCellulars());
   printf("Number Interior Cellulars: %d (%0.2f%%)\n", (nd->NCellulars() - nboundary_cellulars), 100.0 * (nd->NCellulars() - nboundary_cellulars) / (RNScalar)nd->NCellulars());
   printf("------------------------------------------------------\n");
   printf("Average Voxels per Cellular: %0.2f\n", nd->NVoxels() / (RNScalar)nd->NCellulars());
   printf("Average Boundaries per Cellular: %0.2f\n", nd->NBoundaries() / (RNScalar)nd->NCellulars());
   printf("Average Human Labels per Cellular: %0.2f\n", nhuman_labels / (RNScalar)nd->NCellulars());
   printf("------------------------------------------------------\n");

   // return success
   return 1;
}



static int VoxelStatistics(void)
{
   // read voxels
   if (!nd->AreVoxelsResident()) nd->ReadVoxels();

   int ncellular_successes = 0;
   int ncellular_failures = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (!voxel->HumanLabel()) continue;

      NeuronSupervoxel *supervoxel = voxel->Supervoxel();
      if (voxel->HumanLabel() == supervoxel->MajorityHumanLabel()) ncellular_successes++;
      else ncellular_failures++;
   }

   // print out stats
   printf("------------------------------------------------------\n");
   printf("Voxel Statistics\n");
   printf("------------------------------------------------------\n");
   printf("Voxel Cellular Assignment Success: %d (%0.2f%%)\n", ncellular_successes, 100.0 * ncellular_successes / (RNScalar)(ncellular_successes + ncellular_failures));
   printf("Voxel Cellular Assignment Failure: %d (%0.2f%%)\n", ncellular_failures, 100.0 * ncellular_failures / (RNScalar)(ncellular_successes + ncellular_failures));
   printf("------------------------------------------------------\n");

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
         else if (!strcmp(*argv, "-affinities")) affinities = 1;
         else if (!strcmp(*argv, "-watershed")) watershed = 1;
         else if (!strcmp(*argv, "-boundaries")) boundaries = 1;
         else if (!strcmp(*argv, "-cellulars")) cellulars = 1;
         else if (!strcmp(*argv, "-voxels")) voxels = 1;
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

int
main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   ///////////////////////////////////
   //// RUN ALL OF THE STATISTICS ////
   ///////////////////////////////////

   if (affinities && !AffinityStatistics()) exit(-1);
   if (watershed && !WatershedStatistics()) exit(-1);
   if (boundaries && !BoundaryStatistics()) exit(-1);
   if (cellulars && !CellularStatistics()) exit(-1);
   if (voxels && !VoxelStatistics()) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}