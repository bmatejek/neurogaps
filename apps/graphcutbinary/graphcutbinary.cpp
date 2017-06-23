// Source file to create graph cut initial positions



// include files 

#include "Neuron/Neuron.h"



// program arguments

static char *input_name = NULL;
static const char *binary_extension = NULL;
static int print_verbose = 0;
static int print_debug = 0;



// program variables

static NeuronData *nd;



// directory variables

static const char *algs_directory = "algs_data/graphcut";



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
   if (!nd->ReadFile(filename, FALSE)) {
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
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-truth")) { binary_extension = "truth"; }
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input argument
   if (!input_name) { fprintf(stderr, "Need to supply input filename.\n"); return 0; }
   if (!binary_extension) { fprintf(stderr, "Need to provide metric.\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
   // Parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in the neuron data
   if (!ReadData(input_name)) exit(-1);

   // create all memory
   float *smoothing_terms = new float[nd->NBoundaries()];
   for (int ib = 0; ib < nd->NBoundaries(); ++ib)
      smoothing_terms[ib] = FLT_MAX;

   // go through every cellular boundary and determining the smoothing term
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get supervoxels
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_one->IsExtracellular() || supervoxel_two->IsExtracellular()) continue;

      if (!strncmp(binary_extension, "truth", 5)) {
         NeuronHumanLabel *human_label_one = supervoxel_one->MajorityHumanLabel();
         NeuronHumanLabel *human_label_two = supervoxel_two->MajorityHumanLabel();

         smoothing_terms[ib] = (human_label_one == human_label_two);
      }
      else { fprintf(stderr, "Unrecognized extension: %s\n", binary_extension); }
   }
  
   // get neuron root filename
   char root_filename[4096];
   strcpy(root_filename, nd->Filename());
   char *extp = strchr(root_filename, '.');
   *extp = '\0';

   // get binary filename
   char binary_filename[4096];
   sprintf(binary_filename, "%s/%s_%s.binary", algs_directory, root_filename, binary_extension);

   // open binary file for writing
   FILE *binary_fp = fopen(binary_filename, "wb");
   if (!binary_fp) { fprintf(stderr, "Failed to open %s\n", binary_filename); exit(-1); }

   // write the smoothing terms
   int binary_nboundaries = nd->NBoundaries();
   fwrite(&binary_nboundaries, sizeof(int), 1, binary_fp);
   fwrite(smoothing_terms, sizeof(float), nd->NBoundaries(), binary_fp);

   // close file
   fclose(binary_fp);

   // free memory
   delete[] smoothing_terms;
   delete nd;

   // return success
   return 0;
}
