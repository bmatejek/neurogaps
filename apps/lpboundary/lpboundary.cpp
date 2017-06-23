// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



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

   //////////////////////////////////
   //// CONSTRUCT LINEAR PROGRAM ////
   //////////////////////////////////

   // get the number of interior and boundary supervoxels
   unsigned long ninteriors = 0;
   unsigned long nboundaries = 0;
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      if (cellular->IsOnBoundary()) nboundaries++;
      else ninteriors++;
   }
   ninteriors = 800;
   nboundaries = 400;
   printf("%lu, %lu\n", ninteriors, nboundaries);
   // get the number of data, smoothing, and boundary variables
   unsigned long ndata_variables = ninteriors * nboundaries * nboundaries;
   unsigned long nsmoothing_variables = (ninteriors * (ninteriors - 1) / 2) * nboundaries * nboundaries;
   unsigned long nboundary_variables = nboundaries * nboundaries;
   printf("%lu, %lu, %lu\n", ndata_variables, nsmoothing_variables, nboundary_variables);

   // create all variable names
   //std::string *data_variable_names = new std::string[]

   // free up memory
   delete nd;

   // return success
   return 0;
}
