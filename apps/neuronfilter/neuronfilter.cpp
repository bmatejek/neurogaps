// source file for the affinity filter functions



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_debug = 0;
static int print_verbose = 0;
static int maximum_filter = 0;
static int mean_filter = 0;
static int median_filter = 0;
static int minimum_filter = 0;
static int minmax_filter = 0;
static int normalize_filter = 0;
static unsigned int minmax[2] = { 0, 0 };



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

static int ParseArgs(int argc, char **argv) 
{
    // parse arguments
    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-v")) print_verbose = 1;
            else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
            // affinity changes
            else if (!strcmp(*argv, "-maximum")) maximum_filter = 1;
            else if (!strcmp(*argv, "-mean")) mean_filter = 1;
            else if (!strcmp(*argv, "-median")) median_filter = 1;
            else if (!strcmp(*argv, "-minimum")) minimum_filter = 1;
            else if (!strcmp(*argv, "-minmax")) { 
                minmax_filter = 1;
                argv++; argc--; minmax[0] = atoi(*argv); 
                argv++; argc--; minmax[1] = atoi(*argv); 
            }
            // normalize is not exclusive
            else if (!strcmp(*argv, "-normalize")) normalize_filter = 1;
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

    // make sure that at most one option is selected
    if ((minmax_filter + maximum_filter + mean_filter + median_filter + minimum_filter) != 1 && normalize_filter != 1) {
        fprintf(stderr, "Please select one filter type (maximum, mean, median, minimum, minmax, gaussian, diffusion, normalize)\n");
        return 0;
    }

    // make sure that at least one option is selected
    if ((minmax_filter + maximum_filter + mean_filter + median_filter + minimum_filter) == 1 && normalize_filter == 1) {
        printf("Normalized filter will be applied after other filters\n");
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
    
    // read in file
    if (!ReadData(input_filename)) exit(-1);    

    
    ////////////////////////////////////////
    //// Boundary affinity manipulation ////
    ////////////////////////////////////////

    // run maximum affinity filter
    if (maximum_filter) {
        if (!nd->FilterAffinity("maximum", ".txt")) exit(-1);
    }
    // run mean affinity filter
    else if (mean_filter) {
        if (!nd->FilterAffinity("mean", ".txt")) exit(-1);
    }
    // run median affinity filter
    else if (median_filter) {
        if (!nd->FilterAffinity("median", ".txt")) exit(-1);
    }
    // run minimum affinity filter
    else if (minimum_filter) {
        if (!nd->FilterAffinity("minimum", ".txt")) exit(-1);
    }

    // run normalize affinity filter
    if (normalize_filter) {
        if (!nd->FilterAffinity("normalize", ".txt")) exit(-1);
    }

    // free up memory
    delete nd;

    // return success
	return 0;
}