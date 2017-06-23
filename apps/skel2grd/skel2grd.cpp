// Source file for the subsection algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static int zcrop = 0;


// global variables

static NeuronData *nd = NULL;



// directory structure

static const char *grid_directory = "grd_data";



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
         else if (!strcmp(*argv, "-debug")) { print_debug = 1; print_verbose = 1; }
         else if (!strcmp(*argv, "-zcrop")) { argv++; argc--; zcrop = atoi(*argv); }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) { fprintf(stderr, "Need to supply input filename.\n");  return 0; }

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

   // get the root filename
   char root_filename[4096];
   strncpy(root_filename, input_filename, 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   char *sans_path = strrchr(root_filename, '/');
   sans_path++;

   // read skeleton file
   RNSkeleton *skeleton = new RNSkeleton();
   if (!skeleton->ReadSkelFile(input_filename)) exit(-1);

   R3Grid *grid = skeleton->CreateR3Grid(zcrop);
   if (!grid) { fprintf(stderr, "Failed to create R3Grid\n"); exit(-1); }

   // get grid output filename
   char grid_filename[4096];
   sprintf(grid_filename, "%s/%s_skeleton.grd", grid_directory, sans_path);

   // write grid file
   grid->WriteFile(grid_filename);

   // free memory
   delete grid;
   delete skeleton;
   delete nd;

   // return success
   return 0;
}