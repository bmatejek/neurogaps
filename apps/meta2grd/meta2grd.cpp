// Source file for the subsection algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static char *input_filename = NULL;
static char *output_filename = NULL;


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
         // number of files to generate
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else if (!output_filename) output_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) { fprintf(stderr, "Need to supply input filename.\n");  return 0; }
   if (!output_filename) { fprintf(stderr, "Need to supply output filename.\n"); return 0; }

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

   // find out dimensionality
   RNMeta meta = RNReadNeuronMetaFile(input_filename);

   if (meta.NGrids() == 1) { 
      R3Grid *grid = RNReadNeuronMetaRawFile(root_filename);
      grid->WriteFile(output_filename);
   }
   else {
      char output_root_filename[4096];
      strncpy(output_root_filename, output_filename, 4096);
      char *extp = strrchr(output_root_filename, '.');
      *extp = '\0';

      R3Grid **grids = RNReadNeuronMetaRawFile(root_filename, TRUE);
      for (int id = 0; id < meta.NGrids(); ++id) {
         char dim_output_filename[4096];
         sprintf(dim_output_filename, "%s_%02d.grd", output_root_filename, id);
         grids[id]->WriteFile(dim_output_filename);
      }
   }

   // return success
   return 0;
}