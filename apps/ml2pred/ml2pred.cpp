// Source file for predicting vessels from machine labels



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
         else if (!output_filename) output_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) {
      fprintf(stderr, "Need to supply input filename.\n");
      return 0;
   }

   // make sure there is an output filename
   if (!output_filename) {
      fprintf(stderr, "Need to supply output filename.\n");
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

   // get input root filename
   char input_root_filename[4096];
   strncpy(input_root_filename, input_filename, 4096);
   char *input_extp = strrchr(input_root_filename, '.');
   *input_extp = '\0';

   // get output root filename
   char output_root_filename[4096];
   strncpy(output_root_filename, output_filename, 4096);
   char *output_extp = strrchr(output_root_filename, '.');
   *output_extp = '\0';

   R3Grid *input_grid = RNReadNeuronMetaRawFile(input_root_filename);
   if (!input_grid) { fprintf(stderr, "Failed to read %s\n", input_filename); exit(-1); }

   R3Grid *output_grid = new R3Grid(input_grid->XResolution(), input_grid->YResolution(), input_grid->ZResolution());

   for (int iv = 0; iv < input_grid->NEntries(); ++iv) {
      int input_value = (int)(input_grid->GridValue(iv) + 0.5);
      if (input_value > 0) output_grid->SetGridValue(iv, 1);
      else output_grid->SetGridValue(iv, 0);
   }

   RNMeta meta = RNMeta("Uint8", output_grid->XResolution(), output_grid->YResolution(), output_grid->ZResolution(), 1);
   if (!RNWriteNeuronMetaRawFile(output_root_filename, meta, output_grid)) exit(-1);

   // free memory
   delete input_grid;
   delete output_grid;

   // return success
   return 0;
}
