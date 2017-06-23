// Source file for the subsection algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"
#include <vector>


// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *output_filename = NULL;


// global variables

static std::vector<const char *> meta_files = std::vector<const char *>();



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
         else if (!strcmp(*argv, "-output_filename")) { argv++; argc--; output_filename = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         meta_files.push_back(*argv);
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (meta_files.size() == 0) { fprintf(stderr, "Need to supply input meta files.\n");  return 0; }

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

   R3Grid **grids = new R3Grid *[meta_files.size()];
   for (unsigned int im = 0; im < meta_files.size(); ++im) {
      RNMeta meta = RNReadNeuronMetaFile(meta_files[im]);

      char root_filename[4096];
      strncpy(root_filename, meta_files[im], 4096);
      char *extp = strrchr(root_filename, '.');
      *extp = '\0';

      if (meta.NGrids() == 1) {
         grids[im] = RNReadNeuronMetaRawFile(root_filename);
      }
      else {
         grids[im] = RNReadNeuronMetaRawFile(root_filename, TRUE)[0];
      }
   }

   R3Grid *merged_grid = new R3Grid(grids[0]->XResolution() * 3, grids[0]->YResolution() * 2, grids[0]->ZResolution() * 2);
   for (unsigned int im = 0; im < meta_files.size(); ++im) {
      R3Grid *grid = grids[im];

      int xoffset = 0;
      int yoffset = 0;
      if (im == 1 || im == 4) xoffset = 256;
      if (im == 2 || im == 5) xoffset = 512;
      
      if (im > 2) yoffset = 256;

      RNScalar max_value = grid->Maximum();

      for (int ix = 0; ix < grid->XResolution(); ++ix) {
         for (int iy = 0; iy < grid->YResolution(); ++iy) {
            for (int iz = 0; iz < grid->ZResolution(); ++iz) {
               RNScalar grid_value;
               if (im != 4) grid_value = grid->GridValue(ix, iy, iz) / max_value;
               else grid_value = 0.2 * grid->GridValue(ix, iy, iz) / max_value + 0.8;
               merged_grid->SetGridValue(ix + xoffset, iy + yoffset, iz, grid_value);
            }
         }
      }
   }

   RNMeta meta = RNMeta("Float32", merged_grid->XResolution(), merged_grid->YResolution(), merged_grid->ZResolution(), 1);
   if (!RNWriteNeuronMetaRawFile(output_filename, meta, merged_grid)) exit(-1);
   

   // free memory
   for (unsigned int im = 0; im < meta_files.size(); ++im) {
      delete grids[im];
   }
   delete[] grids;
   delete merged_grid;

   // return success
   return 0;
}