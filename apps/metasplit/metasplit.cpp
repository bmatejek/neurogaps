// Source file for the subsection algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
// proportion along all dimensions
static int xsplit = 0;
static int ysplit = 0;
static int zsplit = 0;



// global variables

static char *input_filename = NULL;



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
         // scale changes
         else if (!strcmp(*argv, "-x")) { argv++; argc--; xsplit = atoi(*argv); }
         else if (!strcmp(*argv, "-y")) { argv++; argc--; ysplit = atoi(*argv); }
         else if (!strcmp(*argv, "-z")) { argv++; argc--; zsplit = atoi(*argv); }
         // number of files to generate
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

   // make sure each split integer is define
   if (!xsplit || !ysplit || !zsplit) { fprintf(stderr, "Must define the number of splits for each dimensions\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

static int
OneDimensionalSplit(RNMeta input_meta, const char root_filename[4096])
{
   // read meta file
   R3Grid *input_grid = RNReadNeuronMetaRawFile(root_filename);

   // make sure that the file is divisible
   rn_assertion(input_grid->XResolution() % xsplit == 0);
   rn_assertion(input_grid->YResolution() % ysplit == 0);
   rn_assertion(input_grid->ZResolution() % zsplit == 0);

   // get size of each dimension
   int nx = input_grid->XResolution() / xsplit;
   int ny = input_grid->YResolution() / ysplit;
   int nz = input_grid->ZResolution() / zsplit;

   int nfiles = 1;
   for (int ix = 0; ix < input_grid->XResolution(); ix += nx) {
      for (int iy = 0; iy < input_grid->YResolution(); iy += ny) {
         for (int iz = 0; iz < input_grid->ZResolution(); iz += nz) {
            // get dimensions
            int xmin = ix;
            int ymin = iy;
            int zmin = iz;
            int xmax = ix + nx - 1;
            int ymax = iy + ny - 1;
            int zmax = iz + nz - 1;

            // create filename
            char output_filename[4096];
            sprintf(output_filename, "%s%d", root_filename, nfiles);

            // file widths
            int xwidth = xmax - xmin + 1;
            int ywidth = ymax - ymin + 1;
            int zwidth = zmax - zmin + 1;
            R3Grid *output_grid = new R3Grid(xwidth, ywidth, zwidth);

            int newi = 0;
            for (int oldi = xmin; oldi <= xmax; ++oldi, ++newi) {
               int newj = 0;
               for (int oldj = ymin; oldj <= ymax; ++oldj, ++newj) {
                  int newk = 0;
                  for (int oldk = zmin; oldk <= zmax; ++oldk, ++newk) {
                     output_grid->SetGridValue(newi, newj, newk, input_grid->GridValue(oldi, oldj, oldk));
                  }
               }
            }

            // write the meta file
            RNMeta meta = RNMeta(input_meta.DataType(), xwidth, ywidth, zwidth, 1);
            RNWriteNeuronMetaRawFile(output_filename, meta, output_grid);

            // free memory
            delete output_grid;
            nfiles++;
         }
      }
   }

   // free memory
   delete input_grid;

   // return success
   return 1;
}



static int
ThreeDimensionalSplit(RNMeta input_meta, const char root_filename[4096])
{
   // read meta file
   R3Grid **input_grid = RNReadNeuronMetaRawFile(root_filename, TRUE);

   // make sure that the file is divisible
   rn_assertion(input_grid[0]->XResolution() % xsplit == 0);
   rn_assertion(input_grid[0]->YResolution() % ysplit == 0);
   rn_assertion(input_grid[0]->ZResolution() % zsplit == 0);

   // get size of each dimension
   int nx = input_grid[0]->XResolution() / xsplit;
   int ny = input_grid[0]->YResolution() / ysplit;
   int nz = input_grid[0]->ZResolution() / zsplit;

   int nfiles = 1;
   for (int ix = 0; ix < input_grid[0]->XResolution(); ix += nx) {
      for (int iy = 0; iy < input_grid[0]->YResolution(); iy += ny) {
         for (int iz = 0; iz < input_grid[0]->ZResolution(); iz += nz) {
            // get dimensions
            int xmin = ix;
            int ymin = iy;
            int zmin = iz;
            int xmax = ix + nx - 1;
            int ymax = iy + ny - 1;
            int zmax = iz + nz - 1;

            // create filename
            char output_filename[4096];
            sprintf(output_filename, "%s%d", root_filename, nfiles);

            // file widths
            int xwidth = xmax - xmin + 1;
            int ywidth = ymax - ymin + 1;
            int zwidth = zmax - zmin + 1;
            R3Grid **output_grid = new R3Grid *[3];
            for (int dim = 0; dim <= 2; ++dim) {
               output_grid[dim] = new R3Grid(xwidth, ywidth, zwidth);
            }

            for (int dim = 0; dim <= 2; ++dim) {
               int newi = 0;
               for (int oldi = xmin; oldi <= xmax; ++oldi, ++newi) {
                  int newj = 0;
                  for (int oldj = ymin; oldj <= ymax; ++oldj, ++newj) {
                     int newk = 0;
                     for (int oldk = zmin; oldk <= zmax; ++oldk, ++newk) {
                        output_grid[dim]->SetGridValue(newi, newj, newk, input_grid[dim]->GridValue(oldi, oldj, oldk));
                     }
                  }
               }
            }

            // write the meta file
            RNMeta meta = RNMeta(input_meta.DataType(), xwidth, ywidth, zwidth, 3);
            RNWriteNeuronMetaRawFile(output_filename, meta, output_grid);

            // free memory
            for (int dim = 0; dim <= 2; ++dim)
               delete output_grid[dim];
            nfiles++;
         }
      }
   }

   // free memory
   for (int dim = 0; dim <= 2; ++dim) {
      delete input_grid[dim];
   }
   delete input_grid;

   // return success
   return 1;
}



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

   if (meta.NGrids() == 1) { OneDimensionalSplit(meta, root_filename); }
   else { ThreeDimensionalSplit(meta, root_filename); }

   // return success
   return 0;
}