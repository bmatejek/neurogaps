// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 1;
static int print_debug = 1;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *output_directory = "filters";



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


struct AffinityTruthPair {
   RNScalar affinity;
   RNBoolean truth;
};



int AffinityTruthSort(AffinityTruthPair a, AffinityTruthPair b)
{
   return a.affinity < b.affinity;
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

static void CopyTmp(R3Grid *grid, R3Grid *tmp)
{
   for (int iv = 0; iv < grid->NEntries(); ++iv)
      tmp->SetGridValue(iv, grid->GridValue(iv));
}



static RNScalar Smoothing(R3Grid *tmp, int iv, int dim)
{
   // get coordinates
   int ix, iy, iz;
   nd->IndexToIndices(iv, ix, iy, iz);

   RNScalar previous_value;
   RNScalar current_value;
   RNScalar next_value;

   if (dim == RN_X) {
      if (ix == 0) previous_value = tmp->GridValue(ix + 1, iy, iz);
      else previous_value = tmp->GridValue(ix - 1, iy, iz);

      current_value = tmp->GridValue(ix, iy, iz);

      if (ix == nd->XResolution() - 1) next_value = tmp->GridValue(ix - 1, iy, iz);
      else next_value = tmp->GridValue(ix + 1, iy, iz);
   }
   else if (dim == RN_Y) {
      if (iy == 0) previous_value = tmp->GridValue(ix, iy + 1, iz);
      else previous_value = tmp->GridValue(ix, iy - 1, iz);

      current_value = tmp->GridValue(ix, iy, iz);

      if (iy == nd->YResolution() - 1) next_value = tmp->GridValue(ix, iy - 1, iz);
      else next_value = tmp->GridValue(ix, iy + 1, iz);
   }
   else if (dim == RN_Z) {
      if (iz == 0) previous_value = tmp->GridValue(ix, iy, iz + 1);
      else previous_value = tmp->GridValue(ix, iy, iz - 1);

      current_value = tmp->GridValue(ix, iy, iz);

      if (iz == nd->ZResolution() - 1) next_value = tmp->GridValue(ix, iy, iz - 1);
      else next_value = tmp->GridValue(ix, iy, iz + 1);
   }
   else rn_assertion(FALSE);

   // return the smoothing value at this location
   return previous_value + 2 * current_value + next_value;
}



static RNScalar Difference(R3Grid *tmp, int iv, int dim)
{
   // get coordinates
   int ix, iy, iz;
   nd->IndexToIndices(iv, ix, iy, iz);

   RNScalar previous_value;
   RNScalar next_value;

   if (dim == RN_X) {
      if (ix == 0) previous_value = tmp->GridValue(ix + 1, iy, iz);
      else previous_value = tmp->GridValue(ix - 1, iy, iz);

      if (ix == nd->XResolution() - 1) next_value = tmp->GridValue(ix - 1, iy, iz);
      else next_value = tmp->GridValue(ix + 1, iy, iz);
   }
   else if (dim == RN_Y) {
      if (iy == 0) previous_value = tmp->GridValue(ix, iy + 1, iz);
      else previous_value = tmp->GridValue(ix, iy - 1, iz);

      if (iy == nd->YResolution() - 1) next_value = tmp->GridValue(ix, iy - 1, iz);
      else next_value = tmp->GridValue(ix, iy + 1, iz);
   }
   else if (dim == RN_Z) {
      if (iz == 0) previous_value = tmp->GridValue(ix, iy, iz + 1);
      else previous_value = tmp->GridValue(ix, iy, iz - 1);

      if (iz == nd->ZResolution() - 1) next_value = tmp->GridValue(ix, iy, iz - 1);
      else next_value = tmp->GridValue(ix, iy, iz + 1);
   }
   else rn_assertion(FALSE);

   // return the smoothing value at this location
   return previous_value - next_value;
}



static void CalculateDerivative(R3Grid *starting[3], R3Grid *gradient[3])
{
   R3Grid *xtmp = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   R3Grid *ytmp = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   R3Grid *ztmp = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());

   // calculate derivative at each point for each dimension
   for (int dim = 0; dim <= 2; ++dim) {
      // calculate x derivative
      R3Grid *xaffinity = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
      R3Grid *yaffinity = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
      R3Grid *zaffinity = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());

      // set the initial affinities
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         xaffinity->SetGridValue(iv, starting[dim]->GridValue(iv));
         yaffinity->SetGridValue(iv, starting[dim]->GridValue(iv));
         zaffinity->SetGridValue(iv, starting[dim]->GridValue(iv));
      }
      CopyTmp(xaffinity, xtmp);
      CopyTmp(yaffinity, ytmp);
      CopyTmp(zaffinity, ztmp);

      // apply x filters
      if (print_verbose) { printf("Applying x filters for dim %d\n  ", dim); fflush(stdout); }
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (print_verbose) RNProgressBar(iv, nd->NVoxels());
         xaffinity->SetGridValue(iv, Difference(xtmp, iv, RN_X));
         yaffinity->SetGridValue(iv, Smoothing(ytmp, iv, RN_X));
         zaffinity->SetGridValue(iv, Smoothing(ztmp, iv, RN_X));
      }
      CopyTmp(xaffinity, xtmp);
      CopyTmp(yaffinity, ytmp);
      CopyTmp(zaffinity, ztmp);
      if (print_verbose) printf("\ndone.\n");

      // apply y filters
      if (print_verbose) { printf("Applying y filters for dim %d\n  ", dim); fflush(stdout); }
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (print_verbose) RNProgressBar(iv, nd->NVoxels());
         xaffinity->SetGridValue(iv, Smoothing(xtmp, iv, RN_Y));
         yaffinity->SetGridValue(iv, Difference(ytmp, iv, RN_Y));
         zaffinity->SetGridValue(iv, Smoothing(ztmp, iv, RN_Y));
      }
      CopyTmp(xaffinity, xtmp);
      CopyTmp(yaffinity, ytmp);
      CopyTmp(zaffinity, ztmp);
      if (print_verbose) printf("\ndone.\n");

      // apply z filters
      if (print_verbose) { printf("Applying z filters for dim %d\n  ", dim); fflush(stdout); }
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (print_verbose) RNProgressBar(iv, nd->NVoxels());
         xaffinity->SetGridValue(iv, Smoothing(xtmp, iv, RN_Z));
         yaffinity->SetGridValue(iv, Smoothing(ytmp, iv, RN_Z));
         zaffinity->SetGridValue(iv, Difference(ztmp, iv, RN_Z));
      }
      if (print_verbose) printf("\ndone.\n");

      // calculate the gradient
      if (print_verbose) { printf("Applying gradients for dim %d\n  ", dim); fflush(stdout); }
      gradient[dim] = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (print_verbose) RNProgressBar(iv, nd->NVoxels());
         RNScalar xsquared = xaffinity->GridValue(iv) * xaffinity->GridValue(iv);
         RNScalar ysquared = yaffinity->GridValue(iv) * yaffinity->GridValue(iv);
         RNScalar zsquared = zaffinity->GridValue(iv) * zaffinity->GridValue(iv);

         gradient[dim]->SetGridValue(iv, sqrt(xsquared + ysquared + zsquared));
      }
      if (print_verbose) printf("\ndone.\n");

      // free memory
      delete xaffinity;
      delete yaffinity;
      delete zaffinity;
   }

   // free memory
   delete xtmp;
   delete ytmp;
   delete ztmp;
}



int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // get the affinities
   R3Grid *starting_affinities[3] = { NULL, NULL, NULL };
   R3Grid *first_derivatives[3] = { NULL, NULL, NULL };
   R3Grid *second_derivatives[3] = { NULL, NULL, NULL };

   for (int dim = 0; dim <= 2; ++dim) {
      starting_affinities[dim] = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         starting_affinities[dim]->SetGridValue(iv, nd->Voxel(iv)->Affinity(dim));
      }
   }

   // calculate the derivative
   CalculateDerivative(starting_affinities, first_derivatives);
   CalculateDerivative(first_derivatives, second_derivatives);

   // save the gradient
   RNMeta meta = RNMeta("Float32", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 3);

   char first_derivative_filename[4096];
   sprintf(first_derivative_filename, "%s/%s_first_derivative", output_directory, root_filename);

   // output file
   if (!RNWriteNeuronMetaRawFile(first_derivative_filename, meta, first_derivatives)) exit(-1);

   char second_derivative_filename[4096];
   sprintf(second_derivative_filename, "%s/%s_second_derivative", output_directory, root_filename);

   if (!RNWriteNeuronMetaRawFile(second_derivative_filename, meta, second_derivatives)) exit(-1);

   // free up memory
   delete nd;
   for (int dim = 0; dim <= 2; ++dim) {
      delete starting_affinities[dim];
      delete first_derivatives[dim];
      delete second_derivatives[dim];
   }

   // return success
   return 0;
}
