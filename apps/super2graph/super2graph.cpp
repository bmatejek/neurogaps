// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *algs_directory = "algs_data/laplacian";



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

   int ncellulars = nd->NCellulars();

   // create binary adjacency matrix
   RNScalar **adjacency_matrix = new RNScalar *[ncellulars];
   for (int ic1 = 0; ic1 < ncellulars; ++ic1) {
      adjacency_matrix[ic1] = new RNScalar[ncellulars];
      for (int ic2 = 0; ic2 < ncellulars; ++ic2) {
         adjacency_matrix[ic1][ic2] = 0.0;
      }
   }

   // go through all cellulars
   for (int ic = 0; ic < ncellulars; ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = cellular->Boundary(ib);
         NeuronSupervoxel *supervoxel = boundary->OtherSupervoxel(cellular);
         if (supervoxel->IsExtracellular()) continue;
         NeuronCellular *neighbor = (NeuronCellular *)supervoxel;

         int cellular_data_index = cellular->DataIndex();
         int neighbor_data_index = neighbor->DataIndex();

         RNScalar boundary_affinity = boundary->Mean();
         adjacency_matrix[cellular_data_index][neighbor_data_index] = boundary_affinity;
      }
   }

   // create diagonal matrix
   RNScalar **diagonal_matrix = new RNScalar *[ncellulars];
   for (int ic1 = 0; ic1 < ncellulars; ++ic1) {
      diagonal_matrix[ic1] = new RNScalar[ncellulars];
      for (int ic2 = 0; ic2 < ncellulars; ++ic2) {
         diagonal_matrix[ic1][ic2] = 0;
      }
   }

   // add in diagonal entries
   for (int i = 0; i < ncellulars; ++i) {
      RNScalar diagonal_sum = 0.0;
      for (int j = 0; j < ncellulars; ++j) {
         diagonal_sum += adjacency_matrix[i][j];
      }
      diagonal_matrix[i][i] = diagonal_sum;
   }

   // create laplacian matrix
   RNScalar **laplacian_matrix = new RNScalar *[ncellulars];
   for (int ic = 0; ic < ncellulars; ++ic) {
      laplacian_matrix[ic] = new RNScalar[ncellulars];
   }

   // add in all entries
   for (int i = 0; i < ncellulars; ++i) {
      for (int j = 0; j < ncellulars; ++j) {
         laplacian_matrix[i][j] = diagonal_matrix[i][j] - adjacency_matrix[i][j];
      }
   }

   // output laplacian matrix
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create the filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_laplacian.graph", algs_directory, root_filename);

   // open file
   FILE *fp = fopen(output_filename, "wb");
   if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   fwrite(&ncellulars, sizeof(int), 1, fp);

   for (int i = 0; i < ncellulars; ++i) {
      for (int j = 0; j < ncellulars; ++j) {
         fwrite(&(laplacian_matrix[i][j]), sizeof(RNScalar), 1, fp);
      }
   }

   // close file
   fclose(fp);

   // free up memory
   for (int ic = 0; ic < ncellulars; ++ic) {
      delete[] adjacency_matrix[ic];
      delete[] diagonal_matrix[ic];
      delete[] laplacian_matrix[ic];
   }
   
   delete nd;
   delete[] adjacency_matrix;
   delete[] diagonal_matrix;
   delete[] laplacian_matrix;

   // return success
   return 0;
}
