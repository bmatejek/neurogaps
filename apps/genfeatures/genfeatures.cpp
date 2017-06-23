// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int prediction = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *features_directory = "algs_data/features";
static const char *postprocess_directory = "postprocess";



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



////////////////////////////////////////////////////////////////////////
// Feature creation functions
////////////////////////////////////////////////////////////////////////

static int
CreateRandomWalkFeatures(const char root_filename[4096])
{
   // see if feature files already exist
   RNBoolean found_random_walk = FALSE;
   RNBoolean found_drunkards_walk = FALSE;

   // get random walk file
   char random_walk_cache_filename[4096];
   sprintf(random_walk_cache_filename, "%s/%s_random_walk.feature", features_directory, root_filename);

   // see if the file exists to read
   FILE *random_walk_test_fp = fopen(random_walk_cache_filename, "rb");
   if (random_walk_test_fp) {
      found_random_walk = TRUE;
      fclose(random_walk_test_fp);
   }

   // get drunkards walk file
   char drunkards_walk_cache_filename[4096];
   sprintf(drunkards_walk_cache_filename, "%s/%s_drunkards_walk.feature", features_directory, root_filename);

   // see if the file exists to read
   FILE *drunkards_walk_test_fp = fopen(drunkards_walk_cache_filename, "rb");
   if (drunkards_walk_test_fp) {
      found_drunkards_walk = TRUE;
      fclose(drunkards_walk_test_fp);
   }

   if (!found_random_walk) {
      // calculate walk
      RNScalar *random_walk = new RNScalar[nd->NBoundaries()];
      if (print_verbose) { printf("Createing random walk boundary features...\n   "); fflush(stdout); }
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         if (print_verbose) RNProgressBar(ib, nd->NBoundaries());
         NeuronBoundary *boundary = nd->Boundary(ib);
         random_walk[ib] = boundary->RandomWalk(FALSE);
      }
      if (print_verbose) printf("\ndone!\n");

      // open file for writing
      FILE *random_fp = fopen(random_walk_cache_filename, "wb");
      if (!random_fp) { fprintf(stderr, "Failed to write to %s\n", random_walk_cache_filename); return 0; }

      // write the number of boundaries
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, random_fp);

      // write the random walk to data
      fwrite(random_walk, sizeof(RNScalar), nd->NBoundaries(), random_fp);

      // close files
      fclose(random_fp);

      // free memory
      delete[] random_walk;
   }

   if (!found_drunkards_walk) {
      // calculate walk
      RNScalar *drunkards_walk = new RNScalar[nd->NBoundaries()];
      if (print_verbose) { printf("Createing drunkards walk boundary features...\n   "); fflush(stdout); }
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         if (print_verbose) RNProgressBar(ib, nd->NBoundaries());
         NeuronBoundary *boundary = nd->Boundary(ib);
         drunkards_walk[ib] = boundary->RandomWalk(TRUE);
      }
      if (print_verbose) printf("\ndone!\n");

      // open file for writing
      FILE *drunkards_fp = fopen(drunkards_walk_cache_filename, "wb");
      if (!drunkards_fp) { fprintf(stderr, "Failed to write to %s\n", drunkards_walk_cache_filename); return 0; }

      // write the number of boundaries
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, drunkards_fp);

      // write the drunkards walk to data
      fwrite(drunkards_walk, sizeof(RNScalar), nd->NBoundaries(), drunkards_fp);

      // close file
      fclose(drunkards_fp);

      // free memory
      delete[] drunkards_walk;
   }

   // return success
   return 1;
}



static int
CreatePredictionRandomWalkFeatures(const char root_filename[4096])
{
   // see if feature files already exist
   RNBoolean found_random_walk = FALSE;
   RNBoolean found_drunkards_walk = FALSE;

   // get random walk file
   char random_walk_cache_filename[4096];
   sprintf(random_walk_cache_filename, "%s/%s_random_walk.feature", postprocess_directory, root_filename);

   // see if the file exists to read
   FILE *random_walk_test_fp = fopen(random_walk_cache_filename, "rb");
   if (random_walk_test_fp) {
      found_random_walk = TRUE;
      fclose(random_walk_test_fp);
   }

   // get drunkards walk file
   char drunkards_walk_cache_filename[4096];
   sprintf(drunkards_walk_cache_filename, "%s/%s_drunkards_walk.feature", postprocess_directory, root_filename);

   // see if the file exists to read
   FILE *drunkards_walk_test_fp = fopen(drunkards_walk_cache_filename, "rb");
   if (drunkards_walk_test_fp) {
      found_drunkards_walk = TRUE;
      fclose(drunkards_walk_test_fp);
   }

   if (!found_random_walk) {
      // calculate walk
      RNScalar *random_walk = new RNScalar[nd->NPredictionBoundaries()];
      if (print_verbose) { printf("Createing random walk boundary features...\n   "); fflush(stdout); }
      for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
         if (print_verbose) RNProgressBar(ipb, nd->NPredictionBoundaries());
         NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);
         random_walk[ipb] = boundary->RandomWalk(FALSE);
      }
      if (print_verbose) printf("\ndone!\n");

      // open file for writing
      FILE *random_fp = fopen(random_walk_cache_filename, "wb");
      if (!random_fp) { fprintf(stderr, "Failed to write to %s\n", random_walk_cache_filename); return 0; }

      // write the number of boundaries
      int nboundaries = nd->NPredictionBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, random_fp);

      // write the random walk to data
      fwrite(random_walk, sizeof(RNScalar), nd->NPredictionBoundaries(), random_fp);

      // close files
      fclose(random_fp);

      // free memory
      delete[] random_walk;
   }

   if (!found_drunkards_walk) {
      // calculate walk
      RNScalar *drunkards_walk = new RNScalar[nd->NPredictionBoundaries()];
      if (print_verbose) { printf("Createing drunkards walk boundary features...\n   "); fflush(stdout); }
      for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
         if (print_verbose) RNProgressBar(ipb, nd->NPredictionBoundaries());
         NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);
         drunkards_walk[ipb] = boundary->RandomWalk(TRUE);
      }
      if (print_verbose) printf("\ndone!\n");

      // open file for writing
      FILE *drunkards_fp = fopen(drunkards_walk_cache_filename, "wb");
      if (!drunkards_fp) { fprintf(stderr, "Failed to write to %s\n", drunkards_walk_cache_filename); return 0; }

      // write the number of boundaries
      int nboundaries = nd->NPredictionBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, drunkards_fp);

      // write the drunkards walk to data
      fwrite(drunkards_walk, sizeof(RNScalar), nd->NPredictionBoundaries(), drunkards_fp);

      // close file
      fclose(drunkards_fp);

      // free memory
      delete[] drunkards_walk;
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
         else if (!strcmp(*argv, "-prediction")) prediction = 1;
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

   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create boundary features
   if (!CreateRandomWalkFeatures(root_filename)) exit(-1);

   if (prediction) {
      // read in prediction file
      char meta_root_filename[4096];
      sprintf(meta_root_filename, "output/ordermerge/%s", root_filename);

      R3Grid *prediction_grid = RNReadNeuronMetaRawFile(meta_root_filename);
      if (!prediction_grid) exit(-1);

      // create boundary predictions
      nd->CreatePredictions(prediction_grid);

      if (!CreatePredictionRandomWalkFeatures(root_filename)) exit(-1);

      // delete prediction grid
      delete prediction_grid;
   }

   // free up memory
   delete nd;

   // return success
   return 0;
}
