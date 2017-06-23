// Source file for the simple file conversion algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *train_filename = NULL;
static char *test_filename = NULL;



// global variables

R3Grid **train_grid = NULL;
R3Grid **test_grid = NULL;



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
         if (!train_filename) train_filename = *argv;
         else if (!test_filename) test_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // verify input
   if (!train_filename || !test_filename) { fprintf(stderr, "Need to supply both training and testing filename\n"); return 0; }

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

   // read in training and testing grids
   if (print_verbose) { printf("Reading training affinities..."); fflush(stdout); }
   RNTime train_time;
   train_time.Read();
   train_grid = RNReadNeuronMetaRawFile(train_filename, TRUE);
   if (print_verbose) { printf("done in %0.2f seconds\n", train_time.Elapsed()); }

   std::vector<float> training_affinities = std::vector<float>();

   // iterate through all training affinities
   for (int dim = 0; dim <= 2; ++dim) {
      rn_assertion(train_grid[dim] != NULL);
      for (int iv = 0; iv < train_grid[dim]->NEntries(); ++iv) {
         training_affinities.push_back(train_grid[dim]->GridValue(iv));
      }
   }

   // delete training affinities
   for (int dim = 0; dim <= 2; ++dim) {
      delete train_grid[dim];
   }
   delete[] train_grid;

   // sort the affinities
   if (print_verbose) { printf("Sorting training affinities..."); fflush(stdout); }
   RNTime train_sort_time;
   train_sort_time.Read();
   sort(training_affinities.begin(), training_affinities.end());
   if (print_verbose) { printf("done in %0.2f seconds\n", train_sort_time.Elapsed()); }

   if (print_verbose){ printf("Reading testing affinities..."); fflush(stdout); }
   RNTime test_time;
   test_time.Read();
   test_grid = RNReadNeuronMetaRawFile(test_filename, TRUE);
   if (print_verbose) { printf("done in %0.2f seconds\n", test_time.Elapsed()); }

   std::vector<float> testing_affinities = std::vector<float>();

   // iteratr through all testing affinities
   for (int dim = 0; dim <= 2; ++dim) {
      rn_assertion(test_grid[dim] != NULL);
      for (int iv = 0; iv < test_grid[dim]->NEntries(); ++iv) {
         testing_affinities.push_back(test_grid[dim]->GridValue(iv));
      }
   }

   // sort the affinities
   if (print_verbose) { printf("Sorting testing affinities..."); fflush(stdout); }
   RNTime test_sort_time;
   test_sort_time.Read();
   sort(testing_affinities.begin(), testing_affinities.end());
   if (print_verbose) { printf("done in %0.2f seconds\n", test_sort_time.Elapsed()); }

   // get mean and stddev of both affinities
   RNScalar train_mean = 0.0;
   RNScalar train_stddev = 0.0;
   RNScalar test_mean = 0.0;
   RNScalar test_stddev = 0.0;

   // get means
   for (unsigned int ie = 0; ie < training_affinities.size(); ++ie) {
      train_mean += training_affinities[ie] / training_affinities.size();
   }
   for (unsigned int ie = 0; ie < testing_affinities.size(); ++ie) {
      test_mean += testing_affinities[ie] / testing_affinities.size();
   }
   printf("Means:\n  Train: %0.6lf\n  Test:  %0.6lf\n", train_mean, test_mean);
   for (unsigned int ie = 0; ie < training_affinities.size(); ++ie) {
      train_stddev += (training_affinities[ie] - train_mean) * (training_affinities[ie] - train_mean);
   }
   for (unsigned int ie = 0; ie < testing_affinities.size(); ++ie) {
      test_stddev += (testing_affinities[ie] - test_mean) * (testing_affinities[ie] - test_mean);
   }
   train_stddev = sqrt(train_stddev / (training_affinities.size() - 1));
   test_stddev = sqrt(test_stddev / (testing_affinities.size() - 1));
   printf("StdDevs:\n  Train: %0.6lf\n  Test:  %0.6lf\n", train_stddev, test_stddev);

   for (int dim = 0; dim <= 2; ++dim) {
      for (int iv = 0; iv < test_grid[dim]->NEntries(); ++iv) {
         float test_affinity = test_grid[dim]->GridValue(iv);
         float new_affinity = (((test_affinity - test_mean) / test_stddev) * train_stddev + train_mean);
         test_grid[dim]->SetGridValue(iv, new_affinity);
      }
   }

   char *extp = strrchr(test_filename, '.');
   *extp = '\0';

   // write new meta/raw file
   char output_filename[4096];
   sprintf(output_filename, "%s_normalized", test_filename);

   RNMeta meta = RNMeta("Float32", test_grid[0]->XResolution(), test_grid[0]->YResolution(), test_grid[0]->ZResolution(), 3);
   if (!RNWriteNeuronMetaRawFile(output_filename, meta, test_grid)) exit(-1);

   // delete testing affinities
   for (int dim = 0; dim <= 2; ++dim) {
      delete test_grid[dim];
   }
   delete[] test_grid;

   // return success
   return 0;
}
