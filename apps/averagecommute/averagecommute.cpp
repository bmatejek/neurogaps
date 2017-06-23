// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static char *input_filename = NULL;
static NeuronData *nd = NULL;
static RNScalar maximum_affinity = 0.80;
static RNScalar minimum_affinity = 0.20;


// directory structure

static const char *distances_directory = "distances";



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
   if (!input_filename) { fprintf(stderr, "Need to supply input filename.\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

struct BoundaryRank {
   BoundaryRank(float score, int label) :
      score(score),
      label(label)
   {}

   float score;
   int label;
};


int BoundaryRankScore(BoundaryRank one, BoundaryRank two)
{
   return one.score < two.score;
}



static RNScalar ScaledAffinity(RNScalar affinity)
{
   if (affinity > maximum_affinity) return 1.0;
   else if (affinity < minimum_affinity) return 0.0;
   else return (affinity - minimum_affinity) / (maximum_affinity - minimum_affinity);
}



static NeuronCellular *NextStep(NeuronCellular *current)
{
   RNScalar affinity_sum = 0.0;
   for (int ib = 0; ib < current->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = current->Boundary(ib);

      // skip extracellulars
      if (boundary->OtherSupervoxel(current)->IsExtracellular()) continue;

      affinity_sum += ScaledAffinity(boundary->Mean());
   }

   RNScalar random_direction = affinity_sum * RNRandomScalar();
   affinity_sum = 0.0;
   for (int ib = 0; ib < current->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = current->Boundary(ib);

      // skip extracellulars
      if (boundary->OtherSupervoxel(current)->IsExtracellular()) continue;

      affinity_sum += ScaledAffinity(boundary->Mean());
      if (affinity_sum > random_direction) return (NeuronCellular *)boundary->OtherSupervoxel(current);
   }

   // should never reach here
   rn_assertion(FALSE);
   return NULL;
}




static RNScalar Entropy(int npositives, int nnegatives)
{
   RNScalar positive_proportion = npositives / (RNScalar)(npositives + nnegatives);
   RNScalar negative_proportion = nnegatives / (RNScalar)(npositives + nnegatives);

   // return the entropy
   return -1 * positive_proportion * log2(positive_proportion) - negative_proportion * log2(negative_proportion);
}





int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // start at a random cellular
   int ncelluars = nd->NCellulars();

   // keep track of first_pass sums
   unsigned int **first_pass_sums = new unsigned int *[nd->NCellulars()];
   int **noccurrences = new int *[nd->NCellulars()];
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      first_pass_sums[ic1] = new unsigned int[nd->NCellulars()];
      noccurrences[ic1] = new int[nd->NCellulars()];
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         first_pass_sums[ic1][ic2] = 0;
         noccurrences[ic1][ic2] = 0;
      }
   }

   // start statistics
   RNTime first_pass_time;
   first_pass_time.Read();

   int niterations = 10;
   for (int it = 0; it < niterations; ++it) {
      RNProgressBar(it, niterations);
      int start_index = (int)(RNRandomScalar() * ncelluars);
      // get the current cellular
      NeuronCellular *current_cellular = nd->Cellular(start_index);

      // keep track of reset counter
      RNBoolean **reset = new RNBoolean *[nd->NCellulars()];
      for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
         reset[ic1] = new RNBoolean[nd->NCellulars()];
         for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
            reset[ic1][ic2] = FALSE;
         }
      }


      // keep track of the last occurrence for every cellular
      int *last_pass = new int[nd->NCellulars()];
      for (int ic = 0; ic < nd->NCellulars(); ++ic)
         last_pass[ic] = -1;

      // do random walk on supervoxels
      int nsteps = 10000;
      for (int is = 0; is < nsteps; ++is) {
         last_pass[current_cellular->DataIndex()] = is;

         // update first pass sums
         for (int ic = 0; ic < nd->NCellulars(); ++ic) {
            if (last_pass[ic] == -1) continue;
            int elapsed_time = is - last_pass[ic];

            if (reset[current_cellular->DataIndex()][ic]) {
               first_pass_sums[current_cellular->DataIndex()][ic] += elapsed_time;
               noccurrences[current_cellular->DataIndex()][ic]++;
               reset[current_cellular->DataIndex()][ic] = FALSE;
            }
         }

         for (int ic = 0; ic < nd->NCellulars(); ++ic) {
            reset[ic][current_cellular->DataIndex()] = TRUE;
         }

         // update the current supervoxel
         current_cellular = NextStep(current_cellular);
      }

      // free memory
      for (int ic = 0; ic < nd->NCellulars(); ++ic)
         delete[] reset[ic];
      delete[] reset;
      delete[] last_pass;
   }
   printf("\ndone in %0.2f seconds\n", first_pass_time.Elapsed());

   RNScalar **average_first_pass = new RNScalar *[nd->NCellulars()];
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      average_first_pass[ic1] = new RNScalar[nd->NCellulars()];
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         average_first_pass[ic1][ic2] = first_pass_sums[ic1][ic2] / (RNScalar)noccurrences[ic1][ic2];
      }
   }

   std::vector<BoundaryRank> boundaries = std::vector<BoundaryRank>();

   // go through all boundaries and print results
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get cellular indices
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      // skip extracellulars
      if (supervoxel_one->IsExtracellular()) continue;
      if (supervoxel_two->IsExtracellular()) continue;

      int label = (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel());

      boundaries.push_back(BoundaryRank(average_first_pass[supervoxel_one->DataIndex()][supervoxel_two->DataIndex()], label));
   }

   sort(boundaries.begin(), boundaries.end(), BoundaryRankScore);

   // get number of positive and negative results
   int nparent_positives = 0;
   int nparent_negatives = 0;
   for (unsigned int ie = 0; ie < boundaries.size(); ++ie) {
      if (boundaries[ie].label) nparent_positives++;
      else nparent_negatives++;
   }

   // find the number of right and left ones and zeros
   int *left_negatives = new int[boundaries.size()];
   int *right_positives = new int[boundaries.size()];
   for (int ie = 0; ie < (int)boundaries.size(); ++ie) {
      left_negatives[ie] = 0;
      right_positives[ie] = 0;
   }

   // go from from left to right to get left negatives array
   for (int ie = 1; ie < (int)boundaries.size(); ++ie) {
      if (boundaries[ie - 1].label) left_negatives[ie] = left_negatives[ie - 1];
      else left_negatives[ie] = left_negatives[ie - 1] + 1;
   }

   // go from right to left to get right positives array
   for (int ie = (int)boundaries.size() - 2; ie >= 0; --ie) {
      if (boundaries[ie + 1].label) right_positives[ie] = right_positives[ie + 1] + 1;
      else right_positives[ie] = right_positives[ie + 1];
   }

   RNScalar best_information_gain = 0.0;
   RNScalar best_value = 0.0;

   // get parent entropy
   RNScalar parent_entropy = Entropy(nparent_positives, nparent_negatives);


   // calculate the entropy at each location
   for (int ie = 0; ie < (int)boundaries.size(); ++ie) {
      int right_positive = right_positives[ie];
      int left_negative = left_negatives[ie];
      int left_positive = (ie - left_negative);
      int right_negative = ((boundaries.size() - 1) - ie - right_positive);

      // find current label
      int current_label = boundaries[ie].label;

      // give benefit of doubt for splitting here
      if (current_label) {
         if (left_positive > right_positive) left_positive++;
         else right_positive++;
      }
      else {
         if (left_negative > right_negative) left_negative++;
         else right_negative++;
      }

      // calculate the left and right entropy
      RNScalar left_entropy = Entropy(left_positive, left_negative);
      RNScalar right_entropy = Entropy(right_positive, right_negative);

      RNScalar left_proportion = (left_positive + left_negative) / (RNScalar)(boundaries.size());
      RNScalar right_proportion = (right_positive + right_negative) / (RNScalar)(boundaries.size());

      RNScalar average_child_entropy = left_proportion * left_entropy + right_proportion * right_entropy;
      RNScalar information_gain = parent_entropy - average_child_entropy;

      if (information_gain > best_information_gain) {
         best_information_gain = information_gain;

         RNScalar split_right = (right_negative + left_positive) / (RNScalar)boundaries.size();
         RNScalar split_left = (left_negative + right_positive) / (RNScalar)boundaries.size();

         if (split_right > split_left)
            best_value = split_right;
         else
            best_value = split_left;
      }

   }
   
   printf("%lf\n", best_information_gain);
   printf("%lf\n", best_value);

   // return success
   return 0;
}
