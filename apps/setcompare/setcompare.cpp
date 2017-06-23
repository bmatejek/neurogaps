// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int only_true_matches = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static unsigned int **set_agreements = NULL;
static unsigned int **set_disagreements = NULL;



// directory structure
static const char *algs_directory = "algs_data/dijkstra";
static const char *results_directory = "results/dijkstrasets";



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



static int ReadSets(const char *root_filename)
{
   // create arrays
   set_agreements = new unsigned int *[nd->NCellulars()];
   set_disagreements = new unsigned int *[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      set_agreements[ic] = new unsigned int[nd->NCellulars()];
      set_disagreements[ic] = new unsigned int[nd->NCellulars()];
   }

   // get filenames
   char agree_filename[4096];
   if (only_true_matches) sprintf(agree_filename, "%s/%s_truth.agree", algs_directory, root_filename);
   else sprintf(agree_filename, "%s/%s.agree", algs_directory, root_filename);
   char disagree_filename[4096];
   if (only_true_matches) sprintf(disagree_filename, "%s/%s_truth.disagree", algs_directory, root_filename);
   else sprintf(disagree_filename, "%s/%s.disagree", algs_directory, root_filename);

   // open files
   FILE *agree_fp = fopen(agree_filename, "rb");
   if (!agree_fp) { fprintf(stderr, "Failed to read %s\n", agree_filename); return 0; }
   FILE *disagree_fp = fopen(disagree_filename, "rb");
   if (!disagree_fp) { fprintf(stderr, "Failed to read %s\n", disagree_filename); fclose(agree_fp); return 0; }

   int nagree_cellulars;
   int ndisagree_cellulars;
   fread(&nagree_cellulars, sizeof(int), 1, agree_fp);
   fread(&ndisagree_cellulars, sizeof(int), 1, disagree_fp);
   rn_assertion(nagree_cellulars == nd->NCellulars());
   rn_assertion(ndisagree_cellulars == nd->NCellulars());

   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      fread(set_agreements[ic], sizeof(unsigned int), nd->NCellulars(), agree_fp);
      fread(set_disagreements[ic], sizeof(unsigned int), nd->NCellulars(), disagree_fp);
   }

   // close file
   fclose(agree_fp);
   fclose(disagree_fp);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Precision and recall functions
////////////////////////////////////////////////////////////////////////


struct SetPair {
   RNScalar agreement;
   RNBoolean match;
};



int SetPairCompare(SetPair a, SetPair b)
{
   return a.agreement > b.agreement;
}



static int CreateAllPrecisionRecall(const char *root_filename){
   std::vector<SetPair> all_sets = std::vector<SetPair>();

   int nall_matches = 0;
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      NeuronCellular *cellular_one = nd->Cellular(ic1);
      for (int ic2 = ic1 + 1; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);

         // do these cellulars belong to the same voxel
         RNBoolean match = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());
         if (match) nall_matches++;

         unsigned int nagreements = set_agreements[ic1][ic2];
         unsigned int ndisagreements = set_disagreements[ic1][ic2] + set_disagreements[ic2][ic1];
         if (nagreements + ndisagreements == 0) continue;

         RNScalar agreement_ratio = (RNScalar)nagreements / ((RNScalar)nagreements + ndisagreements);
         rn_assertion(nagreements == set_agreements[ic2][ic1]);

         SetPair pair;
         pair.agreement = agreement_ratio;
         pair.match = match;
         all_sets.push_back(pair);
      }
   }

   std::sort(all_sets.begin(), all_sets.end(), SetPairCompare);
   
   char all_filename[4096];
   if (only_true_matches) sprintf(all_filename, "%s/%s_truth_all_pairs.pandr", results_directory, root_filename);
   else sprintf(all_filename, "%s/%s_all_pairs.pandr", results_directory, root_filename);

   FILE *all_fp = fopen(all_filename, "w");
   if (!all_fp) { fprintf(stderr, "Failed to read file %s\n", all_filename); return 0; }

   // create precision and recall curves
   int ncorrect = 0;
   int nseen = 0;
   for (unsigned int ie = 0; ie < all_sets.size(); ++ie) {
      if (all_sets[ie].match) {
         nseen++;
         ncorrect++;
      }
      if (ie % 100 == 0) {
         RNScalar recall = (RNScalar)nseen / (RNScalar)nall_matches;
         RNScalar precision = (RNScalar)ncorrect / (ie + 1);
         fprintf(all_fp, "%lf,%lf\n", recall, precision);
      }
   }

   // close file   
   fclose(all_fp);

   // return success
   return 1;
}



static int CreateNeighborPrecisionRecall(const char *root_filename)
{
   std::vector<SetPair> neighbor_sets = std::vector<SetPair>();

   int nmatches = 0;
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_one->IsExtracellular()) continue;
      if (supervoxel_two->IsExtracellular()) continue;
      NeuronCellular *cellular_one = (NeuronCellular *)supervoxel_one;
      NeuronCellular *cellular_two = (NeuronCellular *)supervoxel_two;

      int ic1 = cellular_one->DataIndex();
      int ic2 = cellular_two->DataIndex();

      // do these cellulars belong to the same voxel
      RNBoolean match = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());
      if (match) nmatches++;

      unsigned int nagreements = set_agreements[ic1][ic2];
      unsigned int ndisagreements = set_disagreements[ic1][ic2] + set_disagreements[ic2][ic1];
      if (ndisagreements + nagreements == 0) { continue; }

      RNScalar agreement_ratio = (RNScalar)nagreements / ((RNScalar)nagreements + ndisagreements);
      rn_assertion(nagreements == set_agreements[ic2][ic1]);

      SetPair pair;
      pair.agreement = agreement_ratio;
      pair.match = match;
      neighbor_sets.push_back(pair);
   }

   std::sort(neighbor_sets.begin(), neighbor_sets.end(), SetPairCompare);

   char neighbor_filename[4096];
   if (only_true_matches) sprintf(neighbor_filename, "%s/%s_truth_neighbors.pandr", results_directory, root_filename);
   else sprintf(neighbor_filename, "%s/%s_neighbors.pandr", results_directory, root_filename);

   FILE *neighbor_fp = fopen(neighbor_filename, "w");
   if (!neighbor_fp) { fprintf(stderr, "Failed to read file %s\n", neighbor_filename); return 0; }

   // create precision and recall curves
   int ncorrect = 0;
   int nseen = 0;
   for (unsigned int ie = 0; ie < neighbor_sets.size(); ++ie) {
      if (neighbor_sets[ie].match) {
         nseen++;  ncorrect++;
      }
      RNScalar recall = (RNScalar)nseen / (RNScalar)nmatches;
      RNScalar precision = (RNScalar)ncorrect / (ie + 1);
      fprintf(neighbor_fp, "%lf,%lf\n", recall, precision);
   }

   // close file   
   fclose(neighbor_fp);

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
         else if (!strcmp(*argv, "-truth")) { only_true_matches = 1; }
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

   // read input neuron file
   if (!ReadData(input_filename)) exit(-1);

   // get the root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in sets
   if (!ReadSets(root_filename)) exit(-1);
   if (!CreateAllPrecisionRecall(root_filename)) exit(-1);
   if (!CreateNeighborPrecisionRecall(root_filename)) exit(-1);

   // free up memory
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      delete[] set_agreements[ic];
      delete[] set_disagreements[ic];
   }
   delete[] set_agreements;
   delete[] set_disagreements;
   delete nd;

   // return success
   return 0;
}
