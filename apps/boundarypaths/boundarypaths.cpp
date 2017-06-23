// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>


// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static int cut_off = 1;



// global variables

static NeuronData *nd;
static float **data_terms = NULL;



// directory structure

static const char *algs_directory = "algs_data/graphcut";
static const char *tmp_algs_directory = "algs_data/features/tmp";
static const char *hierarchical_directory = "algs_data/hierarchical";
static const char *results_directory = "results/boundaries";
static const char *visuals_directory = "visuals/boundaries";
static const char *cache_directory = "cache";



////////////////////////////////////////////////////////////////////////
// Useful structs
////////////////////////////////////////////////////////////////////////

struct BPBoundary
{
   BPBoundary(NeuronBoundary *boundary, float raw_agreement_score, RNBoolean label) :
      boundary(boundary),
      raw_agreement_score(raw_agreement_score),
      label(label)
   {
   }

   NeuronBoundary *boundary;
   float raw_agreement_score;
   RNBoolean label;
};



// separate ranking functions
int BPBoundaryRawRank(BPBoundary one, BPBoundary two)
{
   return one.raw_agreement_score > two.raw_agreement_score;
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int
ReadData(const char *filename)
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



static int ReadUnaryFile(const char unary_filename[4096])
{
   // read in unary file
   FILE *unary_fp = fopen(unary_filename, "rb");
   if (!unary_fp) { fprintf(stderr, "Failed to read unary filename: %s\n", unary_filename); return 0; }

   // create array of data terms
   int unary_ncellulars;
   fread(&unary_ncellulars, sizeof(int), 1, unary_fp);
   rn_assertion(unary_ncellulars == nd->NCellulars());

   data_terms = new float *[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      data_terms[ic] = new float[nd->NCellulars()];
   }

   for (int i = 0; i < nd->NCellulars(); ++i) {
      fread(data_terms[i], sizeof(float), unary_ncellulars, unary_fp);
   }

   // close file
   fclose(unary_fp);

   // return success
   return 1;
}



static int
MergeResults(char root_filename[4096])
{
   // allocate memory for counter arrays
   float **cellular_matches = new float *[nd->NCellulars()];
   float *cellular_counter = new float[nd->NCellulars()];
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      cellular_matches[ic1] = new float[nd->NCellulars()];
      cellular_counter[ic1] = 0.0;
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         cellular_matches[ic1][ic2] = 0.0;
      }
   }

   // get cache filename
   char cache_filename[4096];
   sprintf(cache_filename, "%s/%s_boundary_paths_%02d.cache", cache_directory, root_filename, cut_off);

   // if cache does not exist, create it
   FILE *cache_fp = fopen(cache_filename, "rb");
   if (!cache_fp) {
      // get unary filename
      char unary_filename[4096];
      sprintf(unary_filename, "%s/%s_random_forest.unary", algs_directory, root_filename);

      // read in data term
      if (!ReadUnaryFile(unary_filename)) return 0;

      // start statistics
      RNTime merge_time;
      merge_time.Read();

      // go through all cellulars
      if (print_verbose) { printf("Reading all cellulars...\n  "); fflush(stdout); }
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (print_verbose) RNProgressBar(ic, nd->NCellulars());
         NeuronCellular *cellular = nd->Cellular(ic);

         // get the input filename
         char input_filename[4096];
         sprintf(input_filename, "%s/%s_%06d.path", tmp_algs_directory, root_filename, ic);

         // open file
         FILE *fp = fopen(input_filename, "rb");
         if (!fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

         // read the number of elements for this cellular
         int npaths;
         fread(&npaths, sizeof(int), 1, fp);
         
	 if (!cellular->IsOnBoundary()) continue;
	 
         // go through all boundary supervoxels
         for (int ip = 0; ip < npaths; ++ip) {
            // read the start index, end index, and number of elements on the path
            int start_index;
            int end_index;
            fread(&start_index, sizeof(int), 1, fp);
            fread(&end_index, sizeof(int), 1, fp);

            if (ic != start_index)  {
               fprintf(stderr, "%s\n", input_filename);
               break;
            }
            rn_assertion((0 <= end_index) && (end_index < nd->NCellulars()));

            int ncellulars_on_path;
            fread(&ncellulars_on_path, sizeof(int), 1, fp);
            rn_assertion(ncellulars_on_path > 0);

            std::vector<int> cellulars_on_path = std::vector<int>();
            for (int icp = 0; icp < ncellulars_on_path; ++icp) {
               int cellular_index;
               fread(&cellular_index, sizeof(int), 1, fp);
               cellulars_on_path.push_back(cellular_index);
            }

	    if (!nd->Cellular(end_index)->IsOnBoundary()) continue;

            // get the data term between these boundary terms
            float data_term = data_terms[start_index][end_index];
            if (data_term > cut_off / 10.0) continue;

            // increment the individual counters
            for (int icp = 0; icp < ncellulars_on_path; ++icp) {
               int cellular_index = cellulars_on_path[icp];
               // divide by two since every path occurs twice
               cellular_counter[cellular_index] += (1.0 - data_term) / 2.0;
            }

            // increment the agreement for all pairs of cellulars on this path
            for (int ic1 = 0; ic1 < ncellulars_on_path; ++ic1) {
               int cellular_one_data_index = cellulars_on_path[ic1];
               for (int ic2 = 0; ic2 < ncellulars_on_path; ++ic2) {
                  int cellular_two_data_index = cellulars_on_path[ic2];

                  // divide by two since every path occurs twice
                  cellular_matches[cellular_one_data_index][cellular_two_data_index] += (1.0 - data_term) / 2.0;
               }
            }
         }

         // close file
         fclose(fp);
      }
      if (print_verbose) { printf("\ndone in %0.2f seconds!\n", merge_time.Elapsed()); }

      FILE *write_cache_fp = fopen(cache_filename, "wb");
      if (!write_cache_fp) { fprintf(stderr, "Failed to write cache file %s\n", cache_filename); return 0; }

      // save cellular_counter and cellular_matches
      fwrite(cellular_counter, sizeof(float), nd->NCellulars(), write_cache_fp);
      for (int ic = 0; ic < nd->NCellulars(); ++ic)
         fwrite(cellular_matches[ic], sizeof(float), nd->NCellulars(), write_cache_fp);

      // close file
      fclose(write_cache_fp);
   }
   else {
      if (fread(cellular_counter, sizeof(float), nd->NCellulars(), cache_fp) != (unsigned int)nd->NCellulars()) {
         fprintf(stderr, "Failed to read cache file %s\n", cache_filename);
         return 0;
      }
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (fread(cellular_matches[ic], sizeof(float), nd->NCellulars(), cache_fp) != (unsigned int)nd->NCellulars()) {
            fprintf(stderr, "Failed to read cache file %s\n", cache_filename);
            return 0;
         }
      }

      // close file
      fclose(cache_fp);
   }

   std::vector<BPBoundary> boundary_pairs = std::vector<BPBoundary>();
   // go through all boundaries
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get the supervoxels
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      // skip extracellulars
      if (supervoxel_one->IsExtracellular() || supervoxel_two->IsExtracellular()) continue;

      // get data indices
      int cellular_one_data_index = supervoxel_one->DataIndex();
      int cellular_two_data_index = supervoxel_two->DataIndex();

      // find the counter and match score
      float match_score = cellular_matches[cellular_one_data_index][cellular_two_data_index];

      // create boundary pairs
      BPBoundary boundary_pair = BPBoundary(boundary, match_score, supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel());

      // add to the vector
      boundary_pairs.push_back(boundary_pair);
   }

   // create PyImage
   RNPyImage *raw_image = new RNPyImage();
   raw_image->SetOutputDirectory(visuals_directory);

   // sort the boundary pairs according to the raw agreement
   std::sort(boundary_pairs.begin(), boundary_pairs.end(), BPBoundaryRawRank);

   // find the number of elements correct to both the left and right of the current position
   int *right_zeros = new int[boundary_pairs.size()];
   int *left_ones = new int[boundary_pairs.size()];
   // intialize to zero
   for (int ib = 0; ib < (int)boundary_pairs.size(); ++ib) {
      right_zeros[ib] = 0;
      left_ones[ib] = 0;
   }
   // go larger to small
   for (int ib = 1; ib < (int)boundary_pairs.size(); ++ib) {
      if (boundary_pairs[ib - 1].label) left_ones[ib] = left_ones[ib - 1] + 1;
      else left_ones[ib] = left_ones[ib - 1];
   }
   // go smaller to larger
   for (int ib = (int)boundary_pairs.size() - 2; ib >= 0; --ib) {
      if (boundary_pairs[ib + 1].label) right_zeros[ib] = right_zeros[ib + 1];
      else right_zeros[ib] = right_zeros[ib + 1] + 1;
   }

   // useful variables
   char title[4096];

   // create plots
   RNPyPointPlot *raw_pandr = new RNPyPointPlot();
   RNPyPointPlot *raw_xscore = new RNPyPointPlot();
   RNPyPointPlot *raw_proportion = new RNPyPointPlot();

   // set x and y labels
   raw_pandr->SetXLabel("Recall");
   raw_pandr->SetYLabel("Precision");
   raw_xscore->SetXLabel("Similarity Score");
   raw_xscore->SetYLabel("Proportion Greater with Positive Label");
   raw_proportion->SetXLabel("Similarity Score");
   raw_proportion->SetYLabel("Proportion Correct at Cut Off");

   // set titles
   sprintf(title, "%s Score Precision and Recall - Cut Off %02d", root_filename, cut_off);
   raw_pandr->SetTitle(title);
   sprintf(title, "%s Score versus Proportion Positive - Cut Off %02d", root_filename, cut_off);
   raw_xscore->SetTitle(title);
   sprintf(title, "%s Score with Threshold - Cut Off %02d", root_filename, cut_off);
   raw_proportion->SetTitle(title);

   // set maxima and legend
   raw_proportion->SetLegend("Proportion Correct");
   raw_proportion->SetExtremaType(MAX_EXTREMA);

   int nraw_correct = 0;
   for (int ib = 0; ib < (int)boundary_pairs.size(); ++ib) {
      int nseen = ib + 1;
      if (boundary_pairs[ib].label) nraw_correct++;

      RNScalar precision = nraw_correct / (RNScalar)nseen;
      RNScalar recall = nseen / (RNScalar)boundary_pairs.size();

      // add to graphs
      raw_pandr->InsertPoint(R2Point(recall, precision));
      raw_xscore->InsertPoint(R2Point(boundary_pairs[ib].raw_agreement_score, precision));

      // add to proportion graph
      int ncorrect = right_zeros[ib] + left_ones[ib];
      if (boundary_pairs[ib].label) ncorrect += 1;
      raw_proportion->InsertPoint(R2Point(boundary_pairs[ib].raw_agreement_score, ncorrect / (RNScalar)boundary_pairs.size()));
   }

   // add plots to images and save
   raw_image->InsertPyPlot(raw_pandr);
   raw_image->InsertPyPlot(raw_xscore);
   raw_image->InsertPyPlot(raw_proportion);

   char raw_output_filename[4096];
   sprintf(raw_output_filename, "%s/%s_raw_%02d.pyimage", results_directory, root_filename, cut_off);

   raw_image->WriteImageFile(raw_output_filename);

   // save the hierarchical merge order
   char hierarchical_filename[4096];
   sprintf(hierarchical_filename, "%s/%s_boundary_paths.hier", hierarchical_directory, root_filename);

   // open file
   FILE *hier_fp = fopen(hierarchical_filename, "wb");
   if (!hier_fp) { fprintf(stderr, "Failed to write hierarchical merge file %s\n", hierarchical_filename); return 0; }

   // write the merge order 
   int nboundaries = boundary_pairs.size();
   fwrite(&nboundaries, sizeof(int), 1, hier_fp);
   for (int ib = 0; ib < nboundaries; ++ib) {
      NeuronBoundary *boundary = boundary_pairs[ib].boundary;
      int boundary_data_index = boundary->DataIndex();
      fwrite(&boundary_data_index, sizeof(int), 1, hier_fp);
   }

   // close file
   fclose(hier_fp);

   // free memory
   delete[] left_ones;
   delete[] right_zeros;
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      delete[] cellular_matches[ic];
   }
   delete[] cellular_matches;
   delete[] cellular_counter;
   delete raw_pandr;
   delete raw_xscore;
   delete raw_proportion;
   delete raw_image;

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
	 else if (!strcmp(*argv, "-cutoff")) { argv++; argc--; cut_off = atoi(*argv); }
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

   if (cut_off < 1 || cut_off > 10) {
      fprintf(stderr, "Cut off input must be between 1 and 10, inclusive\n");
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

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';
   
   // merge all of the results
   if (!MergeResults(root_filename)) exit(-1);
   
   // free memory
   delete nd;

   // return success
   return 0;
}
