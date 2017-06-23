// Source file for the subsection algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int min_merge = 0;



// global variables

static char *root_filename = NULL;
static const char *output_directory = "output";
static char *prediction_method = NULL;
static R3Grid *truth_grid = NULL;



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
         else if (!strcmp(*argv, "-min_merge")) { argv++; argc--; min_merge = atoi(*argv); }
         // number of files to generate
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!root_filename) root_filename = *argv;
         else if (!prediction_method) prediction_method = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!root_filename) { fprintf(stderr, "Need to supply root filename.\n");  return 0; }
   if (!prediction_method) { fprintf(stderr, "Need to supply prediction method.\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

static R3Grid *
HorizontalMerge(R3Grid *left, R3Grid *right)
{
   // get the resolutions of merged grid
   rn_assertion(left->YResolution() == right->YResolution());
   int yresolution = left->YResolution();
   rn_assertion(left->ZResolution() == right->ZResolution());
   int zresolution = left->ZResolution();
   int xresolution = left->XResolution() + right->XResolution();

   // form the merged grid
   R3Grid *merged = new R3Grid(xresolution, yresolution, zresolution);
   if (!merged) { fprintf(stderr, "Failed to allocate memory for merged grid.\n"); return NULL; }

   // get maximum values /* add one since grids are zero indexed */
   int left_max = (int)(left->Maximum() + 1.5);
   int right_max = (int)(right->Maximum() + 1.5);

   // create counter for boundary values
   int **boundary_counter = new int *[right_max];
   for (int i = 0; i < right_max; ++i) {
      boundary_counter[i] = new int[left_max];
      for (int j = 0; j < left_max; ++j) {
         boundary_counter[i][j] = 0;
      }
   }

   int leftx = left->XResolution() - 1;
   int rightx = 0;

   for (int iy = 0; iy < yresolution; ++iy) {
      for (int iz = 0; iz < zresolution; ++iz) {
         // get the bottom and top values
         int left_value = (int)(left->GridValue(leftx, iy, iz) + 0.5);
         int right_value = (int)(right->GridValue(rightx, iy, iz) + 0.5);

         boundary_counter[right_value][left_value]++;
      }
   }

   // keep track of the number of unique neurons in top
   int nunstitched = 0;
   int *mapping = new int[right_max];
   for (int i = 0; i < right_max; ++i) {
      int ntotal_voxels = 0;
      // find best way to stitch top to bottom
      int best_match = -1;
      int best_match_score = min_merge;
      for (int j = 0; j < left_max; ++j) {
         if (boundary_counter[i][j] > best_match_score) {
            best_match = j;
            best_match_score = boundary_counter[i][j];
         }
         ntotal_voxels += boundary_counter[i][j];
      }

      // set the mapping value
      if (best_match != -1) {
         mapping[i] = best_match;
      }
      else {
         mapping[i] = right_max + nunstitched;
         nunstitched++;
      }
   }

   // set values for merged
   for (int ix = 0; ix < xresolution; ++ix) {
      for (int iy = 0; iy < yresolution; ++iy) {
         for (int iz = 0; iz < zresolution; ++iz) {
            if (ix < left->XResolution()) {
               merged->SetGridValue(ix, iy, iz, left->GridValue(ix, iy, iz));
            }
            else {
               int right_value = (int)(right->GridValue(ix - left->XResolution(), iy, iz) + 0.5);
               merged->SetGridValue(ix, iy, iz, mapping[right_value]);
            }
         }
      }
   }

   // free memory
   for (int i = 0; i < right_max; ++i)
      delete[] boundary_counter[i];
   delete[] boundary_counter;
   delete[] mapping;

   // return success
   return merged;
}



static R3Grid *
VerticalMerge(R3Grid *bottom, R3Grid *top)
{
   // get the resolutions of merged grid
   rn_assertion(bottom->XResolution() == top->XResolution());
   int xresolution = bottom->XResolution();
   rn_assertion(bottom->ZResolution() == top->ZResolution());
   int zresolution = bottom->ZResolution();
   int yresolution = bottom->YResolution() + top->YResolution();

   // form the merged grid
   R3Grid *merged = new R3Grid(xresolution, yresolution, zresolution);
   if (!merged) { fprintf(stderr, "Failed to allocate memory for merged grid.\n"); return NULL; }

   // get maximum values /* add one since grids are zero indexed */
   int bottom_max = (int)(bottom->Maximum() + 1.5);
   int top_max = (int)(top->Maximum() + 1.5);

   // create counter for boundary values
   int **boundary_counter = new int *[top_max];
   for (int i = 0; i < top_max; ++i) {
      boundary_counter[i] = new int[bottom_max];
      for (int j = 0; j < bottom_max; ++j) {
         boundary_counter[i][j] = 0;
      }
   }

   int bottomy = bottom->YResolution() - 1;
   int topy = 0;

   for (int ix = 0; ix < xresolution; ++ix) {
      for (int iz = 0; iz < zresolution; ++iz) {
         // get the bottom and top values
         int bottom_value = (int)(bottom->GridValue(ix, bottomy, iz) + 0.5);
         int top_value = (int)(top->GridValue(ix, topy, iz) + 0.5);

         boundary_counter[top_value][bottom_value]++;
      }
   }

   // keep track of the number of unique neurons in top
   int nunstitched = 0;
   int *mapping = new int[top_max];
   for (int i = 0; i < top_max; ++i) {
      // find best way to stitch top to bottom
      int best_match = -1;
      int best_match_score = min_merge;
      int ntotal_voxels = 0;
      for (int j = 0; j < bottom_max; ++j) {
         if (boundary_counter[i][j] > best_match_score) {
            best_match = j;
            best_match_score = boundary_counter[i][j];
         }
         ntotal_voxels += boundary_counter[i][j];
      }

      // set the mapping value
      if (best_match != -1) {
         mapping[i] = best_match;
      }
      else {
         mapping[i] = bottom_max + nunstitched;
         nunstitched++;
      }
   }

   // set values for merged
   for (int ix = 0; ix < xresolution; ++ix) {
      for (int iy = 0; iy < yresolution; ++iy) {
         for (int iz = 0; iz < zresolution; ++iz) {
            if (iy < bottom->YResolution()) {
               merged->SetGridValue(ix, iy, iz, bottom->GridValue(ix, iy, iz));
            }
            else {
               int top_value = (int)(top->GridValue(ix, iy - bottom->YResolution(), iz) + 0.5);
               merged->SetGridValue(ix, iy, iz, mapping[top_value]);
            }
         }
      }
   }

   // free memory
   for (int i = 0; i < top_max; ++i)
      delete[] boundary_counter[i];
   delete[] boundary_counter;
   delete[] mapping;

   // return success
   return merged;
}




int
main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // get all four grids
   char one_filename[4096];
   sprintf(one_filename, "%s/%s/%s1", output_directory, prediction_method, root_filename);
   char two_filename[4096];
   sprintf(two_filename, "%s/%s/%s2", output_directory, prediction_method, root_filename);
   char three_filename[4096];
   sprintf(three_filename, "%s/%s/%s3", output_directory, prediction_method, root_filename);
   char four_filename[4096];
   sprintf(four_filename, "%s/%s/%s4", output_directory, prediction_method, root_filename);

   // read grids
   R3Grid *one = RNReadNeuronMetaRawFile(one_filename);
   if (!one) exit(-1);
   R3Grid *two = RNReadNeuronMetaRawFile(two_filename);
   if (!two) exit(-1);
   R3Grid *three = RNReadNeuronMetaRawFile(three_filename);
   if (!three) exit(-1);
   R3Grid *four = RNReadNeuronMetaRawFile(four_filename);
   if (!four) exit(-1);

   // perform vertical merges
   if (print_verbose) { printf("Merging one and two..."); fflush(stdout); }
   R3Grid *left = VerticalMerge(one, two);
   if (!left) exit(-1);
   if (print_verbose) printf("done!\n");
   if (print_verbose) { printf("Merging three and four..."); fflush(stdout); }
   R3Grid *right = VerticalMerge(three, four);
   if (!right) exit(-1);
   if (print_verbose) printf("done!\n");

   // perform horizontal merge
   if (print_verbose) { printf("Merging left and right..."); fflush(stdout); }
   R3Grid *merged = HorizontalMerge(left, right);
   if (!merged) exit(-1);
   if (print_verbose) printf("done!\n");

   // free memory
   delete one;
   delete two;
   delete three;
   delete four;
   delete left;
   delete right;

   RNMeta meta = RNMeta("Uint16", merged->XResolution(), merged->YResolution(), merged->ZResolution(), 1);

   char output_filename[4096];
   sprintf(output_filename, "%s/%s/%s", output_directory, prediction_method, root_filename);

   // write meta/raw file
   if (!RNWriteNeuronMetaRawFile(output_filename, meta, merged)) exit(-1);

   // get truth
   char truth_filename[4096];
   sprintf(truth_filename, "neuron_data/human_labels/human_labels");
   truth_grid = RNReadNeuronMetaRawFile(truth_filename);
   if (!truth_grid) exit(-1);

   // see how stitching performed
   int ncorrect_boundaries = 0;
   int nincorrect_boundaries = 0;
   for (int ix = 0; ix < merged->XResolution(); ++ix) {
      for (int iz = 0; iz < merged->ZResolution(); ++iz) {
         int bottomy = merged->YResolution() / 2 - 1;
         int topy = merged->YResolution() / 2;

         // get the values on the boundary
         int bottom_value = (int)(merged->GridValue(ix, bottomy, iz) + 0.5);
         int top_value = (int)(merged->GridValue(ix, topy, iz) + 0.5);

         if (bottom_value == top_value) {
            int bottom_truth = (int)(truth_grid->GridValue(ix, bottomy, iz) + 0.5);
            int top_truth = (int)(truth_grid->GridValue(ix, topy, iz) + 0.5);

            if (bottom_truth == top_truth) {
               ncorrect_boundaries++;
            }
            else {
               nincorrect_boundaries++;
            }
         }
      }
   }
   for (int iy = 0; iy < merged->YResolution(); ++iy) {
      for (int iz = 0; iz < merged->ZResolution(); ++iz) {
         int leftx = merged->XResolution() / 2 - 1;
         int rightx = merged->XResolution() / 2;

         // get the values on the boundary
         int left_value = (int)(merged->GridValue(leftx, iy, iz) + 0.5);
         int right_value = (int)(merged->GridValue(rightx, iy, iz) + 0.5);

         if (left_value == right_value) {
            int left_truth = (int)(truth_grid->GridValue(leftx, iy, iz) + 0.5);
            int right_truth = (int)(truth_grid->GridValue(rightx, iy, iz) + 0.5);

            if (left_truth == right_truth) {
               ncorrect_boundaries++;
            }
            else {
               nincorrect_boundaries++;
            }
         }
      }
   }

   printf("%d %d (%lf)\n", ncorrect_boundaries, nincorrect_boundaries, ncorrect_boundaries / (RNScalar)(ncorrect_boundaries + nincorrect_boundaries));

   // free memory
   delete merged;
   delete truth_grid;

   // return success
   return 0;
}
