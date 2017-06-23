// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include "RNDataStructures/RNDataStructures.h"
#include "RNML/RNML.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int threshold = 50;
static int validation = 0;
static int merge_validation = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static int *cellular_mapping = NULL;
static float *forest_results = NULL;
static int *voxel_proposals = NULL;



// directory structure

static const char *graphcut_directory = "algs_data/graphcut";
static const char *boundaries_directory = "boundaries";
static const char *features_directory = "algs_data/features";
static const char *output_directory = "output/boundarymerge";
static const char *results_directory = "results/boundarymerge";
static const char *visuals_directory = "visuals/boundarymerge";



// string constants for image output

static const int ntest_images = 4;
static const char *test_suffixes[ntest_images] = {
   "rand_error", "fscore", "variation_fscore", "variation_information"
};



static const int nsegmentation_images = 2;
static const char *segmentation_suffixes[nsegmentation_images] = {
   "supervoxel_occurrences", "supervoxel_precision_recall"
};



static RNPyImage *test_images[ntest_images];
static RNPyPointPlot **test_plots = NULL;
static RNPyImage *segmentation_images[nsegmentation_images];
static RNPyPointPlot **segmentation_plots = NULL;



// useful structures
struct BMBoundaryRank {
   BMBoundaryRank(NeuronBoundary *boundary, RNScalar forest_score) :
      boundary(boundary),
      forest_score(1.0 - forest_score)
   {
   }

   // instance variabels
   NeuronBoundary *boundary;
   RNScalar forest_score;
};



// vector of sorted boundaries
static std::vector<BMBoundaryRank> boundaries = std::vector<BMBoundaryRank>();



int BoundarySortFunction(BMBoundaryRank boundary_one, BMBoundaryRank boundary_two)
{
   return boundary_one.forest_score > boundary_two.forest_score;
}



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



static int ReadRandomForestScores(char root_filename[4096])
{
   // get filename
   char graphcut_filename[4096];
   sprintf(graphcut_filename, "%s/%s_random_forest_smoothing.binary", graphcut_directory, root_filename);

   // open file
   FILE *graph_fp = fopen(graphcut_filename, "rb");
   if (!graph_fp) { fprintf(stderr, "Failed to read %s\n", graphcut_filename); return 0; }

   // read in merge order
   int ninput_boundaries;
   fread(&ninput_boundaries, sizeof(int), 1, graph_fp);
   rn_assertion(ninput_boundaries == nd->NBoundaries());

   // allocate memory for merge order
   forest_results = new float[ninput_boundaries];
   fread(forest_results, sizeof(float), nd->NBoundaries(), graph_fp);

   // close file
   fclose(graph_fp);

   // return success
   return 1;
}



static int CreateMetricImages(const char root_filename[4096])
{
   // create all output images
   for (int ii = 0; ii < ntest_images; ++ii) {
      test_images[ii] = new RNPyImage();
      test_images[ii]->SetOutputDirectory(visuals_directory);
   }
   for (int ii = 0; ii < nsegmentation_images; ++ii) {
      segmentation_images[ii] = new RNPyImage();
      segmentation_images[ii]->SetOutputDirectory(visuals_directory);
   }

   // create all of the individual test plots
   int ntest_metrics = nd->NTestMetrics();
   test_plots = new RNPyPointPlot *[ntest_metrics];
   for (int im = 0; im < ntest_metrics; ++im) {
      test_plots[im] = new RNPyPointPlot();

      // set title and labels
      char title[4096];
      sprintf(title, "%s %s", root_filename, nd->TestMetricName(im));
      test_plots[im]->SetTitle(title);
      test_plots[im]->SetXLabel("Number of Merges");
      test_plots[im]->SetYLabel(nd->TestMetricName(im));

      // create the legend
      test_plots[im]->SetLegend(nd->TestMetricName(im));

      // add these plots to the proper images
      test_images[im / 3]->InsertPyPlot(test_plots[im]);

      // set extrema labels for certain full plots
      if (im == 0 || im == 9)
         test_plots[im]->SetExtremaType(MIN_EXTREMA);
      else if (im == 3 || im == 6)
         test_plots[im]->SetExtremaType(MAX_EXTREMA);
   }

   // create all of the individual segmentation plots
   int nsegmentation_metrics = nd->NSegmentationMetrics();
   segmentation_plots = new RNPyPointPlot *[nsegmentation_metrics];
   for (int im = 0; im < nsegmentation_metrics; ++im) {
      segmentation_plots[im] = new RNPyPointPlot();

      // set title, labels, and legends
      char title[4096];
      char ylabel[4096];
      sprintf(title, "%s %s", root_filename, nd->SegmentationMetricName(im));
      sprintf(ylabel, "%s", nd->SegmentationMetricName(im));

      segmentation_plots[im]->SetTitle(title);
      segmentation_plots[im]->SetXLabel("Number of Merges");
      segmentation_plots[im]->SetYLabel(ylabel);

      // create the legend
      segmentation_plots[im]->SetLegend(nd->SegmentationMetricName(im));

      // add this plot to the correct image
      segmentation_images[im / 4]->InsertPyPlot(segmentation_plots[im]);

      // set extrema labels for certain full plots
      if (im == 1 || im == 2)
         segmentation_plots[im]->SetExtremaType(MIN_EXTREMA);
      else
         segmentation_plots[im]->SetExtremaType(MAX_EXTREMA);
   }

   // return success
   return 1;
}



static int WriteMetricImages(const char root_filename[4096])
{
   // save the images
   for (int ii = 0; ii < ntest_images; ++ii) {
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_%s.pyimage", results_directory, root_filename, test_suffixes[ii]);
      test_images[ii]->WriteImageFile(output_filename);
   }
   for (int ii = 0; ii < nsegmentation_images; ++ii) {
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_%s.pyimage", results_directory, root_filename, segmentation_suffixes[ii]);
      segmentation_images[ii]->WriteImageFile(output_filename);
   }

   // delete images
   for (int ii = 0; ii < ntest_images; ++ii)
      delete test_images[ii];
   for (int ii = 0; ii < nsegmentation_images; ++ii)
      delete segmentation_images[ii];

   // delete plots
   for (int im = 0; im < nd->NTestMetrics(); ++im)
      delete test_plots[im];
   delete[] test_plots;
   for (int im = 0; im < nd->NSegmentationMetrics(); ++im)
      delete segmentation_plots[im];
   delete[] segmentation_plots;

   // return success
   return 1;
}



static int WriteMetaRawFile(char root_filename[4096])
{
   // get filename
   char meta_filename[4096];
   sprintf(meta_filename, "%s/%s", output_directory, root_filename);

   // create output grid
   R3Grid *output_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      output_grid->SetGridValue(iv, voxel_proposals[iv]);
   }

   // write file
   RNMeta meta = RNMeta("Uint16", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 1);
   if (!RNWriteNeuronMetaRawFile(meta_filename, meta, output_grid)) return 0;

   // free memory
   delete output_grid;

   // return success
   return 1;
}



static int WriteBoundaryFile(char root_filename[4096])
{
   // output results to boundary folder
   char boundary_filename[4096];
   if (validation) sprintf(boundary_filename, "%s/%s_boundarymerge_%03d.boundary", boundaries_directory, root_filename, threshold);
   else sprintf(boundary_filename, "%s/%s_boundarymerge.boundary", boundaries_directory, root_filename);

   // open file
   FILE *boundaries_fp = fopen(boundary_filename, "wb");
   if (!boundaries_fp) { fprintf(stderr, "Failed to write to %s\n", boundary_filename); return 0; }

   int nboundaries = nd->NBoundaries();
   fwrite(&nboundaries, sizeof(int), 1, boundaries_fp);

   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get supervoxels
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      // get data indices
      int supervoxel_one_data_index = supervoxel_one->DataIndex();
      int supervoxel_two_data_index = supervoxel_two->DataIndex();

      RNBoolean merged = FALSE;
      if (!supervoxel_one->IsExtracellular() && !supervoxel_two->IsExtracellular()) {
         // get labels for both cellulars
         int cellular_one_label = cellular_mapping[supervoxel_one_data_index];
         int cellular_two_label = cellular_mapping[supervoxel_two_data_index];

         merged = (cellular_one_label == cellular_two_label);
      }

      // write the result
      fwrite(&merged, sizeof(RNBoolean), 1, boundaries_fp);
   }

   // close file
   fclose(boundaries_fp);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Algorithmic functions
////////////////////////////////////////////////////////////////////////

static int GenerateRanking(void)
{
   // sort smoothing terms
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      if (boundary->SupervoxelOne()->IsExtracellular()) continue;
      if (boundary->SupervoxelTwo()->IsExtracellular()) continue;

      //rn_assertion(boundary_path_order[ib] != -1);
      BMBoundaryRank boundary_rank = BMBoundaryRank(boundary, forest_results[ib]);
      boundaries.push_back(boundary_rank);
   }

   // sort vector
   std::sort(boundaries.begin(), boundaries.end(), BoundarySortFunction);

   // return success
   return 1;
}



static int MergeResults(void)
{
   // count number of cellular boundaries
   int ncellular_boundaries = 0;
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);
      if (boundary->SupervoxelOne()->IsExtracellular()) continue;
      if (boundary->SupervoxelTwo()->IsExtracellular()) continue;

      ncellular_boundaries++;
   }

   RNBoolean *boundary_hash = new RNBoolean[nd->NBoundaries()];
   for (int ib = 0; ib < nd->NBoundaries(); ++ib)
      boundary_hash[ib] = FALSE;

   RNBoolean *cellular_hash = new RNBoolean[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic)
      cellular_hash[ic] = FALSE;

   BMBoundaryRank tmp = BMBoundaryRank(NULL, -1.0);
   RNMinBinaryHeap<BMBoundaryRank *> boundary_ranks = RNMinBinaryHeap<BMBoundaryRank *>(&tmp, &(tmp.forest_score), nd->NBoundaries());

   // generate data
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      if (!cellular->IsOnBoundary()) continue;

      // go through all cellular boundaries
      for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = cellular->Boundary(ib);

         // make sure this belongs to another cellular
         if (boundary->OtherSupervoxel(cellular)->IsExtracellular()) continue;

         int boundary_data_index = boundary->DataIndex();

         if (boundary_hash[boundary_data_index]) continue;

         float forest_result = forest_results[boundary_data_index];

         // add to vector of boundary ranks
         BMBoundaryRank *boundary_rank = new BMBoundaryRank(boundary, forest_result);

         boundary_ranks.Insert(boundary_data_index, boundary_rank);

         // set boundary hash
         boundary_hash[boundary_data_index] = TRUE;
      }
      // set cellular hash to true since all boundaries of this have been added
      cellular_hash[ic] = TRUE;
   }

   // set initial cellular mapping
   cellular_mapping = new int[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      cellular_mapping[ic] = ic;
   }

   // start statistics
   RNTime start_time;
   start_time.Read();

   // go through all merge
   if (print_verbose) { printf("Merging cellulars for threshold %0.2lf...\n  ", (RNScalar)threshold / 100.0); fflush(stdout); }
   int nmerged = 0;
   while (!boundary_ranks.IsEmpty()) {
      if (print_verbose) { RNProgressBar(nmerged, ncellular_boundaries); }
      nmerged++;

      // see which cellulars belong to this boundary
      BMBoundaryRank *boundary_rank = boundary_ranks.DeleteMin();

      // get the boundary for merging
      NeuronBoundary *boundary = boundary_rank->boundary;

      // just checking
      rn_assertion(boundary->SupervoxelOne()->IsCellular());
      rn_assertion(boundary->SupervoxelTwo()->IsCellular());

      // get the cellular indicies
      int cellular_one_data_index = boundary->SupervoxelOne()->DataIndex();
      int cellular_two_data_index = boundary->SupervoxelTwo()->DataIndex();

      NeuronCellular *cellular_one = nd->Cellular(cellular_one_data_index);
      NeuronCellular *cellular_two = nd->Cellular(cellular_two_data_index);

      if (!cellular_hash[cellular_one_data_index]) {
         // go through all cellular boundaries
         for (int ib = 0; ib < cellular_one->NBoundaries(); ++ib) {
            NeuronBoundary *boundary = cellular_one->Boundary(ib);
            int boundary_data_index = boundary->DataIndex();

            // make sure this belongs to another cellular
            if (boundary->OtherSupervoxel(cellular_one)->IsExtracellular()) continue;

            if (boundary_hash[boundary_data_index]) continue;

            float forest_result = forest_results[boundary_data_index];

            // add to vector of boundary ranks
            BMBoundaryRank *boundary_rank = new BMBoundaryRank(boundary, forest_result);

            boundary_ranks.Insert(boundary_data_index, boundary_rank);

            // set boundary hash
            boundary_hash[boundary_data_index] = TRUE;
         }
         // set cellular hash to true since all boundaries of this have been added
         cellular_hash[cellular_one_data_index] = TRUE;
      }
      if (!cellular_hash[cellular_two_data_index]) {
         // go through all cellular boundaries
         for (int ib = 0; ib < cellular_two->NBoundaries(); ++ib) {
            NeuronBoundary *boundary = cellular_two->Boundary(ib);
            int boundary_data_index = boundary->DataIndex();

            // make sure this belongs to another cellular
            if (boundary->OtherSupervoxel(cellular_two)->IsExtracellular()) continue;

            if (boundary_hash[boundary_data_index]) continue;

            float forest_result = forest_results[boundary_data_index];

            // add to vector of boundary ranks
            BMBoundaryRank *boundary_rank = new BMBoundaryRank(boundary, forest_result);

            boundary_ranks.Insert(boundary_data_index, boundary_rank);

            // set boundary hash
            boundary_hash[boundary_data_index] = TRUE;
         }
         // set cellular hash to true since all boundaries of this have been added
         cellular_hash[cellular_two_data_index] = TRUE;
      }

      // get cellular labels
      int cellular_one_label = cellular_mapping[cellular_one_data_index];
      int cellular_two_label = cellular_mapping[cellular_two_data_index];
      if (cellular_one_label == cellular_two_label) continue;

      // find the number of agreeing and disagreeing boundaries
      int nagreements = 0;
      int ndisagreements = 0;

      // get the score
      float merge_score = boundary_rank->forest_score;
      if (merge_score > 0.5) nagreements++;
      else ndisagreements++;

      // go through all boundaries between these two agglomerated supervoxels
      for (int cellular_data_index = 0; cellular_data_index < nd->NCellulars(); ++cellular_data_index) {
         int cellular_label = cellular_mapping[cellular_data_index];
         // only consider if this cellular has the same mapping as cellular one
         if (cellular_label == cellular_one_label) {
            NeuronCellular *cellular = nd->Cellular(cellular_data_index);
            // iterate through all boundaries
            for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
               NeuronBoundary *cellular_boundary = cellular->Boundary(ib);

               NeuronSupervoxel *neighbor_cellular = cellular_boundary->OtherSupervoxel(cellular);
               if (neighbor_cellular->IsExtracellular()) continue;

               // get the label for this cellular index
               int neighbor_cellular_data_index = neighbor_cellular->DataIndex();
               int neighbor_cellular_label = cellular_mapping[neighbor_cellular_data_index];

               // see if this matches label two
               if (neighbor_cellular_label == cellular_two_label) {
                  if (forest_results[cellular_boundary->DataIndex()] > 0.5) nagreements++;
                  else ndisagreements++;
               }
            }
         }
      }

      // see if the cellulars should be merged
      RNScalar ratio = nagreements / (RNScalar)(nagreements + ndisagreements);
      if (ratio > threshold / (RNScalar)100) {
         for (int ic = 0; ic < nd->NCellulars(); ++ic) {
            if (cellular_mapping[ic] == cellular_two_label)
               cellular_mapping[ic] = cellular_one_label;
         }
      }
   }
   if (print_verbose) printf("\ndone in %0.2f seconds!\n", start_time.Elapsed());

   // deflate cellular mapping
   RNDeflateIntegerArray(cellular_mapping, nd->NCellulars());

   // return success
   return 1;
}



static int GenerateStatistics(FILE *output_fp = NULL)
{
   // get voxel proposals
   int nvoxels = nd->NVoxels();
   voxel_proposals = new int[nvoxels];
   for (int iv = 0; iv < nvoxels; ++iv) {
      int voxel_supervoxel = nd->VoxelMapping(iv);
      if (voxel_supervoxel >= nd->NCellulars()) voxel_proposals[iv] = 0;
      else voxel_proposals[iv] = cellular_mapping[voxel_supervoxel] + 1;
   }

   RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
   RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];

   // run the test metric
   if (!nd->TestMetric(voxel_proposals, test_metrics, !validation)) return 0;
   if (!nd->SegmentationMetric(cellular_mapping, segmentation_metrics)) return 0;


   if (validation) {
      rn_assertion(output_fp != NULL);
      // print out the full rand score
      fprintf(output_fp, "%0.2lf", (RNScalar)threshold / 100.0);
      for (int it = 0; it < nd->NTestMetrics(); ++it) {
         fprintf(output_fp, ",%0.6lf", test_metrics[it]);
      }
      for (int is = 0; is < nd->NSegmentationMetrics(); ++is) {
         fprintf(output_fp, ",%0.6lf", segmentation_metrics[is]);
      }
      fprintf(output_fp, "\n");
      fflush(output_fp);
   }

   // free memory
   delete[] test_metrics;
   delete[] segmentation_metrics;
   if (validation) {
      delete[] cellular_mapping;
      delete[] voxel_proposals;
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
         else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
         else if (!strcmp(*argv, "-merge")) merge_validation = 1;
         else if (!strcmp(*argv, "-validation")) validation = 1;
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

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in random forest scores
   if (!ReadRandomForestScores(root_filename)) exit(-1);

   // read voxels
   nd->ReadVoxels();

   // cellular agreements
   RNBoolean **cellular_agreement = new RNBoolean *[nd->NCellulars()];
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      cellular_agreement[ic1] = new RNBoolean[nd->NCellulars()];
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2)
         cellular_agreement[ic1][ic2] = FALSE;
   }

   // go through all boundary cellulars
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      if (!cellular->IsOnBoundary()) continue;

      BMBoundaryRank tmp = BMBoundaryRank(NULL, 0.0);

      // create hash table for cellulars that belong
      RNBoolean *cellular_hash = new RNBoolean[nd->NCellulars()];
      for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1)
         cellular_hash[ic1] = FALSE;
      cellular_hash[ic] = TRUE;

      // start merging cellulars
      RNMinBinaryHeap<BMBoundaryRank *> heap = RNMinBinaryHeap<BMBoundaryRank *>(&tmp, &(tmp.forest_score), nd->NBoundaries());

      // add all boundaries to this cellular
      for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = cellular->Boundary(ib);
         heap.Insert(boundary->DataIndex(), new BMBoundaryRank(boundary, forest_results[boundary->DataIndex()]));
      }

      // take all boundaries off
      while (!heap.IsEmpty()) {
         BMBoundaryRank *boundary_rank = heap.DeleteMin();

         // should this cellular be added
         NeuronBoundary *boundary = boundary_rank->boundary;
         NeuronCellular *cellular_one = (NeuronCellular *)boundary->SupervoxelOne();
         NeuronCellular *cellular_two = (NeuronCellular *)boundary->SupervoxelTwo();

         if (cellular_hash[cellular_one->DataIndex()]) {
            NeuronCellular *neighbor = cellular_two;
            if (boundary_rank->forest_score < 0.5) {
               // add all of these boundaries
               for (int ib = 0; ib < neighbor->NBoundaries(); ++ib) {
                  NeuronCellular *other = (NeuronCellular *)neighbor->Boundary(ib)->OtherSupervoxel(neighbor);
                  if (cellular_hash[other->DataIndex()]) continue;
                  heap.Insert(neighbor->Boundary(ib)->DataIndex(), new BMBoundaryRank(neighbor->Boundary(ib), forest_results[neighbor->Boundary(ib)->DataIndex()]));
               }
            }
         }
         else if (cellular_hash[cellular_two->DataIndex()]) {
            NeuronCellular *neighbor = cellular_one;
            if (boundary_rank->forest_score < 0.5) {
               // add all of these boundaries
               for (int ib = 0; ib < neighbor->NBoundaries(); ++ib) {
                  NeuronCellular *other = (NeuronCellular *)neighbor->Boundary(ib)->OtherSupervoxel(neighbor);
                  if (cellular_hash[other->DataIndex()]) continue;
                  heap.Insert(neighbor->Boundary(ib)->DataIndex(), new BMBoundaryRank(neighbor->Boundary(ib), forest_results[neighbor->Boundary(ib)->DataIndex()]));
               }
            }
         }
      }

      for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
         if (cellular_hash[ic1]) cellular_agreement[ic1][ic] = TRUE;
      }

      // free memory
      delete[] cellular_hash;
   }

   // how many supervoxels are never wrong

   // free up memory
   delete nd;
   delete[] forest_results;

   // return success
   return 0;
}
