// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int threshold = 50;
static int validation = 0;
static int merge_validation = 0;
static int smoothing_global = 0;



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
static const char *output_directory = "output/ordermerge";
static const char *results_directory = "results/ordermerge";
static const char *visuals_directory = "visuals/ordermerge";



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
struct OMBoundaryRank {
   OMBoundaryRank(NeuronBoundary *boundary, float forest_score) :
      boundary(boundary),
      forest_score(forest_score)
   {
   }

   // instance variabels
   NeuronBoundary *boundary;
   float forest_score;
};



// vector of sorted boundaries
static std::vector<OMBoundaryRank> boundaries = std::vector<OMBoundaryRank>();



// useful sorting functions
int BoundarySortFunction(OMBoundaryRank boundary_one, OMBoundaryRank boundary_two)
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
   if (smoothing_global) sprintf(graphcut_filename, "%s/%s_random_forest_smoothing_global.binary", graphcut_directory, root_filename);
   else sprintf(graphcut_filename, "%s/%s_random_forest_smoothing.binary", graphcut_directory, root_filename);

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
      test_plots[im]->SetXLabel("Voting Threshold");
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
      segmentation_plots[im]->SetXLabel("Voting Threshold");
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
      if (smoothing_global) sprintf(output_filename, "%s/%s_%s_smoothing_global.pyimage", results_directory, root_filename, test_suffixes[ii]);
      else sprintf(output_filename, "%s/%s_%s.pyimage", results_directory, root_filename, test_suffixes[ii]);
      test_images[ii]->WriteImageFile(output_filename);
   }
   for (int ii = 0; ii < nsegmentation_images; ++ii) {
      char output_filename[4096];
      if (smoothing_global) sprintf(output_filename, "%s/%s_%s_smoothing_global.pyimage", results_directory, root_filename, segmentation_suffixes[ii]);
      else sprintf(output_filename, "%s/%s_%s.pyimage", results_directory, root_filename, segmentation_suffixes[ii]);
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
   if (smoothing_global) sprintf(meta_filename, "%s/%s_smoothing_global", output_directory, root_filename);
   else sprintf(meta_filename, "%s/%s", output_directory, root_filename);

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
   if (smoothing_global) {
      if (validation) sprintf(boundary_filename, "%s/%s_ordermerge_smoothing_global_%03d.boundary", boundaries_directory, root_filename, threshold);
      else sprintf(boundary_filename, "%s/%s_ordermerge_smoothing_global.boundary", boundaries_directory, root_filename);
   }
   else {
      if (validation) sprintf(boundary_filename, "%s/%s_ordermerge_%03d.boundary", boundaries_directory, root_filename, threshold);
      else sprintf(boundary_filename, "%s/%s_ordermerge.boundary", boundaries_directory, root_filename);
   }

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
      OMBoundaryRank boundary_rank = OMBoundaryRank(boundary, forest_results[ib]);
      boundaries.push_back(boundary_rank);
   }

   // sort vector
   std::sort(boundaries.begin(), boundaries.end(), BoundarySortFunction);

   // return success
   return 1;
}



static int MergeResults(void)
{
   // set initial cellular mapping
   cellular_mapping = new int[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      cellular_mapping[ic] = ic;
   }

   // start statistics
   RNTime start_time;
   start_time.Read();

   // go through all merge
   if (print_verbose) { printf("Merging cellulars...\n  "); fflush(stdout); }
   for (int im = 0; im < (int)boundaries.size(); ++im) {
      if (print_verbose) RNProgressBar(im, boundaries.size());
      int boundary_data_index = boundaries[im].boundary->DataIndex();

      // get the boundary ready for merging
      NeuronBoundary *boundary = nd->Boundary(boundary_data_index);

      // just checking
      rn_assertion(boundary->SupervoxelOne()->IsCellular());
      rn_assertion(boundary->SupervoxelTwo()->IsCellular());

      // get the cellular indices
      int cellular_one_data_index = boundary->SupervoxelOne()->DataIndex();
      int cellular_two_data_index = boundary->SupervoxelTwo()->DataIndex();

      int cellular_one_label = cellular_mapping[cellular_one_data_index];
      int cellular_two_label = cellular_mapping[cellular_two_data_index];

      // find the number of agreeing and disagreeing boundaries
      int nagreements = 0;
      int ndisagreements = 0;

      // get the score
      float merge_score = boundaries[im].forest_score;
      if (merge_score > 0.5) nagreements++;
      else ndisagreements++;

      // go through all boundaries between these two agglomerated supervoxels
      for (int cellular_data_index = 0; cellular_data_index < nd->NCellulars(); ++cellular_data_index) {
         int cellular_label = cellular_mapping[cellular_data_index];
         // only consider if this cellular has the same mapping has cellular_one
         if (cellular_label == cellular_one_label) {
            NeuronCellular *cellular = nd->Cellular(cellular_data_index);
            // iterate through all boundaries
            for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
               NeuronBoundary *cellular_boundary = cellular->Boundary(ib);
               NeuronSupervoxel *neighbor_cellular = cellular_boundary->OtherSupervoxel(cellular);
               if (neighbor_cellular->IsExtracellular()) continue;

               // get the label for this cellular
               int neighbor_cellular_data_index = neighbor_cellular->DataIndex();
               int neighbor_cellular_label = cellular_mapping[neighbor_cellular_data_index];

               // see if this matches cellular label two
               if (neighbor_cellular_label == cellular_two_label) {
                  if (forest_results[cellular_boundary->DataIndex()] > 0.5) nagreements++;
                  else ndisagreements++;
               }
            }
         }
      }

      // see if these cellulars should be merged
      RNScalar ratio = nagreements / (RNScalar)(nagreements + ndisagreements);
      if (ratio > threshold / (RNScalar)100) {
         for (int ic = 0; ic < nd->NCellulars(); ++ic) {
            if (cellular_mapping[ic] == cellular_two_label) {
               cellular_mapping[ic] = cellular_one_label;
            }
         }
      }
   }
   if (print_verbose) printf("\ndone with threshold %0.2lf in %0.2lf seconds!\n", (RNScalar)threshold / 100.0, start_time.Elapsed());

   // deflate the cellular mapping
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

      printf("Score of %0.6lf\n", test_metrics[3]);
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
         else if (!strcmp(*argv, "-smoothing_global")) smoothing_global = 1;
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

   if (smoothing_global) printf("Using smoothing terms with global features.\n");
   else printf("Using only smoothing terms.\n");

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

   if (merge_validation) {
      // create the metric image
      if (!CreateMetricImages(root_filename)) exit(-1);

      // open results file
      char input_result_filename[4096];
      if (smoothing_global) sprintf(input_result_filename, "%s/%s_smoothing_global.results", results_directory, root_filename);
      else sprintf(input_result_filename, "%s/%s.results", results_directory, root_filename);

      // open the file 
      FILE *input_fp = fopen(input_result_filename, "r");
      if (!input_filename) { fprintf(stderr, "Failed to read %s\n", input_result_filename); exit(-1); }

      for (int it = 0; it <= 100; ++it) {
         RNScalar threshold;
         fscanf(input_fp, "%lf", &threshold);
         for (int im = 0; im < nd->NTestMetrics(); ++im) {
            RNScalar test_metric;
            fscanf(input_fp, ",%lf", &test_metric);
            test_plots[im]->InsertPoint(R2Point(threshold, test_metric));
         }
         for (int im = 0; im < nd->NSegmentationMetrics(); ++im) {
            RNScalar segmenation_metric;
            fscanf(input_fp, ",%lf", &segmenation_metric);
            segmentation_plots[im]->InsertPoint(R2Point(threshold, segmenation_metric));
         }
         fscanf(input_fp, "\n");
      }

      // close file
      fclose(input_fp);

      // write the image
      if (!WriteMetricImages(root_filename)) exit(-1);

      // keep track of when boundaries are split
      int *first_split = new int[nd->NBoundaries()];
      int *last_split = new int[nd->NBoundaries()];
      float *proportion_split = new float[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         first_split[ib] = -1;
         last_split[ib] = -1;
         proportion_split[ib] = 0.0;
      }

      // go through all ordermerge results
      for (int it = 0; it <= 100; ++it) {
         RNProgressBar(it, 101);

         // get boundary filename
         char boundaries_filename[4096];
         if (smoothing_global) sprintf(boundaries_filename, "%s/%s_ordermerge_smoothing_global_%03d.boundary", boundaries_directory, root_filename, it);
         else sprintf(boundaries_filename, "%s/%s_ordermerge_%03d.boundary", boundaries_directory, root_filename, it);

         // open file
         FILE *boundary_fp = fopen(boundaries_filename, "rb");
         if (!boundary_fp) { fprintf(stderr, "Failed to read %s\n", boundaries_filename); return 0; }

         // read the number of boundaries
         int input_nboundaries;
         fread(&input_nboundaries, sizeof(int), 1, boundary_fp);
         rn_assertion(input_nboundaries == nd->NBoundaries());

         RNBoolean *boundary_results = new RNBoolean[nd->NBoundaries()];
         if (fread(boundary_results, sizeof(RNBoolean), nd->NBoundaries(), boundary_fp) != (unsigned int)nd->NBoundaries()) {
            fprintf(stderr, "Failed to read %s\n", boundaries_filename);
            return 0;
         }

         // close file
         fclose(boundary_fp);

         // go through all boundaries
         for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
            if (!boundary_results[ib]) {
               if (first_split[ib] == -1) first_split[ib] = it;
               last_split[ib] = it;
               proportion_split[ib] += 0.01;
            }
         }

         // free memory
         delete[] boundary_results;
      }
      printf("\n");

      // find variables that are not split ever
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         if (first_split[ib] == -1) first_split[ib] = 101;
         if (last_split[ib] == -1) last_split[ib] = 101;
      }

      // create dataset to test entropy
      RNDataset *dataset = new RNDataset();
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);

         // get supervoxels
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

         // skip extracellulars
         if (supervoxel_one->IsExtracellular()) continue;
         if (supervoxel_two->IsExtracellular()) continue;

         std::vector<float> attributes = std::vector<float>();
         attributes.push_back(first_split[ib]);
         attributes.push_back(last_split[ib]);
         attributes.push_back(proportion_split[ib]);

         dataset->InsertDatapoint(attributes, supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel(), ib);
      }

      std::vector<std::string> names = std::vector<std::string>();
      names.push_back("first_split");
      names.push_back("last_split");
      names.push_back("proportion_split");

      // calculate entropy
      dataset->SetNames(names);
      dataset->CalculateEntropy();

      // create output feature files
      char first_split_filename[4096];
      char last_split_filename[4096];
      char proportion_filename[4096];

      if (smoothing_global) {
         sprintf(proportion_filename, "%s/%s_ordermerge_smoothing_global_proportion.feature", features_directory, root_filename);
         sprintf(first_split_filename, "%s/%s_ordermerge_smoothing_global_first_split.feature", features_directory, root_filename);
         sprintf(last_split_filename, "%s/%s_ordermerge_smoothing_global_last_split.feature", features_directory, root_filename);
      }
      else {
         sprintf(proportion_filename, "%s/%s_ordermerge_proportion.feature", features_directory, root_filename);
         sprintf(first_split_filename, "%s/%s_ordermerge_first_split.feature", features_directory, root_filename);
         sprintf(last_split_filename, "%s/%s_ordermerge_last_split.feature", features_directory, root_filename);
      }

      // open files
      FILE *first_split_fp = fopen(first_split_filename, "wb");
      if (!first_split_fp) { fprintf(stderr, "Failed to write to %s\n", first_split_filename); return 0; }

      FILE *last_split_fp = fopen(last_split_filename, "wb");
      if (!last_split_fp) { fprintf(stderr, "Failed to write to %s\n", last_split_filename); return 0; }

      FILE *proportion_fp = fopen(proportion_filename, "wb");
      if (!proportion_fp) { fprintf(stderr, "Failed to write to %s\n", proportion_filename); return 0; }

      // write number of boundaries to files
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, first_split_fp);
      fwrite(&nboundaries, sizeof(int), 1, last_split_fp);
      fwrite(&nboundaries, sizeof(int), 1, proportion_fp);

      // write the results
      fwrite(first_split, sizeof(int), nd->NBoundaries(), first_split_fp);
      fwrite(last_split, sizeof(int), nd->NBoundaries(), last_split_fp);
      fwrite(proportion_split, sizeof(float), nd->NBoundaries(), proportion_fp);

      // close files
      fclose(first_split_fp);
      fclose(last_split_fp);
      fclose(proportion_fp);


      // free memory
      delete[] first_split;
      delete[] last_split;
      delete[] proportion_split;
   }
   else if (validation) {
      // read voxels
      nd->ReadVoxels();

      // create output filename
      char output_filename[4096];
      if (smoothing_global) sprintf(output_filename, "%s/%s_smoothing_global.results", results_directory, root_filename);
      else sprintf(output_filename, "%s/%s.results", results_directory, root_filename);

      // open the file
      FILE *output_fp = fopen(output_filename, "w");
      if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

      // generate data
      if (!GenerateRanking()) exit(-1);

      for (threshold = 0; threshold <= 100; ++threshold) {
         // start merging
         if (!MergeResults()) exit(-1);

         // see how well this mapping performs
         if (!GenerateStatistics(output_fp)) exit(-1);

         // write the boundary results
         if (!WriteBoundaryFile(root_filename)) exit(-1);
      }

      // close file
      fclose(output_fp);
   }
   else {
      // read voxels
      nd->ReadVoxels();

      // generate data
      if (!GenerateRanking()) exit(-1);

      // start merging
      if (!MergeResults()) exit(-1);

      // see how well this mapping performs
      if (!GenerateStatistics()) exit(-1);

      // output the meta/raw results
      if (!WriteMetaRawFile(root_filename)) exit(-1);

      // write the boundary results
      if (!WriteBoundaryFile(root_filename)) exit(-1);
   }

   // free up memory
   delete nd;
   delete[] forest_results;
   if (!validation) {
      delete[] cellular_mapping;
      delete[] voxel_proposals;
   }

   // return success
   return 0;
}
