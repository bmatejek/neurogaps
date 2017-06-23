// Source file for the graph cut algorithm



// include files

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"



// program arguments

static char *input_name = NULL;
static const char *unary_extension = NULL;
static const char *binary_extension = NULL;
static int main_boundary_file = 1080;
static int print_verbose = 0;
static int print_debug = 0;



// program variables

static NeuronData *nd = NULL;



// directory structure

static const char *results_directory = "results/graphcut";
static const char *tmp_results_directory = "results/graphcut/tmp";
static const char *visuals_directory = "visuals/graphcut";
static const char *boundaries_directory = "boundaries";
static const char *features_directory = "algs_data/features";



// constants for graphcut

static const int NCUTS = 47;
static const int cuts[NCUTS] = {
   1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60,
   72, 90, 120, 180, 360, 720, 1080, 1440, 1800, 2160, 2880, 3240, 3600,
   4320, 5400, 6480, 7200, 8640, 10800, 12960, 14400, 16200, 21600, 25920,
   32400, 43200, 64800, 129600
};



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



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate neuron data\n");
      return 0;
   }

   nd->ReadFile(input_name);

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }

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
      test_plots[im]->SetXLabel("log2(Unary Weight)");
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
      segmentation_plots[im]->SetXLabel("log2(Unary Weight)");
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
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-unary_random_forest")) unary_extension = "random_forest";
         else if (!strcmp(*argv, "-unary_conservative")) unary_extension = "conservative";
         else if (!strcmp(*argv, "-unary_truth")) unary_extension = "truth";
         else if (!strcmp(*argv, "-binary_random_forest")) binary_extension = "random_forest_smoothing";
         else if (!strcmp(*argv, "-binary_random_forest_global")) binary_extension = "random_forest_smoothing_global";
         else if (!strcmp(*argv, "-binary_truth")) binary_extension = "truth";
	 else if (!strcmp(*argv, "-unary_ordermerge")) unary_extension = "ordermerge";
	 else if (!strcmp(*argv, "-main_boundary_file")) { argv++; argc--; main_boundary_file = atoi(*argv); }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // check filenames
   if (!input_name) { fprintf(stderr, "Must include input neuron data filename\n"); return 0; }
   if (!unary_extension) { fprintf(stderr, "Must include input unary extension\n"); return 0; }
   if (!binary_extension) { fprintf(stderr, "Must include input binary extension\n"); return 0; }

   // return OK status
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in neuron file
   if (!ReadData()) exit(-1);

   // get init, unary, and binary files
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strchr(root_filename, '.');
   *extp = '\0';

   // number of output metrics
   int ntest_metrics = nd->NTestMetrics();
   int nsegmentation_metrics = nd->NSegmentationMetrics();

   // pyimage root filename
   char output_filename[4096];
   sprintf(output_filename, "%s_%s_%s_aggregate", root_filename, unary_extension, binary_extension);

   // create the images
   if (!CreateMetricImages(root_filename)) return 0;

   // keep track of when boundaries are split
   int *first_split = new int[nd->NBoundaries()];
   int *last_split = new int[nd->NBoundaries()];
   float *proportion_split = new float[nd->NBoundaries()];
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      first_split[ib] = -1;
      last_split[ib] = -1;
      proportion_split[ib] = 0.0;
   }

   for (int ic = 0; ic < NCUTS; ++ic) {
      int unary_weight = cuts[ic];

      // read in tmp file
      char metric_filename[4096];
      sprintf(metric_filename, "%s/%s_%s_%s_%06d.rand", tmp_results_directory, root_filename, unary_extension, binary_extension, unary_weight);

      // read file
      FILE *rand_fp = fopen(metric_filename, "rb");
      if (!rand_fp) { fprintf(stderr, "Failed to open %s\n", metric_filename); exit(-1); }

      int ninput_test_metrics;
      int ninput_segmentation_metrics;
      RNScalar input_unary_weight;
      RNScalar input_binary_weight;
      RNScalar *test_metrics = new RNScalar[ntest_metrics];
      RNScalar *segmentation_metrics = new RNScalar[nsegmentation_metrics];

      // read in results
      fread(&ninput_test_metrics, sizeof(int), 1, rand_fp);
      rn_assertion(ninput_test_metrics == ntest_metrics);
      fread(&ninput_segmentation_metrics, sizeof(int), 1, rand_fp);
      rn_assertion(ninput_segmentation_metrics == nsegmentation_metrics);
      fread(&input_unary_weight, sizeof(RNScalar), 1, rand_fp);
      fread(&input_binary_weight, sizeof(RNScalar), 1, rand_fp);
      fread(test_metrics, sizeof(RNScalar), ntest_metrics, rand_fp);
      fread(segmentation_metrics, sizeof(RNScalar), nsegmentation_metrics, rand_fp);

      for (int im = 0; im < ntest_metrics; ++im) {
         test_plots[im]->InsertPoint(R2Point(log2(cuts[ic]), test_metrics[im]));
      }
      for (int im = 0; im < nsegmentation_metrics; ++im) {
         segmentation_plots[im]->InsertPoint(R2Point(log2(cuts[ic]), segmentation_metrics[im]));
      }

      // close file
      fclose(rand_fp);

      // free memory
      delete[] test_metrics;
      delete[] segmentation_metrics;


      // get boundary filename
      char boundaries_filename[4096];
      sprintf(boundaries_filename, "%s/%s_graphcut_%06d.boundary", boundaries_directory, root_filename, cuts[ic]);

      // open file
      FILE *boundary_fp = fopen(boundaries_filename, "rb");
      if (!boundary_fp) { fprintf(stderr, "Failed to read %s\n", boundaries_filename); return 0; }

      // read the number of boundaries
      int ninput_boundaries;
      fread(&ninput_boundaries, sizeof(int), 1, boundary_fp);
      rn_assertion(ninput_boundaries == nd->NBoundaries());

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
            if (first_split[ib] == -1) first_split[ib] = cuts[ic];
            last_split[ib] = cuts[ic];
            proportion_split[ib] += 1.0 / (RNScalar)NCUTS;
         }
      }

      // save if main boundary file
      if (cuts[ic] == main_boundary_file) {
	 // get boundary filename
	 char main_boundaries_filename[4096];
	 sprintf(main_boundaries_filename, "%s/%s_graphcut.boundary", boundaries_directory, root_filename);

	 // open file
	 FILE *main_boundary_fp = fopen(main_boundaries_filename, "wb");
	 if (!main_boundary_fp) { fprintf(stderr, "Failed to write to %s\n", main_boundaries_filename); return 0; }

	 // write the number of boundaries
	 fwrite(&ninput_boundaries, sizeof(int), 1, main_boundary_fp);
	 
	 // write the boundary results
	 fwrite(boundary_results, sizeof(RNBoolean), nd->NBoundaries(), main_boundary_fp);
	 
	 // close file
	 fclose(main_boundary_fp);
      }

      // free memory
      delete[] boundary_results;
   }

   // write the images
   if (!WriteMetricImages(output_filename)) return 0;

   // find variables that are not split ever
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      if (first_split[ib] == -1) first_split[ib] = cuts[NCUTS - 1] + 1;
      if (last_split[ib] == -1) last_split[ib] = cuts[NCUTS - 1] + 1;
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
   sprintf(first_split_filename, "%s/%s_graphcut_first_split.feature", features_directory, root_filename);
   char last_split_filename[4096];
   sprintf(last_split_filename, "%s/%s_graphcut_last_split.feature", features_directory, root_filename);
   char proportion_filename[4096];
   sprintf(proportion_filename, "%s/%s_graphcut_proportion.feature", features_directory, root_filename);

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
   delete nd;

   // return success
   return 0;
}
