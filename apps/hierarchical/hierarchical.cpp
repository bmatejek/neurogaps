// Source file to create graph cut initial positions



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int merge_results = 0;
static char *input_filename = NULL;
static const char *extension = NULL;

// get the starting index and the number of cellulars merges to happen

static int start_index = -1;
static int nsteps = -1;



// program variables

static NeuronData *nd;



// directory variables

static const char *hierarchical_directory = "algs_data/hierarchical";
static const char *tmp_hierarchical_directory = "algs_data/hierarchical/tmp";
static const char *results_directory = "results/hierarchical";
static const char *visuals_directory = "visuals/hierarchical";



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
      sprintf(title, "%s %s %s", root_filename, extension, nd->TestMetricName(im));
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
      sprintf(title, "%s %s %s", root_filename, extension, nd->SegmentationMetricName(im));
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
      sprintf(output_filename, "%s/%s_%s_%s.pyimage", results_directory, root_filename, extension, test_suffixes[ii]);
      test_images[ii]->WriteImageFile(output_filename);
   }
   for (int ii = 0; ii < nsegmentation_images; ++ii) {
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_%s_%s.pyimage", results_directory, root_filename, extension, segmentation_suffixes[ii]);
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
// Hierarchical functions
////////////////////////////////////////////////////////////////////////

static int
HierarchicalSubsection(const char root_filename[4096])
{
   // get merge filename
   char merge_filename[4096];
   sprintf(merge_filename, "%s/%s_%s.hier", hierarchical_directory, root_filename, extension);

   // open file
   FILE *merge_fp = fopen(merge_filename, "rb");
   if (!merge_fp) { fprintf(stderr, "Failed to read %s\n", merge_filename); return 0; }

   int nboundaries;
   fread(&nboundaries, sizeof(int), 1, merge_fp);

   int *boundary_merge_order = new int[nboundaries];
   fread(boundary_merge_order, sizeof(int), nboundaries, merge_fp);

   // close file
   fclose(merge_fp);

   // read voxels
   nd->ReadVoxels();

   // cellulars map to themselves before merging
   int ncellulars = nd->NCellulars();
   int *cellular_mapping = new int[ncellulars];
   for (int ic = 0; ic < ncellulars; ++ic)
      cellular_mapping[ic] = ic;

   int nmerges_occurred = 0;
   for (int ib = 0; ib < nboundaries; ++ib) {
      //if (start_index + nsteps == nmerges_occurred) break;

      // get the boundary to merge
      int boundary_to_merge = boundary_merge_order[ib];
      NeuronBoundary *boundary = nd->Boundary(boundary_to_merge);

      // get the cellular indices
      int cellular_one_data_index = boundary->SupervoxelOne()->DataIndex();
      int cellular_two_data_index = boundary->SupervoxelTwo()->DataIndex();
      rn_assertion(cellular_one_data_index < nd->NCellulars());
      rn_assertion(cellular_two_data_index < nd->NCellulars());

      // see if a merge actually occurred
      int cellular_one_label = cellular_mapping[cellular_one_data_index];
      int cellular_two_label = cellular_mapping[cellular_two_data_index];
      if (cellular_one_label == cellular_two_label) continue;


      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (cellular_mapping[ic] == cellular_two_label)
            cellular_mapping[ic] = cellular_one_label;
      }

      if (start_index <= nmerges_occurred) {
         // deflate the cellular mapping
         RNDeflateIntegerArray(cellular_mapping, ncellulars);

         // get the voxel proposals
         int nvoxels = nd->NVoxels();
         int *voxel_proposals = new int[nvoxels];
         for (int iv = 0; iv < nvoxels; ++iv) {
            int voxel_supervoxel = nd->VoxelMapping(iv);
            if (voxel_supervoxel >= nd->NCellulars()) voxel_proposals[iv] = 0;
            else voxel_proposals[iv] = cellular_mapping[voxel_supervoxel] + 1;
         }

         RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
         RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];

         // run the test metric
         if (!nd->TestMetric(voxel_proposals, test_metrics, TRUE)) return 0;
         if (!nd->SegmentationMetric(cellular_mapping, segmentation_metrics)) return 0;

         // get output filename
         char output_filename[4096];
         sprintf(output_filename, "%s/%s_%s_%06d.hier", tmp_hierarchical_directory, root_filename, extension, nmerges_occurred);

         // open file
         FILE *output_fp = fopen(output_filename, "wb");
         if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

         // write the number of entries
         int ncellulars = nd->NCellulars();
         fwrite(&ncellulars, sizeof(int), 1, output_fp);

         // write the results to file
         fwrite(test_metrics, sizeof(RNScalar), nd->NTestMetrics(), output_fp);
         fwrite(segmentation_metrics, sizeof(RNScalar), nd->NSegmentationMetrics(), output_fp);

         // close file
         fclose(output_fp);

         // free memory
         delete[] voxel_proposals;
         delete[] test_metrics;
         delete[] segmentation_metrics;
      }

      nmerges_occurred++;
   }

   // free memory
   delete[] cellular_mapping;
   delete[] boundary_merge_order;

   // return success
   return 1;
}



static int
MergeHierarchicalResults(const char root_filename[4096])
{   
   // create the images
   if (!CreateMetricImages(root_filename)) return 0;

   // if a file cannot be opened make sure none after can
   RNBoolean failure_check = FALSE;
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      char metrics_filename[4096];
      sprintf(metrics_filename, "%s/%s_%s_%06d.hier", tmp_hierarchical_directory, root_filename, extension, ic);

      // open file
      FILE *metric_fp = fopen(metrics_filename, "rb");
      if (!metric_fp) {
         fprintf(stderr, "Failed to read %s...trying to recover...\n", metrics_filename);
         failure_check = TRUE;
         continue;
      }

      if (failure_check) {
         fprintf(stderr, "Recovery failed...aborting...\n"); exit(-1);
      }

      // read the number of entries
      int ncellulars;
      fread(&ncellulars, sizeof(int), 1, metric_fp);
      rn_assertion(ncellulars == nd->NCellulars());

      // read in the metrics
      for (int im = 0; im < nd->NTestMetrics(); ++im) {
         RNScalar metric_result;
         fread(&metric_result, sizeof(RNScalar), 1, metric_fp);
         test_plots[im]->InsertPoint(R2Point(ic, metric_result));
      }
      for (int im = 0; im < nd->NSegmentationMetrics(); ++im) {
         RNScalar metric_result;
         fread(&metric_result, sizeof(RNScalar), 1, metric_fp);
         segmentation_plots[im]->InsertPoint(R2Point(ic, metric_result));
      }

      // close file
      fclose(metric_fp);
   }

   // write the images
   if (!WriteMetricImages(root_filename)) return 0;

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
   // Parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-start_index")) { argv++; argc--; start_index = atoi(*argv); }
         else if (!strcmp(*argv, "-merge")) merge_results = 1;
         else if (!strcmp(*argv, "-nsteps")) { argv++; argc--; nsteps = atoi(*argv); }
         else if (!strcmp(*argv, "-random_forest")) extension = "random_forest_smoothing_global";
         else if (!strcmp(*argv, "-boundary_max")) extension = "boundary_max";
         else if (!strcmp(*argv, "-boundary_mean")) extension = "boundary_mean";
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input argument
   if (!input_filename) { fprintf(stderr, "Need to supply input filename\n"); return 0; }

   if (!merge_results && (start_index == -1 || nsteps == -1)) { fprintf(stderr, "Either need to merge existing results or choose a starting index and nsteps to complete\n"); return 0; }

   if (!extension) { fprintf(stderr, "Need to supply ordering extension (random_forest, boundary_max, boundary_mean)\n"); return 0; }

   // Return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in the neuron data
   if (!ReadData(input_filename)) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // either run the subsection or merge results
   if (!merge_results && !HierarchicalSubsection(root_filename)) exit(-1);
   if (merge_results && !MergeHierarchicalResults(root_filename)) exit(-1);

   // free memory
   delete nd;

   // return success
   return 0;
}
