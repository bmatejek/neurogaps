// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>
#include "RNML/RNML.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int maximum_medians = -1;
static int create_feature = 0;
static int starting_K = -1;
static int nsteps = 20;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static float **dijkstra_distances = NULL;



// directory structure

static const char *distances_directory = "distances";
static const char *results_directory = "results/kmedians";
static const char *visuals_directory = "visuals/kmedians";
static const char *boundaries_directory = "boundaries";
static const char *features_directory = "algs_data/features";



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



static int ReadDijkstraDistances(const char root_filename[4096])
{
   // get the dijkstra distances filename
   char dijkstra_filename[4096];
   sprintf(dijkstra_filename, "%s/%s_dijkstra.distance", distances_directory, root_filename);

   // open the distance file
   FILE *dijkstra_fp = fopen(dijkstra_filename, "rb");
   if (!dijkstra_fp) { fprintf(stderr, "Failed to read %s\n", dijkstra_filename); return 0; }

   // read the number of cellulars
   int ndijkstra_cellulars;
   fread(&ndijkstra_cellulars, sizeof(int), 1, dijkstra_fp);
   rn_assertion(ndijkstra_cellulars == nd->NCellulars());

   // read in dijkstra distances
   dijkstra_distances = new float *[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      dijkstra_distances[ic] = new float[nd->NCellulars()];
      if (fread(dijkstra_distances[ic], sizeof(float), nd->NCellulars(), dijkstra_fp) != (unsigned int)nd->NCellulars()) {
         fprintf(stderr, "Failed to read dijkstra distance file: %s\n", dijkstra_filename);
         return 0;
      }
   }

   // close file
   fclose(dijkstra_fp);

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
      sprintf(title, "%s K-Median %s", root_filename, nd->TestMetricName(im));
      test_plots[im]->SetTitle(title);
      test_plots[im]->SetXLabel("Number of Medians");
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
      sprintf(title, "%s K-Medians %s", root_filename, nd->SegmentationMetricName(im));
      sprintf(ylabel, "%s", nd->SegmentationMetricName(im));

      segmentation_plots[im]->SetTitle(title);
      segmentation_plots[im]->SetXLabel("Number of Medians");
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
// E/M functions
////////////////////////////////////////////////////////////////////////

static void ExpectationStep(int K, int *clusters, int *medians)
{
   // just checking
   rn_assertion(clusters != NULL);
   rn_assertion(medians != NULL);

   // find the new median for each cluster
   for (int ik = 0; ik < K; ++ik) {
      // get all of the cellulars in this cluster
      std::vector<int> cluster_cellulars = std::vector<int>();
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (clusters[ic] == ik) {
            cluster_cellulars.push_back(ic);
         }
      }

      // find the center with the minimum distance
      int best_median = -1;
      float minimum_distance = FLT_MAX;

      // go through every cellular in this cluster
      for (unsigned int ic1 = 0; ic1 < cluster_cellulars.size(); ++ic1) {
         float distance = 0.0;

         // get this cellular information
         NeuronCellular *cellular_one = nd->Cellular(cluster_cellulars[ic1]);
         int cellular_one_data_index = cellular_one->DataIndex();

         // go through every other cellular in this cluster
         for (unsigned int ic2 = 0; ic2 < cluster_cellulars.size(); ++ic2) {

            // get the other cellular information
            NeuronCellular *cellular_two = nd->Cellular(cluster_cellulars[ic2]);
            int cellular_two_data_index = cellular_two->DataIndex();

            // how far is cellular one from cellular two
            distance += dijkstra_distances[cellular_one_data_index][cellular_two_data_index];
         }

         // if a new minimum is found, update information
         if (distance < minimum_distance) {
            minimum_distance = distance;
            best_median = ic1;
         }
      }

      // update the median for cluster ik
      if (best_median == -1) medians[ik] = -1;
      else medians[ik] = cluster_cellulars[best_median];
   }
}



static RNScalar MaximizationStep(int K, int *clusters, int *medians)
{
   // just checking
   rn_assertion(clusters != NULL);
   rn_assertion(medians != NULL);

   // the total distance
   RNScalar median_distance = 0.0;

   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      // find the closest median
      float minimum_median = FLT_MAX;

      for (int ik = 0; ik < K; ++ik) {
         if (medians[ik] == -1) continue;
         rn_assertion(medians[ik] < nd->NCellulars());
         float distance = dijkstra_distances[ic][medians[ik]];
         if (distance < minimum_median) {
            minimum_median = distance;
            clusters[ic] = ik;
         }
      }

      median_distance += minimum_median;
   }

   // return the distance
   return median_distance;
}




static int RunKMedians(char root_filename[4096])
{
   // read voxels
   nd->ReadVoxels();

   // read dijkstra distances
   if (!ReadDijkstraDistances(root_filename)) exit(-1);

   // create a vector from [0...ncellulars - 1]
   std::vector<int> cellular_labels = std::vector<int>();
   for (int ic = 0; ic < nd->NCellulars(); ++ic)
      cellular_labels.push_back(ic);

   for (int K = starting_K; K < starting_K + nsteps; ++K) {
      if (print_verbose) { printf("Generating results for %d medians...", K); fflush(stdout); }

      // start statistics
      RNTime median_time;
      median_time.Read();

      // array mapping cellulars to clusters, initially 0
      int *clusters = new int[nd->NCellulars()];

      int *best_medians = new int[K];
      RNScalar best_cost = FLT_MAX;

      const int niterations = 100;
      for (int i = 0; i < niterations; ++i) {
         // initialize all clusters to -1
         for (int ic = 0; ic < nd->NCellulars(); ++ic)
            clusters[ic] = -1;

         // randomly choose starting median locations
         std::random_shuffle(cellular_labels.begin(), cellular_labels.end());

         // the initial starting location for the centers
         int *medians = new int[K];
         for (int ik = 0; ik < K; ++ik)
            medians[ik] = cellular_labels[ik];

         RNScalar cost = MaximizationStep(K, clusters, medians);
         RNScalar prev_cost;

         do {
            // update the previous cost
            prev_cost = cost;

            // find the next best medians
            ExpectationStep(K, clusters, medians);

            // update the cost
            cost = MaximizationStep(K, clusters, medians);
         } while (prev_cost - cost > 10e-8);

         // update if best cost so far
         if (cost < best_cost) {
            best_cost = cost;
            for (int ik = 0; ik < K; ++ik) {
               best_medians[ik] = medians[ik];
            }
         }

         delete[] medians;
      }

      // update clusters with the best medians
      RNScalar tmp_cost = MaximizationStep(K, clusters, best_medians);
      rn_assertion(abs(tmp_cost - best_cost) < 10e-6);

      // save the boundary information
      char boundaries_filename[4096];
      sprintf(boundaries_filename, "%s/%s_kmedian_%04d.boundary", boundaries_directory, root_filename, K);

      // open file
      FILE *boundaries_fp = fopen(boundaries_filename, "wb");
      if (!boundaries_fp) { fprintf(stderr, "Failed to open %s\n", boundaries_filename); return 0; }

      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, boundaries_fp);

      // find results of k-median algorithm
      RNBoolean *boundaries = new RNBoolean[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);

         // get supervoxels
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

         // extracellulars are never merged
         if (supervoxel_one->IsExtracellular() || supervoxel_two->IsExtracellular()) {
            boundaries[ib] = FALSE;
         }
         else {
            int cellular_one_data_index = supervoxel_one->DataIndex();
            int cellular_two_data_index = supervoxel_two->DataIndex();

            // are these cellulars merged?
            if (clusters[cellular_one_data_index] == clusters[cellular_two_data_index])
               boundaries[ib] = TRUE;
            else boundaries[ib] = FALSE;
         }
      }

      // write all of the boundaries
      fwrite(boundaries, sizeof(RNBoolean), nboundaries, boundaries_fp);

      // free memory
      delete[] boundaries;

      // create the proposal arrays
      int *proposals = new int[nd->NVoxels()];
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         NeuronVoxel *voxel = nd->Voxel(iv);
         // find the cluster this label belongs to
         if (voxel->IsCellular()) {
            NeuronCellular *cellular = (NeuronCellular *)voxel->Supervoxel();
            int cellular_data_index = cellular->DataIndex();
            proposals[iv] = clusters[cellular_data_index] + 1;
         }
         // otherwise label it as extracellular
         else proposals[iv] = 0;
      }

      // deflate proposals
      RNDeflateIntegerArray(proposals, nd->NVoxels());

      // run test metric
      RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
      nd->TestMetric(proposals, test_metrics, FALSE);

      // write these metric values to the boundary file
      fwrite(test_metrics, sizeof(RNScalar), nd->NTestMetrics(), boundaries_fp);

      // run supervoxel semgentation metrics
      RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];
      nd->SegmentationMetric(clusters, segmentation_metrics);

      // write the metric values to the boundary file
      fwrite(segmentation_metrics, sizeof(RNScalar), nd->NSegmentationMetrics(), boundaries_fp);

      // close file
      fclose(boundaries_fp);

      // free memory
      delete[] test_metrics;
      delete[] segmentation_metrics;
      delete[] proposals;
      delete[] clusters;
      delete[] best_medians;

      // print statistics
      if (print_verbose) printf("done in %0.2f seconds\n", median_time.Elapsed());
   }

   // free memory
   for (int ic = 0; ic < nd->NCellulars(); ++ic)
      delete[] dijkstra_distances[ic];
   delete[] dijkstra_distances;

   // return success
   return 1;
}



static int CreateFeature(char root_filename[4096])
{
   if (maximum_medians == -1) { fprintf(stderr, "Need to supply maximum medians\n"); return 0; }

   // create the structure for all of the metric images
   if (!CreateMetricImages(root_filename)) exit(-1);

   // keep track of when boundaries are split
   int *first_split = new int[nd->NBoundaries()];
   int *last_split = new int[nd->NBoundaries()];
   float *proportion_split = new float[nd->NBoundaries()];
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      first_split[ib] = -1;
      last_split[ib] = -1;
      proportion_split[ib] = 0.0;
   }

   // go through all k median results
   for (int K = 2; K < maximum_medians; ++K) {
      RNProgressBar(K - 2, maximum_medians - 2);
      // get boundary filename
      char boundaries_filename[4096];
      sprintf(boundaries_filename, "%s/%s_kmedian_%04d.boundary", boundaries_directory, root_filename, K);

      // open file
      FILE *boundary_fp = fopen(boundaries_filename, "rb");
      if (!boundary_fp) { fprintf(stderr, "Failed to read %s\n", boundaries_filename); return 0; }

      // read the number of boundaries
      int input_nboundaries;
      fread(&input_nboundaries, sizeof(int), 1, boundary_fp);
      rn_assertion(input_nboundaries == nd->NBoundaries());

      RNBoolean *boundary_results = new RNBoolean[nd->NBoundaries()];
      if (fread(boundary_results, sizeof(RNBoolean), nd->NBoundaries(), boundary_fp) != (unsigned int)nd->NBoundaries()) {
         fprintf(stderr, "Failed to read boundary results from %s\n", boundaries_filename);
         return 0;
      }

      // read in test and segmentation scores
      RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
      RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];
      if (fread(test_metrics, sizeof(RNScalar), nd->NTestMetrics(), boundary_fp) != (unsigned int)nd->NTestMetrics()) {
         fprintf(stderr, "Failed to read test metric results from %s\n", boundaries_filename);
         return 0;
      }
      if (fread(segmentation_metrics, sizeof(RNScalar), nd->NSegmentationMetrics(), boundary_fp) != (unsigned int)nd->NSegmentationMetrics()) {
         fprintf(stderr, "Failed to read segmentation metric results from %s\n", boundaries_filename);
         return 0;
      }

      for (int it = 0; it < nd->NTestMetrics(); ++it)
         test_plots[it]->InsertPoint(R2Point(K, test_metrics[it]));
      for (int is = 0; is < nd->NSegmentationMetrics(); ++is)
         segmentation_plots[is]->InsertPoint(R2Point(K, segmentation_metrics[is]));

      // close file
      fclose(boundary_fp);

      // go through all boundaries
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         if (!boundary_results[ib]) {
            if (first_split[ib] == -1) first_split[ib] = K;
            last_split[ib] = K;
            proportion_split[ib] += 1.0 / maximum_medians;
         }
         else {
            if (last_split[ib] != -1) last_split[ib] = -1;
         }
      }

      // free memory
      delete[] boundary_results;
   }
   printf("\n");

   // find variables that are not split ever
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      if (first_split[ib] == -1) first_split[ib] = maximum_medians + 1;
      if (last_split[ib] == -1) last_split[ib] = maximum_medians + 1;
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

   // write the images
   if (!WriteMetricImages(root_filename)) return 0;

   // create output feature files
   char first_split_filename[4096];
   sprintf(first_split_filename, "%s/%s_kmedian_first_split.feature", features_directory, root_filename);
   char last_split_filename[4096];
   sprintf(last_split_filename, "%s/%s_kmedian_last_split.feature", features_directory, root_filename);
   char proportion_filename[4096];
   sprintf(proportion_filename, "%s/%s_kmedian_proportion.feature", features_directory, root_filename);

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
         else if (!strcmp(*argv, "-maximum_medians")) { argv++; argc--; maximum_medians = atoi(*argv); }
         else if (!strcmp(*argv, "-create_feature")) create_feature = 1;
         else if (!strcmp(*argv, "-starting_K")) { argv++; argc--; starting_K = atoi(*argv); }
         else if (!strcmp(*argv, "-nsteps")) { argv++; argc--; nsteps = atoi(*argv); }
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
   if (starting_K == -1 && !create_feature) { fprintf(stderr, "Either need a starting k or to create feature.\n"); return 0; }

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

   if (create_feature && !CreateFeature(root_filename)) exit(-1);
   if (!create_feature && !RunKMedians(root_filename)) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}
