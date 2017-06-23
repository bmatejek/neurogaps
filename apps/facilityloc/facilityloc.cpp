// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include "RNML/RNML.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int maximum_facilities = -1;
static int create_feature = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static float **dijkstra_distances = NULL;



// directory structure

static const char *features_directory = "algs_data/features";
static const char *distances_directory = "distances";
static const char *results_directory = "results/facilities";
static const char *visuals_directory = "visuals/facilities";
static const char *boundaries_directory = "boundaries";



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



static int CreateFeature(char root_filename[4096])
{
   if (maximum_facilities == -1) { fprintf(stderr, "Need to supply maximum medians\n"); return 0; }

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
   for (int K = 2; K < maximum_facilities; ++K) {
      RNProgressBar(K - 2, maximum_facilities - 2);
      // get boundary filename
      char boundaries_filename[4096];
      sprintf(boundaries_filename, "%s/%s_facility_%04d.boundary", boundaries_directory, root_filename, K);

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
            if (first_split[ib] == -1) first_split[ib] = K;
	    last_split[ib] = K;
	    proportion_split[ib] += 1.0 / maximum_facilities;
         }
      }

      // free memory
      delete[] boundary_results;
   }
   printf("\n");

   // find variables that are not split ever
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      if (first_split[ib] == -1) first_split[ib] = maximum_facilities + 1;
      if (last_split[ib] == -1) last_split[ib] = maximum_facilities + 1;
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
   sprintf(first_split_filename, "%s/%s_facility_location_first_split.feature", features_directory, root_filename);
   char last_split_filename[4096];
   sprintf(last_split_filename, "%s/%s_facility_location_last_split.feature", features_directory, root_filename);
   char proportion_filename[4096];
   sprintf(proportion_filename, "%s/%s_facility_location_proportion.feature", features_directory, root_filename);

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
      sprintf(title, "%s Facility Location %s", root_filename, nd->TestMetricName(im));
      test_plots[im]->SetTitle(title);
      test_plots[im]->SetXLabel("Number of Facilities on the Boundary");
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
      sprintf(title, "%s Facility Location %s", root_filename, nd->SegmentationMetricName(im));
      sprintf(ylabel, "%s", nd->SegmentationMetricName(im));

      segmentation_plots[im]->SetTitle(title);
      segmentation_plots[im]->SetXLabel("Number of Facilities on the Boundary");
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
// Farthest-first traversal
////////////////////////////////////////////////////////////////////////

static int *FarthestFirstTraversal(int K)
{
   // choose K centers from the boundaries
   int *facility_locations = new int[K];

   // get a list of boundary cellulars
   std::vector<int> possible_facilities = std::vector<int>();
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      if (cellular->IsOnBoundary()) possible_facilities.push_back(ic);
   }

   // find facility initially farthest from all other
   float farthest_distance = 0.0;
   float farthest_index = -1;
   for (unsigned int ic1 = 0; ic1 < possible_facilities.size(); ++ic1) {
      NeuronCellular *cellular_one = nd->Cellular(possible_facilities[ic1]);

      // find the closest facility to this cellular
      float nearest_location = FLT_MAX;

      for (unsigned int ic2 = 0; ic2 < possible_facilities.size(); ++ic2) {
         if (ic1 == ic2) continue;
         NeuronCellular *cellular_two = nd->Cellular(possible_facilities[ic2]);

         float distance = dijkstra_distances[cellular_one->DataIndex()][cellular_two->DataIndex()];
         if (distance < nearest_location)
            nearest_location = distance;
      }

      if (nearest_location > farthest_distance) {
         farthest_distance = nearest_location;
         farthest_index = ic1;
      }
   }

   // randomly choose a facility location for the first
   facility_locations[0] = possible_facilities[farthest_index];
   
   // choose the furthest points remaining
   for (int ik = 1; ik < K; ++ik) {
      int farthest_facility = -1;
      float farthest_distance = 0.0;

      // go through all possible facility locations
      for (unsigned int ip = 0; ip < possible_facilities.size(); ++ip) {
         float minimum_distance = FLT_MAX;
         // go through all current facilities
         for (int ic = 0; ic < ik; ++ic) {
            float distance = dijkstra_distances[possible_facilities[ip]][facility_locations[ic]];
            if (distance < minimum_distance) minimum_distance = distance;
         }

         // see if this is the new farthest facility
         if (minimum_distance > farthest_distance) {
            farthest_distance = minimum_distance;
            farthest_facility = ip;
         }
      }

      // add the farthest facility to the list
      rn_assertion(farthest_facility != -1);
      facility_locations[ik] = possible_facilities[farthest_facility];
   }

   // return the end facility locations
   return facility_locations;
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
         else if (!strcmp(*argv, "-maximum_facilities")) { argv++; argc--; maximum_facilities = atoi(*argv); }
	 else if (!strcmp(*argv, "-create_feature")) create_feature = 1;
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
   else if (!create_feature) {
      // read voxels
      nd->ReadVoxels();

      // read dijkstra distances
      if (!ReadDijkstraDistances(root_filename)) exit(-1);
      
      // count the number of boundary cellulars
      int nboundary_cellulars = 0;
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
	 NeuronCellular *cellular = nd->Cellular(ic);
	 if (cellular->IsOnBoundary()) nboundary_cellulars++;
      }
      
      // unless the user specified the maximum number of facilities
      if (maximum_facilities != -1) nboundary_cellulars = maximum_facilities + 1;
      
      // create the structure for all of the metric images
      if (!CreateMetricImages(root_filename)) exit(-1);
      
      // run multiple iterations of the facility location problem
      for (int K = 2; K < nboundary_cellulars; ++K) {
	 // start statistics
	 if (print_verbose) { printf("Generating results for %d facilities...", K); fflush(stdout); }
	 RNTime location_time;
	 location_time.Read();
	 
	 // get the facilities for this k
	 int *facility_locations = FarthestFirstTraversal(K);
	 
	 // assign all cellulars to the nearest facility
	 int *cellular_labeling = new int[nd->NCellulars()];
	 for (int ic = 0; ic < nd->NCellulars(); ++ic) {
	    NeuronCellular *cellular = nd->Cellular(ic);
	    int cellular_data_index = cellular->DataIndex();
	    
	    // get the minimum facility distance
	    float minimum_distance = FLT_MAX;
	    int minimum_facility_location = -1;
	    
	    // go through all facilities
	    for (int ik = 0; ik < K; ++ik) {
	       int facility_data_index = facility_locations[ik];
	       float distance = dijkstra_distances[cellular_data_index][facility_data_index];
	       if (distance < minimum_distance) {
		  minimum_distance = distance;
		  minimum_facility_location = ik;
	       }
	    }
	    
	    rn_assertion(minimum_facility_location != -1);
	    cellular_labeling[ic] = minimum_facility_location;
	 }
	 
	 // save the boundary information
	 char boundaries_filename[4096];
	 sprintf(boundaries_filename, "%s/%s_facility_%04d.boundary", boundaries_directory, root_filename, K);
	 
	 // open file
	 FILE *boundaries_fp = fopen(boundaries_filename, "wb");
	 if (!boundaries_fp) { fprintf(stderr, "Failed to open %s\n", boundaries_filename); exit(-1); }
	 
	 int nboundaries = nd->NBoundaries();
	 fwrite(&nboundaries, sizeof(int), 1, boundaries_fp);

	 // find results of the facility location algorithm
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
	       if (cellular_labeling[cellular_one_data_index] == cellular_labeling[cellular_two_data_index])
		  boundaries[ib] = TRUE;
	       else boundaries[ib] = FALSE;
	    }
	 }
	 
	 // write all of the boundaries
	 fwrite(boundaries, sizeof(RNBoolean), nboundaries, boundaries_fp);
	 
	 // free memory
	 delete[] boundaries;
      
	 // close file
	 fclose(boundaries_fp);

	 // create proposal array for voxels
	 int *proposals = new int[nd->NVoxels()];
	 for (int iv = 0; iv < nd->NVoxels(); ++iv) {
	    NeuronVoxel *voxel = nd->Voxel(iv);

	    // cellular voxels are labeled according to facility location problem
	    if (voxel->IsCellular()) {
	       NeuronCellular *cellular = (NeuronCellular *)voxel->Supervoxel();
	       int cellular_data_index = cellular->DataIndex();
	       proposals[iv] = cellular_labeling[cellular_data_index] + 1;
	    }
	    // otherwise labeled extracellular
	    else proposals[iv] = 0;
	 }
	 
	 // run test metrics
	 RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
	 nd->TestMetric(proposals, test_metrics, FALSE);

	 // add these metric points to the plots
	 for (int im = 0; im < nd->NTestMetrics(); ++im) {
	    test_plots[im]->InsertPoint(R2Point(K, test_metrics[im]));
	 }
	 
	 // run supervoxel segmentation metrics
	 RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];
	 nd->SegmentationMetric(cellular_labeling, segmentation_metrics);
	 
	 // add these metric points to the plots
	 for (int im = 0; im < nd->NSegmentationMetrics(); ++im) {
	    segmentation_plots[im]->InsertPoint(R2Point(K, segmentation_metrics[im]));
	 }

	 // free memory
	 delete[] test_metrics;
	 delete[] segmentation_metrics;
	 delete[] proposals;
	 delete[] facility_locations;
	 delete[] cellular_labeling;

	 // print statistics
	 if (print_verbose) printf("done in %0.2f seconds\n", location_time.Elapsed());
      }
      
      // write the images
      if (!WriteMetricImages(root_filename)) exit(-1);
      
      for (int ic = 0; ic < nd->NCellulars(); ++ic)
	 delete[] dijkstra_distances[ic];
      delete[] dijkstra_distances;
   }
   
   // free up memory
   delete nd;

   // return success
   return 0;
}
