// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static float *distances = NULL;
static int *prev_node = NULL;

// directory structure

static const char *distances_directory = "distances";
static const char *results_directory = "results/dijkstra";
static const char *visuals_directory = "visuals/dijkstra";



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
   if (!nd->ReadFile(filename, TRUE)) {
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



static int ReadDijkstraFiles(const char *root_filename)
{
   // get distances filename
   char distances_filename[4096];
   sprintf(distances_filename, "%s/%s_dijkstra_distances_boundary.distance", distances_directory, root_filename);

   // open file
   FILE *distances_fp = fopen(distances_filename, "rb");
   if (!distances_fp) { fprintf(stderr, "Failed to read %s\n", distances_filename); return 0; }

   int distances_nvoxels;
   fread(&distances_nvoxels, sizeof(int), 1, distances_fp);
   rn_assertion(distances_nvoxels == nd->NVoxels());
   distances = new float[nd->NVoxels()];
   if (!distances) { fprintf(stderr, "Failed to allocate memory for distances array\n"); return 0; }
   fread(distances, sizeof(float), nd->NVoxels(), distances_fp);

   // close file
   fclose(distances_fp);

   // get prev node filename
   char prev_node_filename[4096];
   sprintf(prev_node_filename, "%s/%s_dijkstra_prev_node_boundary.distance", distances_directory, root_filename);

   // open file
   FILE *prev_node_fp = fopen(prev_node_filename, "rb");
   if (!prev_node_fp) { fprintf(stderr, "Failed to read %s\n", prev_node_filename); return 0; }

   int prev_node_nvoxels;
   fread(&prev_node_nvoxels, sizeof(int), 1, prev_node_fp);
   rn_assertion(prev_node_nvoxels == nd->NVoxels());
   prev_node = new int[nd->NVoxels()];
   if (!prev_node) { fprintf(stderr, "Failed to allocate memory for prev node array\n"); return 0; }
   fread(prev_node, sizeof(int), nd->NVoxels(), prev_node_fp);

   // close file
   fclose(prev_node_fp);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Path generation functions
////////////////////////////////////////////////////////////////////////

// store useful information for dijkstra path
struct DijkstraPath {
   // constructor
   DijkstraPath(RNBoolean in_ground_truth, RNScalar maximum_step, RNScalar distance, int ncorrect_steps, int nincorrect_steps, RNBoolean correct_path) :
      in_ground_truth(in_ground_truth),
      maximum_step(maximum_step),
      distance(distance),
      ncorrect_steps(ncorrect_steps),
      nincorrect_steps(nincorrect_steps),
      correct_path(correct_path)
   {}

   // instance variables
   RNBoolean in_ground_truth;
   RNScalar maximum_step;
   RNScalar distance;
   int ncorrect_steps;
   int nincorrect_steps;
   RNBoolean correct_path;
};



// create a vector for all of the dijkstra path attributes
static std::vector<DijkstraPath> voxel_paths = std::vector<DijkstraPath>();



// physically generate all of the dijkstra paths
static int GenerateDijkstraPaths(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // go through all voxels
   if (print_verbose) { printf("Analyzing statistics for all voxels...\n  "); fflush(stdout); }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (print_verbose) RNProgressBar(iv, nd->NVoxels());

      // consider only voxels in ground truth     
      NeuronVoxel *voxel = nd->Voxel(iv);

      RNScalar distance = distances[iv];
      int prev_node_index = iv;

      // get this human label
      NeuronHumanLabel *human_label = voxel->HumanLabel();
      RNBoolean correct_path = TRUE;

      // largest step along the path
      RNScalar maximum_step = 0.0;

      // keep track of the number of invidual correct steps
      int nindividual_correct_steps = 0;
      int nindividual_incorrect_steps = 0;

      // iterate through the path
      while (prev_node_index != -1) {
         NeuronVoxel *current_voxel = nd->Voxel(prev_node_index);
         if (current_voxel->HumanLabel() != human_label) {
            correct_path = FALSE;
            nindividual_incorrect_steps++;
         }
         else nindividual_correct_steps++;

         // update the previous node
         prev_node_index = prev_node[prev_node_index];

         // see if this step along the path was the greatest
         if (prev_node_index != -1) {
            NeuronVoxel *previous_voxel = nd->Voxel(prev_node_index);
            RNScalar affinity = 1.0 - current_voxel->AffinityToNeighbor(previous_voxel);
            if (affinity > maximum_step) maximum_step = affinity;
         }
      }

      // add this to voxel path list
      voxel_paths.push_back(DijkstraPath(human_label != NULL, maximum_step, distance, nindividual_correct_steps, nindividual_incorrect_steps, correct_path));
   }
   if (print_verbose) printf("\ndone in %0.2f seconds.\n\n", start_time.Elapsed());

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Generate statistics
////////////////////////////////////////////////////////////////////////

static int GenerateGroundTruthPathStatistics(const char root_filename[4096])
{
   // counter variables
   int ncorrect_paths = 0;
   int nincorrect_paths = 0;
   unsigned long long ncorrect_steps = 0;
   unsigned long long nincorrect_steps = 0;

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      DijkstraPath path = voxel_paths[iv];

      // only consider ground truth
      if (!path.in_ground_truth) continue;

      // is this path correct
      if (path.correct_path) ncorrect_paths++;
      else nincorrect_paths++;

      // how many steps along this path are incorrect
      ncorrect_steps += path.ncorrect_steps;
      nincorrect_steps += path.nincorrect_steps;
   }

   // what percent of paths were correct
   printf("Correct Paths: %d\n", ncorrect_paths);
   printf("Incorrect Paths: %d\n", nincorrect_paths);
   printf("Proportion: %lf\n\n", ncorrect_paths / ((RNScalar)ncorrect_paths + nincorrect_paths));
   printf("Correct Steps: %llu\n", ncorrect_steps);
   printf("Incorrect Steps: %llu\n", nincorrect_steps);
   printf("Proportion: %lf\n\n", ncorrect_steps / ((RNScalar)ncorrect_steps + nincorrect_steps));
   printf("Incorrect Steps On Average: %lf\n", nincorrect_steps / (RNScalar)voxel_paths.size());

   // create output charts
   RNPyHistogram *error_histogram = new RNPyHistogram();
   error_histogram->SetTitle("Path Errors over all Dijkstra Paths");
   error_histogram->SetYLabel("Number of Occurences");
   error_histogram->SetXLabel("Number of Path Errors");

   // set the maximum number of errors as bin size
   int maximum_errors = 0;
   for (unsigned int iv = 0; iv < voxel_paths.size(); ++iv) {
      DijkstraPath path = voxel_paths[iv];

      // only consider voxels that are in the ground truth
      if (!path.in_ground_truth) continue;

      // add the point
      if (path.nincorrect_steps > maximum_errors) maximum_errors = path.nincorrect_steps;
      error_histogram->InsertPoint(path.nincorrect_steps);
   }
   error_histogram->SetNBins(maximum_errors);
   error_histogram->SetLegend("# Voxel Errors");

   // create the output image
   RNPyImage *histograms = new RNPyImage();
   histograms->SetOutputDirectory(visuals_directory);
   histograms->InsertPyPlot(error_histogram);

   // output filename
   char histogram_output_filename[4096];
   sprintf(histogram_output_filename, "%s/%s_nerrors.pyimage", results_directory, root_filename);

   // write the image file
   histograms->WriteImageFile(histogram_output_filename);

   // free memory
   delete error_histogram;
   delete histograms;

   // return success
   return 1;
}



static int GenerateCellularPathStatistics(const char root_filename[4096])
{
   RNPyHistogram *vessel_histogram = new RNPyHistogram();
   RNPyHistogram *nonvessel_histogram = new RNPyHistogram();
   vessel_histogram->SetTitle("Maximum Dijkstra Step for Voxels in Ground Truth");
   vessel_histogram->SetYLabel("Number of Occurences");
   vessel_histogram->SetXLabel("Maximum Dijkstra Step");
   nonvessel_histogram->SetTitle("Maximum Dijkstra Step for Voxels Not in Ground Truth");
   nonvessel_histogram->SetYLabel("Number of Occurences");
   nonvessel_histogram->SetXLabel("Maximum Dijkstra Step");

   // set the maximum number of errors as bin size
   for (unsigned int iv = 0; iv < voxel_paths.size(); ++iv) {
      // only consider the voxels that the affinity map considers cellular
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->IsExtracellular()) continue;

      // get the path variable
      DijkstraPath path = voxel_paths[iv];

      if (voxel->HumanLabel()) vessel_histogram->InsertPoint(path.maximum_step);
      else nonvessel_histogram->InsertPoint(path.maximum_step);
   }
   vessel_histogram->SetNBins(100);
   nonvessel_histogram->SetNBins(100);
   vessel_histogram->SetLegend("Maximum Dijkstra Step");
   nonvessel_histogram->SetLegend("Maximum Dijkstra Step");

   RNPyPointPlot *proportion_plot = new RNPyPointPlot();
   proportion_plot->SetTitle("Maximum Step Proportions");
   proportion_plot->SetXLabel("Maximum Step");
   proportion_plot->SetYLabel("Proportion of Paths For This Maximum");
   proportion_plot->SetLegend("Proportion Not in Ground Truth");

   int nbins = 1000;
   unsigned int *ncorrect = new unsigned int[nbins];
   unsigned int *nincorrect = new unsigned int[nbins];
   for (int ib = 0; ib < nbins; ++ib) {
      ncorrect[ib] = 0;
      nincorrect[ib] = 0;
   }

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      // only consider voxels that affinity map thinks are vessels
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->IsExtracellular()) continue;

      DijkstraPath path = voxel_paths[iv];
      RNScalar maximum_step = path.maximum_step;
      int bin = (int)(maximum_step * nbins + 0.5);

      // increment if voxel is in vessel or not
      if (path.in_ground_truth) ncorrect[bin]++;
      else nincorrect[bin]++;
   }

   // insert plot points
   for (int ib = 0; ib < nbins; ++ib) {
      proportion_plot->InsertPoint(R2Point(ib / 1000.0, nincorrect[ib] / (RNScalar)(ncorrect[ib] + nincorrect[ib])));
   }

   // free memory
   delete[] ncorrect;
   delete[] nincorrect;
   
   // create the output image
   RNPyImage *histograms = new RNPyImage();
   histograms->SetOutputDirectory(visuals_directory);
   histograms->InsertPyPlot(vessel_histogram);
   histograms->InsertPyPlot(nonvessel_histogram);

   // output filename
   char histograms_output_filename[4096];
   sprintf(histograms_output_filename, "%s/%s_maximum_step.pyimage", results_directory, root_filename);

   // write the image file
   histograms->WriteImageFile(histograms_output_filename);

   // create the output image
   RNPyImage *plot = new RNPyImage();
   plot->SetOutputDirectory(visuals_directory);
   plot->InsertPyPlot(proportion_plot);

   // output filename
   char plot_output_filename[4096];
   sprintf(plot_output_filename, "%s/%s_maximum_step_proportion.pyimage", results_directory, root_filename);

   // write the image file
   plot->WriteImageFile(plot_output_filename);

   // free memory
   delete vessel_histogram;
   delete nonvessel_histogram;
   delete proportion_plot;
   delete histograms;
   delete plot;

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

   // generate all the statistics
   if (!ReadDijkstraFiles(root_filename)) exit(-1);
   if (!GenerateDijkstraPaths()) exit(-1);
   if (!GenerateGroundTruthPathStatistics(root_filename)) exit(-1);
   if (!GenerateCellularPathStatistics(root_filename)) exit(-1);

   // free up memory
   delete nd;
   delete[] distances;
   delete[] prev_node;

   // return success
   return 0;
}
