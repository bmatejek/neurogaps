// Source file for the simple file conversion algorithm



// include files

#include "RNML/RNML.h"
#include "Neuron/Neuron.h"
#include <vector>
#include <limits>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int predict_all = 0;
// string arguments
static const char *input_filename = NULL;
static const char *extension = "post";



// directory structure

static const char *dataset_directory = "algs_data/dataset";
static const char *postprocess_directory = "postprocess";



// useful constants

static const int NSTATISTICS = 26;
static const char *statistic_names[NSTATISTICS] = {
   "distances", "nsteps", "kurtosis", "maximum", "mean",
   "median", "minimum", "skew", "stddev", "first_quintile",
   "second_quintile", "third_quintile", "fourth_quintile",
   "distances_rank", "nsteps_rank", "kurtosis_rank",
   "maximum_rank", "mean_rank", "median_rank", "minimum_rank",
   "skew_rank", "stddev_rank", "first_quintile_rank",
   "second_quintile_rank", "third_quintile_rank", "fourth_quintile_rank"
};



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static NeuronData *ReadData(const char *filename)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   NeuronData *nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate memory for neuron data\n");
      return NULL;
   }

   // read in the file
   if (!nd->ReadFile(filename, TRUE)) {
      fprintf(stderr, "Failed to read file\n");
      return NULL;
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

   // return the neuron structure
   return nd;
}



////////////////////////////////////////////////////////////////////////
// Create data sets
////////////////////////////////////////////////////////////////////////

static RNDataset *CreatePostDataset(NeuronData *nd)
{
   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create the feature vectors from this neuron section
   std::vector<std::string> names = std::vector<std::string>();
   std::vector<std::vector<float> > features = std::vector<std::vector<float> >();
   std::vector<unsigned int> labels = std::vector<unsigned int>();
   std::vector<unsigned int> identifications = std::vector<unsigned int>();


   /////////////////////////////////
   //// CREATE VECTOR STRUCTURE ////
   /////////////////////////////////

   int nentries = 0;
   for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
      NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);

      // make sure the boundary is between two cellulars
      NeuronPrediction *prediction_one = boundary->PredictionOne();
      NeuronPrediction *prediction_two = boundary->PredictionTwo();
      if (!predict_all && !(prediction_one->IsOnBoundary() ^ prediction_two->IsOnBoundary())) continue;

      // create a new feature vector for this boundary
      features.push_back(std::vector<float>());

      // save the label
      //labels.push_back(prediction_one->MajorityHumanLabel() == prediction_two->MajorityHumanLabel());

      // save the identification of the boundary
      identifications.push_back(ipb);

      // increment the number of entries in the data set
      nentries++;
   }


   ////////////////////////////////////
   //// POPULATE BOUNDARY FEATURES ////
   ////////////////////////////////////

   // create the list of features
   names.push_back("boundary_mean");
   names.push_back("boundary_min");
   names.push_back("boundary_max");
   names.push_back("boundary_median");
   names.push_back("boundary_stddev");
   names.push_back("boundary_rank");
   names.push_back("boundary_scaled_ranking");
   names.push_back("nvoxels");
   names.push_back("nvoxels_proportion");

   nentries = 0;
   if (print_verbose) printf("Populating database with boundary features...\n  ");
   RNTime boundary_time;
   boundary_time.Read();
   for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
      if (print_verbose) { RNProgressBar(ipb, nd->NPredictionBoundaries()); }
      NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);

      NeuronPrediction *prediction_one = boundary->PredictionOne();
      NeuronPrediction *prediction_two = boundary->PredictionTwo();
      if (!predict_all && !(prediction_one->IsOnBoundary() ^ prediction_two->IsOnBoundary())) continue;

      //  get list of features
      features[nentries].push_back(boundary->Mean());
      features[nentries].push_back(boundary->Minimum());
      features[nentries].push_back(boundary->Maximum());
      features[nentries].push_back(boundary->Median());
      features[nentries].push_back(boundary->StdDev());
      features[nentries].push_back(boundary->BoundaryRank());
      features[nentries].push_back(boundary->BoundaryScaledRanking());
      features[nentries].push_back(boundary->CombinedVoxels());
      features[nentries].push_back(boundary->CombinedVoxelsProportion());

      // increment the number of entries in the data set
      nentries++;
   }
   if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());



   /////////////////////////////////
   //// POPULATE SHAPE FEATURES ////
   /////////////////////////////////

   names.push_back("random_walk");
   names.push_back("drunkards_walk");

   // output filenames
   nentries = 0;
   if (print_verbose) printf("Populating database with boundary random walk features...\n  ");
   boundary_time.Read();
   for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
      if (print_verbose) { RNProgressBar(ipb, nd->NPredictionBoundaries()); }
      NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);

      NeuronPrediction *prediction_one = boundary->PredictionOne();
      NeuronPrediction *prediction_two = boundary->PredictionTwo();
      if (!predict_all && !(prediction_one->IsOnBoundary() ^ prediction_two->IsOnBoundary())) continue;

      //  get list of features
      features[nentries].push_back(boundary->RandomWalk(FALSE));
      features[nentries].push_back(boundary->RandomWalk(TRUE));


      // increment the number of entries in the data set
      nentries++;
   }
   if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());

   // add all of the dijkstra features
   for (int id = 0; id < NSTATISTICS; ++id) {
      if (print_verbose) { printf("Creating features for %s...", statistic_names[id]); fflush(stdout); }
      RNTime feature_time;
      feature_time.Read();
      // get the feature filename
      char features_filename[4096];
      sprintf(features_filename, "%s/%s_dijkstra_%s.feature", postprocess_directory, root_filename, statistic_names[id]);

      // open file
      FILE *features_fp = fopen(features_filename, "rb");
      if (!features_fp) { fprintf(stderr, "Failed to read %s\n", features_filename); return 0; }

      // make sure ncellulars is correct
      int features_nboundaries;
      fread(&features_nboundaries, sizeof(int), 1, features_fp);
      rn_assertion(features_nboundaries == nd->NPredictionBoundaries());

      // read the dijkstra features
      float *dijkstra_features = new float[nd->NPredictionBoundaries()];
      fread(dijkstra_features, sizeof(float), nd->NPredictionBoundaries(), features_fp);

      // close file
      fclose(features_fp);

      // add feature to list
      names.push_back(statistic_names[id]);

      nentries = 0;
      for (int ipb = 0; ipb < nd->NPredictionBoundaries(); ++ipb) {
         NeuronPredictionBoundary *boundary = nd->PredictionBoundary(ipb);

	 NeuronPrediction *prediction_one = boundary->PredictionOne();
	 NeuronPrediction *prediction_two = boundary->PredictionTwo();
	 if (!predict_all && !(prediction_one->IsOnBoundary() ^ prediction_two->IsOnBoundary())) continue;

         float dijkstra_feature = dijkstra_features[boundary->DataIndex()];
         features[nentries].push_back(dijkstra_feature);

         // increment the number of entries in the data set
         nentries++;
      }

      // free memory
      delete[] dijkstra_features;

      if (print_verbose) printf("done in %0.2f seconds.\n", feature_time.Elapsed());
   }

   // to ease memory issues, create dataset like this
   if (print_verbose) { printf("Creating dataset..."); fflush(stdout); }
   RNTime dataset_time;
   dataset_time.Read();
   RNDataset *dataset = new RNDataset(features, labels, names, identifications, 2);
   if (print_verbose) printf("done in %0.2f seconds\n", dataset_time.Elapsed());

   return dataset;
}



static int CreateData(NeuronData *nd)
{
   // create dataset
   RNDataset *dataset = CreatePostDataset(nd);
   if (!dataset) return 0;

   // save the training data set
   char root_filename[4096];
   sprintf(root_filename, "%s", nd->Filename());
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   char dataset_filename[4096];
   sprintf(dataset_filename, "%s/%s_%s.mldb", dataset_directory, root_filename, extension);

   // write the dataset file to disk
   dataset->WriteFile(dataset_filename);

   // free up memory
   delete dataset;

   //  return success
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
	 else if (!strcmp(*argv, "-all")) { predict_all = 1; extension = "post_all"; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) { input_filename = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there are training files
   if (!input_filename) { fprintf(stderr, "Must supply neuron file\n"); return 0; }

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

   // read in data
   NeuronData *nd = ReadData(input_filename);
   if (!nd) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in prediction file
   char meta_root_filename[4096];
   sprintf(meta_root_filename, "output/ordermerge/%s", root_filename);

   R3Grid *prediction_grid = RNReadNeuronMetaRawFile(meta_root_filename);
   if (!prediction_grid) exit(-1);

   // create boundary predictions
   nd->CreatePredictions(prediction_grid);

   if (print_verbose) printf("Creating datasets for %s\n", input_filename);
   if (!CreateData(nd)) exit(-1);

   // free memory
   delete nd;
   delete prediction_grid;

   // return success
   return 0;
}
