// Source file for the simple file conversion algorithm



// include files

#include "RNML/RNML.h"
#include "Neuron/Neuron.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static std::vector<const char *> neuron_files = std::vector<const char *>();



// directory structure

static const char *distances_directory = "distances";
static const char *dataset_directory = "algs_data/dataset";
static const char *features_directory = "algs_data/features";
static const char *boundaries_directory = "boundaries";



// useful constants

static const int NDIJKSTRA_FEATURES = 12;
static const char *dijkstra_features_names[NDIJKSTRA_FEATURES] = {
   "nsteps", "kurtosis", "maximum", "mean",
   "median", "minimum", "skew", "stddev", "first_quintile",
   "second_quintile", "third_quintile", "fourth_quintile"
};



static const int NDISTANCE_FEATURES = 3;
static const char *distance_features_names[NDISTANCE_FEATURES] = {
   "average_commute", "dijkstra", "first_pass"
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
   if (!nd->ReadFile(filename)) {
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

static RNDataset *CreateDataset(NeuronData *nd)
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
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // make sure the boundary is between two cellulars
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      if (supervoxel_one->IsExtracellular()) continue;
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_two->IsExtracellular()) continue;

      // create a new feature vector for this boundary
      features.push_back(std::vector<float>());

      // save the label
      labels.push_back(supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel());

      // save the identification of the boundary
      identifications.push_back(ib);

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
   names.push_back("random_walk");
   names.push_back("drunkards_walk");
   names.push_back("nvoxels");
   names.push_back("nvoxels_proportion");

   nentries = 0;
   if (print_verbose) printf("Populating database with boundary features...\n  ");
   RNTime boundary_time;
   boundary_time.Read();
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      if (print_verbose) { RNProgressBar(ib, nd->NBoundaries()); }
      NeuronBoundary *boundary = nd->Boundary(ib);

      // make sure the boundary is between two cellulars
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      if (supervoxel_one->IsExtracellular()) continue;
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_two->IsExtracellular()) continue;

      //  get list of features
      features[nentries].push_back(boundary->Mean());
      features[nentries].push_back(boundary->Minimum());
      features[nentries].push_back(boundary->Maximum());
      features[nentries].push_back(boundary->Median());
      features[nentries].push_back(boundary->StdDev());
      features[nentries].push_back(boundary->BoundaryRank());
      features[nentries].push_back(boundary->BoundaryScaledRanking());
      features[nentries].push_back(boundary->RandomWalk(FALSE));
      features[nentries].push_back(boundary->RandomWalk(TRUE));
      features[nentries].push_back(boundary->CombinedVoxels());
      features[nentries].push_back(boundary->CombinedVoxelsProportion());

      // increment the number of entries in the data set
      nentries++;
   }
   if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());


   ////////////////////////////////////
   //// POPULATE DIJKSTRA FEATURES ////
   ////////////////////////////////////

   // add all of the dijkstra features
   for (int id = 0; id < NDIJKSTRA_FEATURES; ++id) {
      if (print_verbose) { printf("Creating features for %s...", dijkstra_features_names[id]); fflush(stdout); }
      RNTime feature_time;
      feature_time.Read();
      // get the feature filename
      char features_filename[4096];
      sprintf(features_filename, "%s/%s_dijkstra_%s.feature", features_directory, root_filename, dijkstra_features_names[id]);

      // open file
      FILE *features_fp = fopen(features_filename, "rb");
      if (!features_fp) { fprintf(stderr, "Failed to read %s\n", features_filename); return 0; }

      // make sure ncellulars is correct
      int features_ncellulars;
      fread(&features_ncellulars, sizeof(int), 1, features_fp);
      rn_assertion(features_ncellulars == nd->NCellulars());

      // read the dijkstra features
      float **dijkstra_features = new float *[nd->NCellulars()];
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         dijkstra_features[ic] = new float[nd->NCellulars()];

         // steps are recoreded as integers
         if (!strcmp(dijkstra_features_names[id], "dijkstra_nsteps")) {
            int *dijkstra_nsteps = new int[nd->NCellulars()];
	    if (fread(dijkstra_nsteps, sizeof(int), nd->NCellulars(), features_fp) != (unsigned int)nd->NCellulars()) {
               fprintf(stderr, "Failed to read %s\n", features_filename);
               return NULL;
            }
            
            // translate integer features to floats      
            for (int j = 0; j < nd->NCellulars(); ++j) {
               dijkstra_features[ic][j] = (float)dijkstra_nsteps[j];
            }

            // free tmp memory
            delete[] dijkstra_nsteps;
         }
	 else {
	    if (fread(dijkstra_features[ic], sizeof(float), nd->NCellulars(), features_fp) != (unsigned int)nd->NCellulars()) {
	       fprintf(stderr, "Failed to read %s\n", features_filename);
	       return NULL;
	    }
	 }
      }

      // close file
      fclose(features_fp);

      // add feature to list
      names.push_back(dijkstra_features_names[id]);

      nentries = 0;
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);

         // make sure the boundary is between two cellulars
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         if (supervoxel_one->IsExtracellular()) continue;
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
         if (supervoxel_two->IsExtracellular()) continue;

         // get data index
         int cellular_one_data_index = supervoxel_one->DataIndex();
         int cellular_two_data_index = supervoxel_two->DataIndex();

         float dijkstra_feature = dijkstra_features[cellular_one_data_index][cellular_two_data_index];
         features[nentries].push_back(dijkstra_feature);

         // increment the number of entries in the data set
         nentries++;
      }

      // free memory
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         delete[] dijkstra_features[ic];
      }
      delete[] dijkstra_features;

      if (print_verbose) printf("done in %0.2f seconds.\n", feature_time.Elapsed());
   }


   ////////////////////////////////////
   //// POPULATE DISTANCE FEATURES ////
   ////////////////////////////////////

   // add all of the distance features
   for (int id = 0; id < NDISTANCE_FEATURES; ++id) {
      if (print_verbose) { printf("Creating features for %s...", distance_features_names[id]); fflush(stdout); }
      RNTime distance_time;
      distance_time.Read();

      // get the distances filename
      char distances_filename[4096];
      sprintf(distances_filename, "%s/%s_%s.distance", distances_directory, root_filename, distance_features_names[id]);

      // open file
      FILE *distances_fp = fopen(distances_filename, "rb");
      if (!distances_fp) { fprintf(stderr, "Failed to read %s\n", distances_filename); return 0; }

      // make sure ncellulars is correct
      int distance_ncellulars;
      fread(&distance_ncellulars, sizeof(int), 1, distances_fp);
      rn_assertion(distance_ncellulars == nd->NCellulars());

      // read the dijkstra features
      float **distances = new float *[nd->NCellulars()];
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         distances[ic] = new float[nd->NCellulars()];
         if (fread(distances[ic], sizeof(float), nd->NCellulars(), distances_fp) != (unsigned int)nd->NCellulars()) {
            fprintf(stderr, "Failed to read %s\n", distances_filename);
            return 0;
         }
      }

      // close file
      fclose(distances_fp);

      // add feature to vectors
      names.push_back(distance_features_names[id]);

      nentries = 0;
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);

         // make sure the boundary is between two cellulars
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         if (supervoxel_one->IsExtracellular()) continue;
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
         if (supervoxel_two->IsExtracellular()) continue;

         // get data index
         int cellular_one_data_index = supervoxel_one->DataIndex();
         int cellular_two_data_index = supervoxel_two->DataIndex();

         float distance_feature = distances[cellular_one_data_index][cellular_two_data_index];
         features[nentries].push_back(distance_feature);

         // increment the number of entries in the data set
         nentries++;
      }

      // free memory
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         delete[] distances[ic];
      }
      delete[] distances;

      if (print_verbose) printf("done in %0.2f seconds.\n", distance_time.Elapsed());
   }


   //////////////////////////////////
   //// POPULATE GLOBAL FEATURES ////
   //////////////////////////////////

   // add k-median and facility location results
   /*for (int K = 2; K <= 1000; ++K) {
      if (K % 20 != 0) continue;
      char kmedian_filename[4096];
      sprintf(kmedian_filename, "%s/%s_kmedian_%04d.boundary", boundaries_directory, root_filename, K);
      char facility_filename[4096];
      sprintf(facility_filename, "%s/%s_facility_%04d.boundary", boundaries_directory, root_filename, K);

      FILE *kmedian_fp = fopen(kmedian_filename, "rb");
      if (!kmedian_fp) { fprintf(stderr, "Failed to read %s\n", kmedian_filename); return NULL; }
      FILE *facility_fp = fopen(facility_filename, "rb");
      if (!facility_fp) { fprintf(stderr, "Failed to read %s\n", facility_filename); return NULL; }

      int kmedian_nboundaries = 0;
      fread(&kmedian_nboundaries, sizeof(int), 1, kmedian_fp);
      int facility_nboundaries = 0;
      fread(&facility_nboundaries, sizeof(int), 1, facility_fp);

      // read in boundaries
      RNBoolean *kmedian_boundaries = new RNBoolean[kmedian_nboundaries];
      if (fread(kmedian_boundaries, sizeof(RNBoolean), kmedian_nboundaries, kmedian_fp) != (unsigned int) kmedian_nboundaries) {
	 fprintf(stderr, "Failed to read %s\n", kmedian_filename);
	 return NULL;
      }
      RNBoolean *facility_boundaries = new RNBoolean[facility_nboundaries];
      if (fread(facility_boundaries, sizeof(RNBoolean), facility_nboundaries, facility_fp) != (unsigned int) facility_nboundaries) {
	 fprintf(stderr, "Failed to read %s\n", facility_filename);
	 return NULL;
      }
   
      // get feature names
      char kmedian_feature_name[4096];
      sprintf(kmedian_feature_name, "kmedian_%04d", K);
      char facility_feature_name[4096];
      sprintf(facility_feature_name, "facility_%04d", K);
   
      // add features to vectors
      names.push_back(kmedian_feature_name);
      names.push_back(facility_feature_name);
      
      nentries = 0;
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
	 NeuronBoundary *boundary = nd->Boundary(ib);
	 
	 // make sure the boundary is between two cellulars
	 NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
	 if (supervoxel_one->IsExtracellular()) continue;
	 NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
	 if (supervoxel_two->IsExtracellular()) continue;

	 // insert feature
	 features[nentries].push_back(kmedian_boundaries[ib]);
	 features[nentries].push_back(facility_boundaries[ib]);

	 // increment the number of entries in the data set
	 nentries++;
      }

      // free memory
      delete[] kmedian_boundaries;
      delete[] facility_boundaries;

      // close files
      fclose(kmedian_fp);
      fclose(facility_fp);
      }*/

   // add in order merge results
   char ordermerge_filename[4096];
   sprintf(ordermerge_filename, "%s/%s_ordermerge.boundary", boundaries_directory, root_filename);

   // open file
   FILE *ordermerge_fp = fopen(ordermerge_filename, "rb");
   if (!ordermerge_fp) { fprintf(stderr, "Failed to read %s\n", ordermerge_filename); return NULL; }

   int ordermerge_nboundaries = 0;
   fread(&ordermerge_nboundaries, sizeof(int), 1, ordermerge_fp);
   rn_assertion(ordermerge_nboundaries == nd->NBoundaries());

   // read in boundaries
   RNBoolean *ordermerge_boundaries = new RNBoolean[ordermerge_nboundaries];
   if (fread(ordermerge_boundaries, sizeof(RNBoolean), ordermerge_nboundaries, ordermerge_fp) != (unsigned int)ordermerge_nboundaries) {
      fprintf(stderr, "Failed to read %s\n", ordermerge_filename);
      return NULL;
   }

   // add feature to vectors
   names.push_back("ordermerge_results");

   nentries = 0;
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // make sure the boundary is between two cellulars
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      if (supervoxel_one->IsExtracellular()) continue;
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_two->IsExtracellular()) continue;

      // insert feature
      features[nentries].push_back(ordermerge_boundaries[ib]);

      // increment the number of entries in the data set
      nentries++;
   }

   // free memory
   delete[] ordermerge_boundaries;

   // close file
   fclose(ordermerge_fp);

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
   RNDataset *dataset = CreateDataset(nd);
   if (!dataset) return 0;

   // save the training data set
   char root_filename[4096];
   sprintf(root_filename, "%s", nd->Filename());
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   char dataset_filename[4096];
   sprintf(dataset_filename, "%s/%s_global.mldb", dataset_directory, root_filename);

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
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         neuron_files.push_back(*argv);
      }
      argv++; argc--;
   }

   if (!neuron_files.size()) { fprintf(stderr, "Must supply at least on neuron file\n"); return 0; }

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

   // create data for all supplied neuron files
   for (unsigned int in = 0; in < neuron_files.size(); ++in) {
      // read in training data
      NeuronData *nd = ReadData(neuron_files[in]);
      if (!nd) exit(-1);

      if (print_verbose) printf("Creating datasets for %s\n", neuron_files[in]);
      if (!CreateData(nd)) exit(-1);

      // free memory
      delete nd;
   }

   // return success
   return 0;
}
