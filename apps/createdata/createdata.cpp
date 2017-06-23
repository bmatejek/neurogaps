// Source file for the simple file conversion algorithm



// include files

#include "RNML/RNML.h"
#include "Neuron/Neuron.h"
#include <vector>
#include <limits>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int smoothing_terms = 0;
static int data_terms = 0;
static int boundary_terms = 0;
static int truncated = 0;
// feature options
static int boundary_features = 0;
static int shape_features = 0;
static int global_features = 0;
static int only_ordermerge = 0;
// string arguments
static const char *extension = NULL;
static const char *input_filename = NULL;



// directory structure

static const char *distances_directory = "distances";
static const char *dataset_directory = "algs_data/dataset";
static const char *features_directory = "algs_data/features";
static const char *boundaries_directory = "boundaries";
static const char *ordermerge_output_directory = "output/ordermerge";
static const char *boundary_data_directory = "algs_data/boundarydata";



// useful constants

static const int NBOUNDARY_FEATURES = 7;
static const char *boundary_features_names[NBOUNDARY_FEATURES] = {
   "skew", "kurtosis", "fewer_40", "fewer_60", "fewer_80",
   "degree_differences", "mutual_neighbors"
};



static const int NSHAPE_FEATURES = 14;
static const char *shape_feature_names[NSHAPE_FEATURES] = {
   "first_derivative_mean", "first_derivative_median", "first_derivative_maximum",
   "first_derivative_minimum", "first_derivative_stddev", "first_derivative_skew", "first_derivative_kurtosis",
   "second_derivative_mean", "second_derivative_median", "second_derivative_maximum",
   "second_derivative_minimum", "second_derivative_stddev", "second_derivative_skew", "second_derivative_kurtosis"
};



static const int NDIJKSTRA_FEATURES = 12;
static const char *dijkstra_features_names[NDIJKSTRA_FEATURES] = {
   "nsteps", "kurtosis", "maximum", "mean",
   "median", "minimum", "skew", "stddev", "first_quintile",
   "second_quintile", "third_quintile", "fourth_quintile"
};



static const int NDISTANCE_FEATURES = 1;
static const char *distance_features_names[NDISTANCE_FEATURES] = {
   "dijkstra"
};


static const int NGLOBAL_ALGORITHMS = 4;
static const char *global_algorithm_names[NGLOBAL_ALGORITHMS] = {
   "facility_location", "kmedian", "ordermerge", "boundarymerge"//, "graphcut"
};



static const int NFEATURES_PER_GLOBAL = 2;
static const char *global_feature_names[NFEATURES_PER_GLOBAL] = {
   "first_split", "proportion"
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

static RNDataset *CreateSmoothingDataset(NeuronData *nd)
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

   if (boundary_features) {
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
         features[nentries].push_back(boundary->CombinedVoxels());
         features[nentries].push_back(boundary->CombinedVoxelsProportion());

         // increment the number of entries in the data set
         nentries++;
      }
      if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());

      for (int is = 0; is < NBOUNDARY_FEATURES; ++is) {
         char filename[4096];
         sprintf(filename, "%s/%s_%s.feature", boundary_data_directory, root_filename, boundary_features_names[is]);

         // open file
         FILE *fp = fopen(filename, "rb");
         if (!fp) { fprintf(stderr, "Failed to read %s\n", filename); return NULL; }

         // read the number of boundaries
         int nboundaries;
         fread(&nboundaries, sizeof(int), 1, fp);
         rn_assertion(nboundaries == nd->NBoundaries());

         RNScalar *values = new RNScalar[nd->NBoundaries()];
         if (fread(values, sizeof(RNScalar), nd->NBoundaries(), fp) != (unsigned int)nd->NBoundaries()) {
            fprintf(stderr, "Failed to read %s\n", filename);
            return 0;
         }

         names.push_back(boundary_features_names[is]);

         nentries = 0;
         RNTime boundary_time;
         boundary_time.Read();
         for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
            NeuronBoundary *boundary = nd->Boundary(ib);

            // make sure the boundary is between two cellulars
            NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
            if (supervoxel_one->IsExtracellular()) continue;
            NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
            if (supervoxel_two->IsExtracellular()) continue;

            //  get list of features
            features[nentries].push_back(values[ib]);

            // increment the number of entries in the data set
            nentries++;
         }

         // free memory
         delete[] values;
      }
   }


   /////////////////////////////////
   //// POPULATE SHAPE FEATURES ////
   /////////////////////////////////

   if (shape_features) {
      for (int is = 0; is < NSHAPE_FEATURES; ++is) {
         char filename[4096];
         sprintf(filename, "%s/%s_%s.feature", boundary_data_directory, root_filename, shape_feature_names[is]);

         // open file
         FILE *fp = fopen(filename, "rb");
         if (!fp) { fprintf(stderr, "Failed to read %s\n", filename); return NULL; }

         // read the number of boundaries
         int nboundaries;
         fread(&nboundaries, sizeof(int), 1, fp);
         rn_assertion(nboundaries == nd->NBoundaries());

         RNScalar *values = new RNScalar[nd->NBoundaries()];
         if (fread(values, sizeof(RNScalar), nd->NBoundaries(), fp) != (unsigned int)nd->NBoundaries()) {
            fprintf(stderr, "Failed to read %s\n", filename);
            return 0;
         }

         names.push_back(shape_feature_names[is]);

         nentries = 0;
         RNTime boundary_time;
         boundary_time.Read();
         for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
            NeuronBoundary *boundary = nd->Boundary(ib);

            // make sure the boundary is between two cellulars
            NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
            if (supervoxel_one->IsExtracellular()) continue;
            NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
            if (supervoxel_two->IsExtracellular()) continue;

            //  get list of features
            features[nentries].push_back(values[ib]);

            // increment the number of entries in the data set
            nentries++;
         }

         // free memory
         delete[] values;
      }

      names.push_back("random_walk");
      names.push_back("drunkards_walk");
      names.push_back("monte_carlo_ratio");
      names.push_back("monte_carlo_raw");

      // output filenames
      char monte_carlo_ratio_filename[4096];
      sprintf(monte_carlo_ratio_filename, "%s/%s_monte_carlo_ratio.feature", features_directory, root_filename);
      char monte_carlo_filename[4096];
      sprintf(monte_carlo_filename, "%s/%s_monte_carlo.feature", features_directory, root_filename);

      // open file
      FILE *monte_carlo_ratio_fp = fopen(monte_carlo_ratio_filename, "rb");
      if (!monte_carlo_ratio_fp) { fprintf(stderr, "Failed to read %s\n", monte_carlo_ratio_filename); return NULL; }

      FILE *monte_carlo_fp = fopen(monte_carlo_filename, "rb");
      if (!monte_carlo_fp) { fprintf(stderr, "Failed to read %s\n", monte_carlo_filename); return NULL; }

      int ratio_nboundaries;
      int raw_nboundaries;
      fread(&ratio_nboundaries, sizeof(int), 1, monte_carlo_ratio_fp);
      rn_assertion(ratio_nboundaries == nd->NBoundaries());
      fread(&raw_nboundaries, sizeof(int), 1, monte_carlo_fp);
      rn_assertion(raw_nboundaries == nd->NBoundaries());

      float *monte_carlo_ratios = new float[ratio_nboundaries];
      if (fread(monte_carlo_ratios, sizeof(float), ratio_nboundaries, monte_carlo_ratio_fp) != (unsigned int)ratio_nboundaries) {
         fprintf(stderr, "Failed to read %s\n", monte_carlo_ratio_filename);
         return NULL;
      }
      unsigned int *monte_carlo_raw = new unsigned int[raw_nboundaries];
      if (fread(monte_carlo_raw, sizeof(unsigned int), raw_nboundaries, monte_carlo_fp) != (unsigned int)raw_nboundaries) {
         fprintf(stderr, "Failed to read %s\n", monte_carlo_filename);
         return NULL;
      }

      // close files
      fclose(monte_carlo_ratio_fp);
      fclose(monte_carlo_fp);

      nentries = 0;
      if (print_verbose) printf("Populating database with boundary random walk features...\n  ");
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
         features[nentries].push_back(boundary->RandomWalk(FALSE));
         features[nentries].push_back(boundary->RandomWalk(TRUE));
         features[nentries].push_back(monte_carlo_ratios[ib]);
         features[nentries].push_back((float)monte_carlo_raw[ib]);

         // increment the number of entries in the data set
         nentries++;
      }
      if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());

      // free memory
      delete[] monte_carlo_ratios;
      delete[] monte_carlo_raw;

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
            if (!strcmp(dijkstra_features_names[id], "nsteps")) {
               int *dijkstra_nsteps = new int[nd->NCellulars()];

               if (fread(dijkstra_nsteps, sizeof(int), nd->NCellulars(), features_fp) != (unsigned int)nd->NCellulars()) {
                  fprintf(stderr, "Failed to read %s\n", features_filename);
                  return NULL;
               }


               for (int j = 0; j < nd->NCellulars(); ++j) {
                  dijkstra_features[ic][j] = (float)dijkstra_nsteps[j];
               }

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
               return NULL;
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
   }


   //////////////////////////////////
   //// POPULATE GLOBAL FEATURES ////
   //////////////////////////////////

   if (global_features) {
      for (int ig = 0; ig < NGLOBAL_ALGORITHMS; ++ig) {
         for (int ifg = 0; ifg < NFEATURES_PER_GLOBAL; ++ifg) {
            char feature_name[4096];
            sprintf(feature_name, "%s_%s", global_algorithm_names[ig], global_feature_names[ifg]);

            char feature_filename[4096];
            sprintf(feature_filename, "%s/%s_%s.feature", features_directory, root_filename, feature_name);

            // open file
            FILE *global_fp = fopen(feature_filename, "rb");
            if (!global_fp) { fprintf(stderr, "Failed to read %s\n", feature_filename); return NULL; }

            int nglobal_boundaries;
            fread(&nglobal_boundaries, sizeof(int), 1, global_fp);
            rn_assertion(nglobal_boundaries == nd->NBoundaries());

            float *global_features = new float[nd->NBoundaries()];

            if (!strcmp(global_feature_names[ifg], "proportion")) {
               if (fread(global_features, sizeof(float), nd->NBoundaries(), global_fp) != (unsigned int)nd->NBoundaries()) {
                  fprintf(stderr, "Failed to read %s\n", feature_filename);
                  return NULL;
               }
            }
            else {
               int *tmp_features = new int[nd->NBoundaries()];
               if (fread(tmp_features, sizeof(int), nd->NBoundaries(), global_fp) != (unsigned int)nd->NBoundaries()) {
                  fprintf(stderr, "Failed to read %s\n", feature_filename);
                  return NULL;
               }
               for (int ib = 0; ib < nd->NBoundaries(); ++ib)
                  global_features[ib] = (float)tmp_features[ib];

               delete[] tmp_features;
            }

            names.push_back(feature_name);

            nentries = 0;
            if (print_verbose) { printf("Creating features for %s...", feature_name); fflush(stdout); }
            RNTime boundary_time;
            boundary_time.Read();
            for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
               NeuronBoundary *boundary = nd->Boundary(ib);

               // make sure the boundary is between two cellulars
               NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
               if (supervoxel_one->IsExtracellular()) continue;
               NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
               if (supervoxel_two->IsExtracellular()) continue;

               //  get list of features
               features[nentries].push_back(global_features[ib]);

               // increment the number of entries in the data set
               nentries++;
            }
            if (print_verbose) printf("done in %0.2f seconds.\n", boundary_time.Elapsed());


            // free memory
            delete[] global_features;

            // close file
            fclose(global_fp);
         }
      }

      // read ordermerge results
      char ordermerge_filename[4096];
      sprintf(ordermerge_filename, "%s/%s_ordermerge.boundary", boundaries_directory, root_filename);
      char boundarymerge_filename[4096];
      sprintf(boundarymerge_filename, "%s/%s_boundarymerge.boundary", boundaries_directory, root_filename);
      //char graphcut_filename[4096];
      //sprintf(graphcut_filename, "%s/%s_graphcut.boundary", boundaries_directory, root_filename);

      // open file
      FILE *ordermerge_fp = fopen(ordermerge_filename, "rb");
      if (!ordermerge_fp) { fprintf(stderr, "Failed to read %s\n", ordermerge_filename); return NULL; }
      FILE *boundarymerge_fp = fopen(boundarymerge_filename, "rb");
      if (!boundarymerge_fp) { fprintf(stderr, "Failed to read %s\n", boundarymerge_filename); return NULL; }
      //FILE *graphcut_fp = fopen(graphcut_filename, "rb");
      //if (!graphcut_fp) { fprintf(stderr, "Failed to read %s\n", graphcut_filename); return NULL; }

      // read number of boundaries
      int nordermerge_boundaries;
      fread(&nordermerge_boundaries, sizeof(int), 1, ordermerge_fp);
      rn_assertion(nordermerge_boundaries == nd->NBoundaries());
      int nboundarymerge_boundaries;
      fread(&nboundarymerge_boundaries, sizeof(int), 1, boundarymerge_fp);
      rn_assertion(nboundarymerge_boundaries == nd->NBoundaries());
      //int ngraphcut_boundaries;
      //fread(&ngraphcut_boundaries, sizeof(int), 1, graphcut_fp);
      //rn_assertion(ngraphcut_boundaries == nd->NBoundaries());

      RNBoolean *ordermerge_boundaries = new RNBoolean[nd->NBoundaries()];
      if (fread(ordermerge_boundaries, sizeof(RNBoolean), nd->NBoundaries(), ordermerge_fp) != (unsigned int)nd->NBoundaries()) {
         fprintf(stderr, "Failed to read %s\n", ordermerge_filename);
         return NULL;
      }
      
      RNBoolean *boundarymerge_boundaries = new RNBoolean[nd->NBoundaries()];
      if (fread(boundarymerge_boundaries, sizeof(RNBoolean), nd->NBoundaries(), boundarymerge_fp) != (unsigned int)nd->NBoundaries()) {
         fprintf(stderr, "Failed to read %s\n", boundarymerge_filename);
         return NULL;
      }

      /*RNBoolean *graphcut_boundaries = new RNBoolean[nd->NBoundaries()];
      if (fread(graphcut_boundaries, sizeof(RNBoolean), nd->NBoundaries(), graphcut_fp) != (unsigned int)nd->NBoundaries()) {
      fprintf(stderr, "Failed to read %s\n", graphcut_filename);
      return NULL;
      }*/

      names.push_back("ordermerge_boundary");
      names.push_back("boundarymerge_boundary");
      //names.push_back("graphcut_boundary");

      nentries = 0;
      if (print_verbose) printf("Populating database with global features...\n  ");
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
         features[nentries].push_back(ordermerge_boundaries[ib]);
         features[nentries].push_back(boundarymerge_boundaries[ib]);
         //features[nentries].push_back(graphcut_boundaries[ib]);

         // increment the number of entries in the data set
         nentries++;
      }
      if (print_verbose) printf("\ndone in %0.2f seconds.\n", boundary_time.Elapsed());

      // free memory
      delete[] ordermerge_boundaries;
      delete[] boundarymerge_boundaries;
      //delete[] graphcut_boundaries;

      // close file
      fclose(ordermerge_fp);
      //fclose(boundarymerge_fp);
      //fclose(graphcut_fp);
   }

   // to ease memory issues, create dataset like this
   if (print_verbose) { printf("Creating dataset..."); fflush(stdout); }
   RNTime dataset_time;
   dataset_time.Read();
   RNDataset *smoothing_dataset = new RNDataset(features, labels, names, identifications, 2);
   if (print_verbose) printf("done in %0.2f seconds\n", dataset_time.Elapsed());

   return smoothing_dataset;
}



static RNDataset *CreateDataDataset(NeuronData *nd)
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

   // count the number of positive and negative examples
   unsigned int positive_examples = 0;
   unsigned int negative_examples = 0;
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      NeuronCellular *cellular_one = nd->Cellular(ic1);
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);

         // onlyy consider pairs where cellular two is on the boundary
         if (!cellular_two->IsOnBoundary()) continue;

         // get the truth label
         unsigned int label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());

         if (label) positive_examples++;
         else negative_examples++;
      }
   }

   // go through all pairs of cellulars
   if (print_verbose) { printf("Creating structure for data terms...\n  "); fflush(stdout); }
   RNTime structure_time;
   structure_time.Read();
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      if (print_verbose) RNProgressBar(ic1, nd->NCellulars());
      NeuronCellular *cellular_one = nd->Cellular(ic1);

      // boundary data set only has pair of boundary supervoxels
      if (boundary_terms) {
         if (!cellular_one->IsOnBoundary()) continue;
      }

      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);

         // only consider pairs where cellular two is on the boundary
         if (!cellular_two->IsOnBoundary()) continue;

         // get the truth label
         unsigned int label = (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());

         // throw out some examples if we want truncated results
         if (truncated) {
            RNScalar proportion_positive = 30000 / (RNScalar)positive_examples;
            RNScalar proportion_negative = 30000 / (RNScalar)negative_examples;

            RNScalar random_value = RNRandomScalar();
            if (label && (random_value > proportion_positive)) continue;
            if (!label && (random_value > proportion_negative)) continue;
         }


         // save the label
         labels.push_back(label);

         // create a new feature vector
         features.push_back(std::vector<float>());

         // save the identification of the boundary
         int id = ic1 * nd->NCellulars() + ic2;
         identifications.push_back(id);
      }
   }
   if (print_verbose) printf("\ndone in %0.2f seconds.\n", structure_time.Elapsed());


   ////////////////////////////////////
   //// POPULATE DIJKSTRA FEATURES ////
   ////////////////////////////////////

   if (!only_ordermerge) {
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
            if (!strcmp(dijkstra_features_names[id], "nsteps")) {
               int *dijkstra_nsteps = new int[nd->NCellulars()];
               if (fread(dijkstra_nsteps, sizeof(int), nd->NCellulars(), features_fp) != (unsigned int)nd->NCellulars()) {
                  fprintf(stderr, "Failed to read %s\n", features_filename);
                  return NULL;
               }
               for (int j = 0; j < nd->NCellulars(); ++j) {
                  dijkstra_features[ic][j] = (float)dijkstra_nsteps[j];
               }
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

         // add features to vectors
         names.push_back(dijkstra_features_names[id]);

         // go through all features
         for (unsigned int ie = 0; ie < features.size(); ++ie) {
            int id = identifications[ie];
            int ic1 = id / nd->NCellulars();
            int ic2 = id % nd->NCellulars();

            float dijkstra_feature = dijkstra_features[ic1][ic2];
            features[ie].push_back(dijkstra_feature);
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
               return NULL;
            }
         }

         // close file
         fclose(distances_fp);

         // add features to vectors
         names.push_back(distance_features_names[id]);

         // go through all features
         for (unsigned int ie = 0; ie < features.size(); ++ie) {
            int id = identifications[ie];
            int ic1 = id / nd->NCellulars();
            int ic2 = id % nd->NCellulars();

            float distance_feature = distances[ic1][ic2];

            features[ie].push_back(distance_feature);
         }

         // free memory
         for (int ic = 0; ic < nd->NCellulars(); ++ic) {
            delete[] distances[ic];
         }
         delete[] distances;

         if (print_verbose) printf("done in %0.2f seconds.\n", distance_time.Elapsed());
      }
   }

   // add ordermerge file
   char ordermerge_filename[4096];
   sprintf(ordermerge_filename, "%s/%s", ordermerge_output_directory, root_filename);
   //char boundarymerge_filename[4096];
   //sprintf(boundarymerge_filename, "%s/%s", ordermerge_output_directory, root_filename);

   R3Grid *ordermerge_prediction = RNReadNeuronMetaRawFile(ordermerge_filename);
   if (!ordermerge_prediction) { fprintf(stderr, "Failed to read %s\n", ordermerge_filename); return NULL; }
   //R3Grid *boundarymerge_prediction = RNReadNeuronMetaRawFile(boundarymerge_filename);
   //if (!boundarymerge_prediction) { fprintf(stderr, "Failed to read %s\n", boundarymerge_filename); return NULL; }

   names.push_back("ordermerge_output");
   //names.push_back("boundarymerge_output");

   // read voxels
   nd->ReadVoxels();

   // go through all features
   for (unsigned int ie = 0; ie < features.size(); ++ie) {
      int id = identifications[ie];
      int ic1 = id / nd->NCellulars();
      int ic2 = id % nd->NCellulars();

      NeuronCellular *cellular_one = nd->Cellular(ic1);
      NeuronCellular *cellular_two = nd->Cellular(ic2);

      NeuronVoxel *voxel_one = cellular_one->CenterVoxel();
      NeuronVoxel *voxel_two = cellular_two->CenterVoxel();

      int ordermerge_prediction_one = ordermerge_prediction->GridValue(voxel_one->DataIndex());
      int ordermerge_prediction_two = ordermerge_prediction->GridValue(voxel_two->DataIndex());
      //int boundarymerge_prediction_one = boundarymerge_prediction->GridValue(voxel_one->DataIndex());
      //int boundarymerge_prediction_two = boundarymerge_prediction->GridValue(voxel_two->DataIndex());

      float ordermerge_feature = (float)(ordermerge_prediction_one == ordermerge_prediction_two);
      //float boundarymerge_feature = (float)(boundarymerge_prediction_one == boundarymerge_prediction_two);

      features[ie].push_back(ordermerge_feature);
      //features[ie].push_back(boundarymerge_feature);
   }

   // free memory
   delete ordermerge_prediction;
   //delete boundarymerge_prediction;

   // release voxels
   nd->ReleaseVoxels();

   // to ease memory issues, create dataset like this
   if (print_verbose) { printf("Creating dataset..."); fflush(stdout); }
   RNTime dataset_time;
   dataset_time.Read();
   RNDataset *data_dataset = new RNDataset(features, labels, names, identifications, 2);
   if (print_verbose) printf("done in %0.2f seconds\n", dataset_time.Elapsed());

   return data_dataset;
}



static int CreateData(NeuronData *nd)
{
   // create dataset
   RNDataset *dataset = NULL;
   if (smoothing_terms) dataset = CreateSmoothingDataset(nd);
   else dataset = CreateDataDataset(nd);
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
         else if (!strcmp(*argv, "-smoothing_term")) { smoothing_terms = 1; }
         else if (!strcmp(*argv, "-data_term")) { data_terms = 1; }
         else if (!strcmp(*argv, "-boundary_term")) { boundary_terms = 1; }
         // feature options
         else if (!strcmp(*argv, "-boundary_features")) { boundary_features = 1; }
         else if (!strcmp(*argv, "-shape_features")) { shape_features = 1; }
         else if (!strcmp(*argv, "-global_features")) { global_features = 1; }
         else if (!strcmp(*argv, "-only_ordermerge")) { only_ordermerge = 1; }
         // truncated option for data terms
         else if (!strcmp(*argv, "-truncated")) { truncated = 1; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) { input_filename = *argv; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there are consistent labeling options
   if (data_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      // check truncated flag
      if (truncated) extension = "data_truncated";
      else extension = "data";
   }
   else if (boundary_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      if (truncated) { fprintf(stderr, "Can only truncate data terms\n"); return 0; }

      data_terms = 1;
      extension = "data_boundary";
   }
   else if (smoothing_terms) {
      if (truncated) { fprintf(stderr, "Can only truncate data terms\n"); return 0; }

      if (global_features) extension = "smoothing_global";
      else extension = "smoothing";

      // turn on all feature options besides global
      boundary_features = 1;
      shape_features = 1;
   }
   else if (boundary_features || shape_features || global_features) {
      if (truncated) { fprintf(stderr, "Can only truncate data terms\n"); return 0; }

      // make sure only one feature is selected
      if (boundary_features + shape_features + global_features > 1) {
         fprintf(stderr, "Can only select one among boundary, shape, and global features.\n");
         return 0;
      }

      // get the correct extension
      if (boundary_features) extension = "boundary";
      if (shape_features) extension = "shape";
      if (global_features) extension = "global";

      // turn on smoothing terms to go to the correct functions
      smoothing_terms = 1;
   }
   else {
      fprintf(stderr, "Must select either data or smoothing terms.\n");
      return 0;
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

   if (print_verbose) printf("Creating datasets for %s\n", input_filename);
   if (!CreateData(nd)) exit(-1);

   // free memory
   delete nd;

   // return success
   return 0;
}
