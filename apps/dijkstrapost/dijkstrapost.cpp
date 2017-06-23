// Source file for the simple file conversion 



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static RNScalar max_affinity = 0.82;
static RNScalar scaling[3];
static int nsteps = 20;
static int starting_prediction_index = -1;
static int merge_results = 0;



// global variables

static NeuronData *nd = NULL;
static R3Grid *prediction_grid = NULL;
static int npredictions = -1;


// directory structure

static const char *tmp_postprocess_directory = "postprocess/tmp";
static const char *postprocess_directory = "postprocess";



// useful character arrays

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



////////////////////////////////////////////////////////////////////////
// Data functions
////////////////////////////////////////////////////////////////////////

static int Ranking(std::vector<RNScalar> &attributes, int index)
{
   RNScalar value = attributes[index];

   int rank = 1;
   for (unsigned int ia = 0; ia < attributes.size(); ++ia) {
      if (attributes[ia] < value) rank++;
   }

   return rank;
}



// struct for dijkstra voxel node
struct DijkstraVoxelNode
{
   NeuronVoxel *voxel;
   DijkstraVoxelNode *prev;
   RNScalar distance;
   RNBoolean visited;
};



static int Dijkstra(int prediction_index, char root_filename[4096])
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   NeuronPrediction *prediction = nd->Prediction(prediction_index);

   // allocate temporary data
   DijkstraVoxelNode *voxel_data = new DijkstraVoxelNode[nd->NVoxels()];
   if (!voxel_data) {
      fprintf(stderr, "Unable to allocate temporary data for geodisic distances\n");
      return 0;
   }

   // initialize all data
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      voxel_data[iv].voxel = nd->Voxel(iv);
      voxel_data[iv].prev = NULL;
      voxel_data[iv].distance = FLT_MAX;
      voxel_data[iv].visited = FALSE;
   }

   // initialize priority queue with all of the boundary sources
   DijkstraVoxelNode tmp;
   RNMinBinaryHeap<DijkstraVoxelNode *> voxel_heap(&tmp, &(tmp.distance), nd->NVoxels());

   // get source index
   int source_index = prediction->CenterVoxel()->DataIndex();
   voxel_data[source_index].voxel = nd->Voxel(source_index);
   voxel_data[source_index].prev = NULL;
   voxel_data[source_index].distance = 0.0;
   voxel_data[source_index].visited = TRUE;
   voxel_heap.Insert(source_index, &(voxel_data[source_index]));

   RNBoolean *removed = new RNBoolean[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv)
      removed[iv] = FALSE;

   // variables to speed up algorithm
   const int nneighbors = 6;

   // visit vertices in increasing order
   int voxels_seen = 0;
   while (!voxel_heap.IsEmpty()) {
      DijkstraVoxelNode *current = voxel_heap.DeleteMin();

      NeuronVoxel *voxel = current->voxel;
      removed[voxel->DataIndex()] = TRUE;

      RNBoolean finished = TRUE;
      for (int ib = 0; ib < prediction->NBoundaries(); ++ib) {
         NeuronPredictionBoundary *boundary = prediction->Boundary(ib);

         // go through all neighbors
         NeuronPrediction *neighbor = boundary->OtherPrediction(prediction);
         int center_voxel_data_index = neighbor->CenterVoxel()->DataIndex();

         if (!removed[center_voxel_data_index]) finished = FALSE;
      }

      if (finished) break;

      // visit all of the neighbors of current voxel
      for (int iv = 0; iv < nneighbors; ++iv) {
         // get neighbor
         NeuronVoxel *neighbor = voxel->Neighbor(iv);
         if (!neighbor) continue;

         // get the neighbor data
         int neighbor_id = neighbor->DataIndex();
         DijkstraVoxelNode *neighbor_data = &(voxel_data[neighbor_id]);

         int ix, iy, iz;
         nd->IndexToIndices(voxel->DataIndex(), ix, iy, iz);
         int ii, ij, ik;
         nd->IndexToIndices(neighbor->DataIndex(), ii, ij, ik);

         RNScalar affinity = voxel->AffinityToNeighbor(neighbor);

         if (abs(ix - ii) == 1)
            affinity = scaling[RN_X] * (max_affinity - affinity);
         else if (abs(iy - ij) == 1)
            affinity = scaling[RN_Y] * (max_affinity - affinity);
         else if (abs(iz - ik) == 1)
            affinity = scaling[RN_Z] * (max_affinity - affinity);
         else
            rn_assertion(FALSE);

         if (affinity < 0.0) affinity = 0.0;
         RNScalar new_distance = current->distance + affinity;
         RNScalar old_distance = neighbor_data->distance;

         if (!neighbor_data->visited) {
            neighbor_data->prev = current;
            neighbor_data->distance = new_distance;
            neighbor_data->visited = TRUE;
            voxel_heap.Insert(neighbor_id, neighbor_data);
         }
         // update distance if needed
         else if (new_distance < old_distance) {
            neighbor_data->prev = current;
            neighbor_data->distance = new_distance;
            voxel_heap.DecreaseKey(neighbor_id, neighbor_data);
         }
      }

      // increment counter variable
      voxels_seen++;
   }

   std::vector<RNScalar> distances = std::vector<RNScalar>();
   std::vector<RNScalar> nsteps = std::vector<RNScalar>();
   std::vector<RNScalar> kurtosis = std::vector<RNScalar>();
   std::vector<RNScalar> maximum = std::vector<RNScalar>();
   std::vector<RNScalar> mean = std::vector<RNScalar>();
   std::vector<RNScalar> median = std::vector<RNScalar>();
   std::vector<RNScalar> minimum = std::vector<RNScalar>();
   std::vector<RNScalar> skew = std::vector<RNScalar>();
   std::vector<RNScalar> stddev = std::vector<RNScalar>();
   std::vector<RNScalar> first_quintile = std::vector<RNScalar>();
   std::vector<RNScalar> second_quintile = std::vector<RNScalar>();
   std::vector<RNScalar> third_quintile = std::vector<RNScalar>();
   std::vector<RNScalar> fourth_quintile = std::vector<RNScalar>();

   // go through all neighbors
   for (int ib = 0; ib < prediction->NBoundaries(); ++ib) {
      NeuronPredictionBoundary *boundary = prediction->Boundary(ib);

      NeuronPrediction *neighbor = boundary->OtherPrediction(prediction);

      rn_assertion(removed[neighbor->CenterVoxel()->DataIndex()]);

      int voxel_index = neighbor->CenterVoxel()->DataIndex();
      DijkstraVoxelNode *data = &voxel_data[voxel_index];

      RNScalar distance = data->distance;
      std::vector<RNScalar> path_steps = std::vector<RNScalar>();

      // recreate the dijkstra path
      while (data->prev != NULL) {
         // get the current voxel
         NeuronVoxel *current_voxel = data->voxel;

         // get the previous voxel
         DijkstraVoxelNode *previous_node = data->prev;
         NeuronVoxel *previous_voxel = previous_node->voxel;

         // get the affinity between the previous voxel and the current voxel
         RNScalar affinity = current_voxel->AffinityToNeighbor(previous_voxel);
         path_steps.push_back(affinity);

         // update the current voxel on the path
         data = data->prev;
      }
      RNDistribution distribution = RNDistribution(path_steps);

      distances.push_back(distance);
      nsteps.push_back(distribution.NPoints());
      kurtosis.push_back(distribution.Kurtosis());
      maximum.push_back(distribution.Maximum());
      mean.push_back(distribution.Mean());
      median.push_back(distribution.Median());
      minimum.push_back(distribution.Minimum());
      skew.push_back(distribution.Skew());
      stddev.push_back(distribution.StdDev());
      first_quintile.push_back(distribution.Quintile(1));
      second_quintile.push_back(distribution.Quintile(2));
      third_quintile.push_back(distribution.Quintile(3));
      fourth_quintile.push_back(distribution.Quintile(4));
   }

   // create output file
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_dijkstra_%04d.feature", tmp_postprocess_directory, root_filename, prediction_index);

   // open file
   FILE *output_fp = fopen(output_filename, "wb");
   if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // write the number of boundaries
   int nboundaries = prediction->NBoundaries();
   fwrite(&nboundaries, sizeof(int), 1, output_fp);

   // go through all neighbors
   for (int ib = 0; ib < prediction->NBoundaries(); ++ib) {
      // write all of the results
      fwrite(&(distances[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(nsteps[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(kurtosis[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(maximum[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(mean[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(median[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(minimum[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(skew[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(stddev[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(first_quintile[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(second_quintile[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(third_quintile[ib]), sizeof(RNScalar), 1, output_fp);
      fwrite(&(fourth_quintile[ib]), sizeof(RNScalar), 1, output_fp);

      int distance_ranking = Ranking(distances, ib);
      int nsteps_ranking = Ranking(nsteps, ib);
      int kurtosis_ranking = Ranking(kurtosis, ib);
      int maximum_ranking = Ranking(maximum, ib);
      int mean_ranking = Ranking(mean, ib);
      int median_ranking = Ranking(median, ib);
      int minimum_ranking = Ranking(minimum, ib);
      int skew_ranking = Ranking(skew, ib);
      int stddev_ranking = Ranking(stddev, ib);
      int first_quintile_ranking = Ranking(first_quintile, ib);
      int second_quintile_ranking = Ranking(second_quintile, ib);
      int third_quintile_ranking = Ranking(third_quintile, ib);
      int fourth_quintile_ranking = Ranking(fourth_quintile, ib);
      fwrite(&distance_ranking, sizeof(int), 1, output_fp);
      fwrite(&nsteps_ranking, sizeof(int), 1, output_fp);
      fwrite(&kurtosis_ranking, sizeof(int), 1, output_fp);
      fwrite(&maximum_ranking, sizeof(int), 1, output_fp);
      fwrite(&mean_ranking, sizeof(int), 1, output_fp);
      fwrite(&median_ranking, sizeof(int), 1, output_fp);
      fwrite(&minimum_ranking, sizeof(int), 1, output_fp);
      fwrite(&skew_ranking, sizeof(int), 1, output_fp);
      fwrite(&stddev_ranking, sizeof(int), 1, output_fp);
      fwrite(&first_quintile_ranking, sizeof(int), 1, output_fp);
      fwrite(&second_quintile_ranking, sizeof(int), 1, output_fp);
      fwrite(&third_quintile_ranking, sizeof(int), 1, output_fp);
      fwrite(&fourth_quintile_ranking, sizeof(int), 1, output_fp);
   }

   // close file
   fclose(output_fp);

   // free memory
   delete[] removed;
   delete[] voxel_data;

   // print statistics
   if (print_verbose) { printf("Completed dijkstra for prediction %d in %0.2f seconds\n", prediction_index, start_time.Elapsed()); }

   // return success
   return 1;
}



static int MergeResults(char root_filename[4096])
{
   // create array for all results
   float **dijkstra_features = new float *[NSTATISTICS];
   for (int id = 0; id < NSTATISTICS; ++id) {
      dijkstra_features[id] = new float[nd->NPredictionBoundaries()];
   }

   // go through all predictions
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      NeuronPrediction *prediction = nd->Prediction(ip);

      // get filename
      char input_filename[4096];
      sprintf(input_filename, "%s/%s_dijkstra_%04d.feature", tmp_postprocess_directory, root_filename, ip);

      // open file
      FILE *input_fp = fopen(input_filename, "rb");
      if (!input_fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

      int nboundaries;
      fread(&nboundaries, sizeof(int), 1, input_fp);
      rn_assertion(nboundaries == prediction->NBoundaries());

      for (int ib = 0; ib < prediction->NBoundaries(); ++ib) {
         NeuronPredictionBoundary *boundary = prediction->Boundary(ib);
         // read in all statistics
         for (int is = 0; is < NSTATISTICS; ++is) {
            // raw values
            if (is < NSTATISTICS / 2) {
               RNScalar attribute_value;
               if (fread(&attribute_value, sizeof(RNScalar), 1, input_fp) != 1) {
                  fprintf(stderr, "Failed to read %s\n", input_filename);
                  return 0;
               }

               // add to arrays
               dijkstra_features[is][boundary->DataIndex()] = (float)attribute_value;
            }
            // rankings
            else {
               int attribute_value;
               if (fread(&attribute_value, sizeof(int), 1, input_fp) != 1) {
                  fprintf(stderr, "Failed to read %s\n", input_filename); 
                  return 0;
               }

               // add to arrays
               dijkstra_features[is][boundary->DataIndex()] = (float)attribute_value;
            }
         }
      }

      // close file
      fclose(input_fp);
   }

   for (int is = 0; is < NSTATISTICS; ++is) {
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_dijkstra_%s.feature", postprocess_directory, root_filename, statistic_names[is]);

      // open file
      FILE *output_fp = fopen(output_filename, "wb");
      if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

      // write the total number of boundaries
      int nboundaries = nd->NPredictionBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, output_fp);

      // write the dikjstra feature
      fwrite(dijkstra_features[is], sizeof(float), nd->NPredictionBoundaries(), output_fp);

      // close file
      fclose(output_fp);
   }

   // free memory
   for (int id = 0; id < NSTATISTICS; ++id)
      delete[] dijkstra_features[id];
   delete[] dijkstra_features;

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
         else if (!strcmp(*argv, "-prediction_index")) { argv++; argc--; starting_prediction_index = atoi(*argv); }
         else if (!strcmp(*argv, "-nsteps")) { argv++; argc--; nsteps = atoi(*argv); }
         else if (!strcmp(*argv, "-merge")) { merge_results = 1; }
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
   if (!merge_results && starting_prediction_index == -1) {
      fprintf(stderr, "Must either merge results or choose a starting prediction index\n");
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

   // set scaling factors
   for (int dim = 0; dim <= 2; ++dim) {
      scaling[dim] = nd->Resolution(RN_X) / nd->Resolution(dim);
   }

   // read in prediction file
   char meta_root_filename[4096];
   sprintf(meta_root_filename, "output/ordermerge/%s", root_filename);

   prediction_grid = RNReadNeuronMetaRawFile(meta_root_filename);
   if (!prediction_grid) exit(-1);

   // get the number of predictions
   npredictions = (int)(prediction_grid->Maximum() + 1.5);

   // create boundary predictions
   nd->CreatePredictions(prediction_grid);

   // run dijkstra algorithm for nsteps at starting prediction index
   if (!merge_results) {
      for (int ip = starting_prediction_index; ip < starting_prediction_index + nsteps; ++ip) {
         if (ip >= nd->NPredictions()) continue;
         if (!Dijkstra(ip, root_filename)) exit(-1);
      }
   }
   else {
      if (!MergeResults(root_filename)) exit(-1);
   }

   // free memory
   delete nd;
   delete prediction_grid;

   // return success
   return 0;
}
