// Source file for the simple dijkstra algorithm



// include files 

#include "Neuron/Neuron.h"
#include "RNDataStructures/RNDataStructures.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static int cellular_index = -1;
static int nsteps = 5;
static int merge_results = 0;



// program variables

static NeuronData *nd;
static RNScalar max_affinity = 0.82;
static RNScalar scaling[3];




// useful directories

static const char *tmp_algs_directory = "algs_data/features/tmp";
static const char *features_directory = "algs_data/features";
static const char *distances_directory = "distances";
static const char *tmp_distances_directory = "distances/tmp";



static const int NFEATURES = 13;
static const char *feature_names[NFEATURES] = {
   "distances", "nsteps", "kurtosis", "maximum", "mean",
   "median", "minimum", "skew", "stddev", "first_quintile",
   "second_quintile", "third_quintile", "fourth_quintile"
};



////////////////////////////////////////////////////////////////////////
// Useful structs
////////////////////////////////////////////////////////////////////////

// struct for dijkstra voxel node
struct DijkstraVoxelNode
{
   NeuronVoxel *voxel;
   DijkstraVoxelNode *prev;
   RNScalar distance;
   RNBoolean visited;
};



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



////////////////////////////////////////////////////////////////////////
// Dijkstra algorithm for voxels
////////////////////////////////////////////////////////////////////////

static int
RunDijkstraVoxels(int source_index, float *distances, int *nsteps, float *kurtosis, float *maximum, float *mean, float *median, float *minimum, float *skew, float *stddev, float *first_quintile, float *second_quintile, float *third_quintile, float *fourth_quintile, FILE *boundary_fp)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

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

   // insert the source into the heap
   voxel_data[source_index].voxel = nd->Voxel(source_index);
   voxel_data[source_index].prev = NULL;
   voxel_data[source_index].distance = 0.0;
   voxel_data[source_index].visited = TRUE;
   voxel_heap.Insert(source_index, &(voxel_data[source_index]));

   // variables to speed up algorithm
   const int nneighbors = 6;

   // visit vertices in increasing order
   int voxels_seen = 0;
   while (!voxel_heap.IsEmpty()) {
      DijkstraVoxelNode *current = voxel_heap.DeleteMin();

      NeuronVoxel *voxel = current->voxel;

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

   // save the information for all of the cellulars
   NeuronCellular *cellular_one = nd->Cellular(cellular_index);
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular_two = nd->Cellular(ic);
      
      // get voxel index
      int voxel_index = cellular_two->CenterVoxel()->DataIndex();
      DijkstraVoxelNode *data = &voxel_data[voxel_index];

      // create a distribution of all steps on the dijkstra path
      RNScalar distance = data->distance;
      std::vector<RNScalar> path_steps = std::vector<RNScalar>();

      // get the data indices
      int cellular_one_data_index = cellular_one->DataIndex();
      int cellular_two_data_index = cellular_two->DataIndex();

      // keep track of which cellular this voxel is in
      NeuronSupervoxel *current_supervoxel = NULL;
      NeuronSupervoxel *previous_supervoxel = cellular_two;

      // create supervoxel path
      std::vector<int> cellular_path = std::vector<int>();
      RNBoolean *cellular_path_hash = new RNBoolean[nd->NCellulars()];
      for (int icp = 0; icp < nd->NCellulars(); ++icp)
         cellular_path_hash[icp] = FALSE;

      // add start and end cellulars to path
      cellular_path_hash[cellular_one_data_index] = TRUE;
      cellular_path.push_back(cellular_one_data_index);
      cellular_path_hash[cellular_two_data_index] = TRUE;
      cellular_path.push_back(cellular_two_data_index);

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

         // if the current supervoxel doesn't equal previous supervoxel, add to path
         current_supervoxel = current_voxel->Supervoxel();
         if (current_supervoxel != previous_supervoxel) {
            // only consider supervoxels not on the stack
            if (current_supervoxel->IsCellular()) {
               if (!cellular_path_hash[current_supervoxel->DataIndex()]) {
                  cellular_path_hash[current_supervoxel->DataIndex()] = TRUE;
                  cellular_path.push_back(current_supervoxel->DataIndex());
               }
            }
         }

         // update the previous supervoxel
         previous_supervoxel = current_voxel->Supervoxel();

         // update the current voxel on the path
         data = data->prev;
      }
      RNDistribution distribution = RNDistribution(path_steps);

      distances[ic] = distance;
      nsteps[ic] = distribution.NPoints();
      kurtosis[ic] = distribution.Kurtosis();
      maximum[ic] = distribution.Maximum();
      mean[ic] = distribution.Mean();
      median[ic] = distribution.Median();
      minimum[ic] = distribution.Minimum();
      skew[ic] = distribution.Skew();
      stddev[ic] = distribution.StdDev();
      first_quintile[ic] = distribution.Quintile(1);
      second_quintile[ic] = distribution.Quintile(2);
      third_quintile[ic] = distribution.Quintile(3);
      fourth_quintile[ic] = distribution.Quintile(4);

      // write the cellular start and end index
      int start_index = cellular_index;
      int end_index = ic;
      fwrite(&start_index, sizeof(int), 1, boundary_fp);
      fwrite(&end_index, sizeof(int), 1, boundary_fp);

      // write the number of cellulars on the path
      int ncellulars_on_path = cellular_path.size();
      fwrite(&ncellulars_on_path, sizeof(int), 1, boundary_fp);

      // write all of the cellulars
      for (int icp = 0; icp < ncellulars_on_path; ++icp) {
         int element_cellular_index = cellular_path[icp];
         fwrite(&element_cellular_index, sizeof(int), 1, boundary_fp);
      }

      // free memory
      delete[] cellular_path_hash;
   }

   // remove voxel data
   delete[] voxel_data;

   // print statistics
   if (print_verbose) printf("Ran dijkstra algorithm in %.2f seconds\n", start_time.Elapsed());

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Dijkstra algorithm for single iteration or merging
////////////////////////////////////////////////////////////////////////


static int
RunDijkstraSingleIteration(char root_filename[4096])
{
   // get output filename for boundary paths
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_%06d.path", tmp_algs_directory, root_filename, cellular_index);

   // open file
   FILE *boundary_fp = fopen(output_filename, "wb");
   if (!boundary_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   int ncellulars = nd->NCellulars();
   fwrite(&ncellulars, sizeof(int), 1, boundary_fp);

   // get the center voxel for this cellular
   NeuronCellular *cellular = nd->Cellular(cellular_index);
   int voxel_index = cellular->CenterVoxel()->DataIndex();

   float *distances = new float[ncellulars];
   int *nsteps = new int[ncellulars];
   float *kurtosis = new float[ncellulars];
   float *maximum = new float[ncellulars];
   float *mean = new float[ncellulars];
   float *median = new float[ncellulars];
   float *minimum = new float[ncellulars];
   float *skew = new float[ncellulars];
   float *stddev = new float[ncellulars];
   float *first_quintile = new float[ncellulars];
   float *second_quintile = new float[ncellulars];
   float *third_quintile = new float[ncellulars];
   float *fourth_quintile = new float[ncellulars];

   // run the actual algorithm
   if (!RunDijkstraVoxels(voxel_index, distances, nsteps, kurtosis, maximum, mean, median, minimum, skew, stddev, first_quintile, second_quintile, third_quintile, fourth_quintile, boundary_fp)) { fprintf(stderr, "Failed to run dijkstra's algorithm on voxels\n"); return 0; }

   // create output filename
   char output_dijkstra_filename[4096];
   sprintf(output_dijkstra_filename, "%s/%s_dijkstra_%05d.distance", tmp_distances_directory, root_filename, cellular_index);

   // write file
   FILE *dijkstra_fp = fopen(output_dijkstra_filename, "wb");
   if (!dijkstra_fp) { fprintf(stderr, "Failed to write to %s\n", output_dijkstra_filename); return 0; }

   // save the number of cellulars
   fwrite(&ncellulars, sizeof(int), 1, dijkstra_fp);

   // save the information stored in path attributes
   fwrite(distances, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(nsteps, sizeof(int), ncellulars, dijkstra_fp);
   fwrite(kurtosis, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(maximum, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(mean, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(median, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(minimum, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(skew, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(stddev, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(first_quintile, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(second_quintile, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(third_quintile, sizeof(float), ncellulars, dijkstra_fp);
   fwrite(fourth_quintile, sizeof(float), ncellulars, dijkstra_fp);
   
   // close files
   fclose(dijkstra_fp);
   fclose(boundary_fp);

   // free memory
   delete[] distances;
   delete[] nsteps;
   delete[] kurtosis;
   delete[] maximum;
   delete[] mean;
   delete[] median;
   delete[] minimum;
   delete[] skew;
   delete[] stddev;
   delete[] first_quintile;
   delete[] second_quintile;
   delete[] third_quintile;
   delete[] fourth_quintile;

   // return success 
   return 1;
}



static int
MergeAllCellulars(char root_filename[4096])
{
   // create all filenames
   char filename[NFEATURES][4096];
   for (int i = 0; i < NFEATURES; ++i) {
      if (i == 0) sprintf(filename[i], "%s/%s_dijkstra.distance", distances_directory, root_filename);
      else sprintf(filename[i], "%s/%s_dijkstra_%s.feature", features_directory, root_filename, feature_names[i]);
   }

   // open all files
   FILE *fp[NFEATURES];
   for (int i = 0; i < NFEATURES; ++i) {
      fp[i] = fopen(filename[i], "wb");
      if (!fp[i]) {
         fprintf(stderr, "Failed to write to %s\n", filename[i]);
         for (int j = 0; j < i; ++j) {
            fclose(fp[j]);
         }
         return 0;
      }
   }

   // write the number of cellulars
   int ncellulars = nd->NCellulars();
   for (int i = 0; i < NFEATURES; ++i)
      fwrite(&ncellulars, sizeof(int), 1, fp[i]);

   // open all cellular files
   if (print_verbose) { printf("Merging all cellulars...\n  "); fflush(stdout); }
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      if (print_verbose) RNProgressBar(ic, nd->NCellulars());

      char cellular_filename[4096];
      sprintf(cellular_filename, "%s/%s_dijkstra_%05d.distance", tmp_distances_directory, root_filename, ic);

      // open file
      FILE *cellular_fp = fopen(cellular_filename, "rb");
      if (!cellular_fp) {
         fprintf(stderr, "Failed to read %s\n", cellular_filename);
         for (int j = 0; j < NFEATURES; ++j)
            fclose(fp[j]);
         return 0;
      }

      // get the number of input cellulars
      int ninput_ncellulars;
      fread(&ninput_ncellulars, sizeof(int), 1, cellular_fp);
      rn_assertion(ninput_ncellulars == nd->NCellulars());

      // read in all of the data
      for (int j = 0; j < NFEATURES; ++j) {
         if (j == 1) {
            // read the results
            int *results = new int[ninput_ncellulars];
            if (fread(results, sizeof(int), ninput_ncellulars, cellular_fp) != (unsigned int)ninput_ncellulars) {
               fprintf(stderr, "Corrupted file: %s\n", cellular_filename);
               return 0;
            }

            // write the results to the merge file
            fwrite(results, sizeof(int), ninput_ncellulars, fp[j]);

            // free memory
            delete[] results;
         }
         else {
            // read the results
            float *results = new float[ninput_ncellulars];
            if (fread(results, sizeof(float), ninput_ncellulars, cellular_fp) != (unsigned int)ninput_ncellulars) {
               fprintf(stderr, "Corrupted file: %s\n", cellular_filename);
               return 0;
            }
	    
            // write the result to the merge file
            fwrite(results, sizeof(float), ninput_ncellulars, fp[j]);

            // free memory
            delete[] results;
         }
      }

      // close file
      fclose(cellular_fp);
   }
   if (print_verbose) printf("\ndone!\n");

   // close all files
   for (int j = 0; j < NFEATURES; ++j) {
      fclose(fp[j]);
   }

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
         else if (!strcmp(*argv, "-max_affinity")) { argv++; argc--; max_affinity = ((RNScalar)atoi(*argv)) / 100; }
         else if (!strcmp(*argv, "-cellular_index")) { argv++; argc--; cellular_index = atoi(*argv); }
         else if (!strcmp(*argv, "-nsteps")) { argv++; argc--; nsteps = atoi(*argv); }
         else if (!strcmp(*argv, "-merge")) { merge_results = 1; }
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input argument
   if (!input_filename) { fprintf(stderr, "Need to supply input argument\n"); return 0; }

   // make sure a cellular index is selected 
   if (cellular_index == -1 && !merge_results) { fprintf(stderr, "Need to specify the source index or merge existing results\n"); return 0; }

   // Return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
   // Parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in the neuron data
   if (!ReadData(input_filename)) exit(-1);

   // set scaling factors
   for (int dim = 0; dim <= 2; ++dim) {
      scaling[dim] = nd->Resolution(RN_X) / nd->Resolution(dim);
   }

   // find root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *endp = strrchr(root_filename, '.');
   *endp = '\0';

   // run on one cellular or merge results
   if (!merge_results) {
      // read all voxels
      nd->ReadVoxels();

      // perform 20 iterations
      for (int ic = 0; ic < nsteps; ++ic) {
         if (cellular_index >= nd->NCellulars()) break;

         if (!RunDijkstraSingleIteration(root_filename)) exit(-1);
         cellular_index++;
      }

      // release all voxels
      nd->ReleaseVoxels();

   }
   if (merge_results && !MergeAllCellulars(root_filename)) exit(-1);

   // free memory
   delete nd;

   // Return success
   return 0;
}
