// Source file for the simple dijkstra algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;



// program variables

static NeuronData *nd;
static RNScalar max_affinity = 1.0;
static RNScalar scaling[3];



// useful directories

static const char *distances_directory = "distances";



// global variables

static float *distances = NULL;
static int *prev_node = NULL;



////////////////////////////////////////////////////////////////////////
// Useful structs
////////////////////////////////////////////////////////////////////////

// struct for dijkstra voxel node
struct DijkstraNode
{
   NeuronVoxel *voxel;
   DijkstraNode *prev;
   RNScalar distance;
   RNBoolean visited;
   int supervoxel_index;
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
// Dijkstra algorithm for voxels
////////////////////////////////////////////////////////////////////////

static int
RunDijkstraAlgorithm(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate temporary data
   DijkstraNode *data = new DijkstraNode[nd->NVoxels()];
   if (!data) {
      fprintf(stderr, "Unable to allocate temporary data for geodisic distance.\n");
      return 0;
   }

   DijkstraNode tmp;
   RNMinBinaryHeap<DijkstraNode *> heap(&tmp, &(tmp.distance), nd->NVoxels());

   // initialize all data
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      data[iv].voxel = voxel;
      data[iv].prev = NULL;
      data[iv].distance = FLT_MAX;
      data[iv].visited = FALSE;
      data[iv].supervoxel_index = -1;
   }

   // all cellulars start as sources
   for (int is = 0; is < nd->NSupervoxels(); ++is) {
      int voxel_index = nd->Supervoxel(is)->CenterVoxel()->DataIndex();
      data[voxel_index].distance = 0.0;
      data[voxel_index].visited = TRUE;
      data[voxel_index].supervoxel_index = is;
      heap.Insert(voxel_index, &(data[voxel_index]));
   }

   // visit vertices in increasing order
   int voxels_seen = 0;
   if (print_verbose) { printf("Run dijkstra algorithm to boundary...\n  "); fflush(stdout); }
   while (!heap.IsEmpty()) {
      if (print_verbose) RNProgressBar(voxels_seen, nd->NVoxels());
      DijkstraNode *current = heap.DeleteMin();

      NeuronVoxel *voxel = current->voxel;

      // visit all of the neighbors of current voxel
      for (int iv = 0; iv < voxel->NNeighbors(); ++iv) {
         // get neighbor
         NeuronVoxel *neighbor = voxel->Neighbor(iv);
         if (!neighbor) continue;

         // get the neighbor data
         int neighbor_id = neighbor->DataIndex();
         DijkstraNode *neighbor_data = &(data[neighbor_id]);

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
            neighbor_data->supervoxel_index = current->supervoxel_index;
            heap.Insert(neighbor_id, neighbor_data);
         }
         else if (new_distance < old_distance) {
            neighbor_data->prev = current;
            neighbor_data->distance = new_distance;
            neighbor_data->supervoxel_index = current->supervoxel_index;
            heap.DecreaseKey(neighbor_id, neighbor_data);
         }
      }

      // incrememt counter variables
      voxels_seen++;
   }

   // how many voxels belong to different cellulars
   int ncellular_successes = 0;
   int ncellular_failures = 0;
   int nprevious_successes = 0;
   int nprevious_failures = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (!voxel->HumanLabel()) continue;

      int old_supervoxel_index = voxel->Supervoxel()->DataIndex();
      int new_supervoxel_index = data[iv].supervoxel_index;

      NeuronSupervoxel *new_supervoxel = nd->Supervoxel(new_supervoxel_index);
      NeuronSupervoxel *old_supervoxel = nd->Supervoxel(old_supervoxel_index);
      if (voxel->HumanLabel() == new_supervoxel->MajorityHumanLabel()) ncellular_successes++;
      else ncellular_failures++;
      if (voxel->HumanLabel() == old_supervoxel->MajorityHumanLabel()) nprevious_successes++;
      else nprevious_failures++;
   }

   // print out stats
   printf("------------------------------------------------------\n");
   printf("Voxel Statistics\n");
   printf("------------------------------------------------------\n");
   printf("Voxel Cellular Assignment Success: %d (%0.2f%%)\n", ncellular_successes, 100.0 * ncellular_successes / (RNScalar)(ncellular_successes + ncellular_failures));
   printf("Voxel Cellular Assignment Failure: %d (%0.2f%%)\n", ncellular_failures, 100.0 * ncellular_failures / (RNScalar)(ncellular_successes + ncellular_failures));
   printf("Prev Voxel Cellular Assignment Success: %d (%0.2f%%)\n", nprevious_successes, 100.0 * nprevious_successes / (RNScalar)(nprevious_successes + nprevious_failures));
   printf("Prev Voxel Cellular Assignment Failure: %d (%0.2f%%)\n", nprevious_failures, 100.0 * nprevious_failures / (RNScalar)(nprevious_successes + nprevious_failures));
   printf("------------------------------------------------------\n");


   // free memory
   delete[] data;

   // print statistics
   if (print_verbose) {
      printf("\nGenerated dijkstra distances...\n");
      printf("  Time = %0.2f seconds\n", start_time.Elapsed());
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
         else if (!strcmp(*argv, "-max_affinity")) { argv++; argc--; max_affinity = atof(*argv); }
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

   // allocate memory for distance arrays
   distances = new float[nd->NVoxels()];
   prev_node = new int[nd->NVoxels()];

   // run dijkstra algorithm
   if (!RunDijkstraAlgorithm()) return 0;

   // save distance and prev node files
   char distances_filename[4096];
   sprintf(distances_filename, "%s/%s_dijkstra_distances_boundary.distance", distances_directory, root_filename);

   // open file
   FILE *distances_fp = fopen(distances_filename, "wb");
   if (!distances_fp) { fprintf(stderr, "Failed to write to %s\n", distances_filename); exit(-1); }

   // write file
   int nvoxels = nd->NVoxels();
   fwrite(&nvoxels, sizeof(int), 1, distances_fp);
   fwrite(distances, sizeof(float), nvoxels, distances_fp);

   // close file
   fclose(distances_fp);

   char prev_node_filename[4096];
   sprintf(prev_node_filename, "%s/%s_dijkstra_prev_node_boundary.distance", distances_directory, root_filename);

   // open file
   FILE *prev_fp = fopen(prev_node_filename, "wb");
   if (!prev_fp) { fprintf(stderr, "Failed to write to %s\n", prev_node_filename); exit(-1); }

   // write file
   fwrite(&nvoxels, sizeof(int), 1, prev_fp);
   fwrite(prev_node, sizeof(int), nvoxels, prev_fp);

   // close file
   fclose(prev_fp);

   // free memory
   delete nd;
   delete[] distances;
   delete[] prev_node;

   // Return success
   return 0;
}
