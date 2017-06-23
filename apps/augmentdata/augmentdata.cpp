// Source file for the neuron statistics algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int affinities = 0;
static int dijkstra = 0;
static int medians = 0;
static int randomized = 0;



// global variables

static NeuronData *nd = NULL;
static char *input_filename = NULL;
static RNScalar max_affinity = 0.82;
static RNScalar scaling[3];


// directory structure

static const char *results_directory = "results/miscellaneous";
static const char *graphcut_directory = "algs_data/graphcut";
static const char *boundaries_directory = "boundaries";
static const char *monte_carlo_directory = "algs_data/montecarlo";



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
      printf("  Bounding Box: (%0.2f, %0.2f, %0.2f) to (%0.2f, %0.2f, %0.2f)\n", nd->WorldBox().XMin(), nd->WorldBox().YMin(), nd->WorldBox().ZMin(), nd->WorldBox().XMax(), nd->WorldBox().YMax(), nd->WorldBox().ZMax());
      printf("  Voxels: %d\n", nd->NVoxels());
      printf("  Supervoxels: %d\n", nd->NSupervoxels());
      printf("  Extracellulars: %d\n", nd->NExtracellulars());
      printf("  Boundaries: %d\n", nd->NBoundaries());
      printf("  Human Labels: %d\n", nd->NHumanLabels());
      printf("  Predictions: %d\n", nd->NPredictions());
      fflush(stdout);
   }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// General statistic functions
////////////////////////////////////////////////////////////////////////

int BoundaryMaxSort(NeuronBoundary *one, NeuronBoundary *two)
{
   return one->Maximum() > two->Maximum();
}



int BoundaryMeanSort(NeuronBoundary *one, NeuronBoundary *two)
{
   return one->Mean() > two->Mean();
}



static int AffinityPlots(char root_filename[4096])
{
   int noversegmented_boundaries = 0;

   // add all boundaries to the vector
   std::vector<NeuronBoundary *> boundary_maxes = std::vector<NeuronBoundary *>();
   std::vector<NeuronBoundary *> boundary_means = std::vector<NeuronBoundary *>();
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // skip extracellulars
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
      if (supervoxel_one->IsExtracellular()) continue;
      if (supervoxel_two->IsExtracellular()) continue;

      // add to vector
      boundary_maxes.push_back(boundary);
      boundary_means.push_back(boundary);

      if (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel())
         noversegmented_boundaries++;
   }

   // sort based on boundary mean
   sort(boundary_maxes.begin(), boundary_maxes.end(), BoundaryMaxSort);
   sort(boundary_means.begin(), boundary_means.end(), BoundaryMeanSort);

   char output_filename[4096];
   sprintf(output_filename, "%s/%s_affinities.results", results_directory, root_filename);

   FILE *output_fp = fopen(output_filename, "wb");
   if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // write the number of oversegmented boundaries
   fwrite(&noversegmented_boundaries, sizeof(int), 1, output_fp);

   // calculate precision and recall for max and then mean
   int ncorrect_merges = 0;
   for (unsigned int ib = 0; ib < boundary_maxes.size(); ++ib) {
      NeuronBoundary *boundary = boundary_maxes[ib];

      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      if (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel()) {
         ncorrect_merges++;

         // output the precision and recall
         RNScalar precision = ncorrect_merges / (RNScalar)(ib + 1);
         RNScalar recall = ncorrect_merges / (RNScalar)noversegmented_boundaries;

         fwrite(&precision, sizeof(RNScalar), 1, output_fp);
         fwrite(&recall, sizeof(RNScalar), 1, output_fp);
      }
   }
   ncorrect_merges = 0;
   for (unsigned int ib = 0; ib < boundary_means.size(); ++ib) {
      NeuronBoundary *boundary = boundary_means[ib];

      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      if (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel()) {
         ncorrect_merges++;

         // output the precision and recall
         RNScalar precision = ncorrect_merges / (RNScalar)(ib + 1);
         RNScalar recall = ncorrect_merges / (RNScalar)noversegmented_boundaries;

         fwrite(&precision, sizeof(RNScalar), 1, output_fp);
         fwrite(&recall, sizeof(RNScalar), 1, output_fp);
      }
   }

   // close file
   fclose(output_fp);

   // return success
   return 1;
}



// struct for dijkstra voxel node
struct DijkstraVoxelNode
{
   NeuronVoxel *voxel;
   DijkstraVoxelNode *prev;
   RNScalar distance;
   RNBoolean visited;
};



static int
RunDijkstraVoxels(int source_cellular_index, int target_cellular_index, FILE *output_fp)
{
   // save the information for all of the cellulars
   NeuronCellular *cellular_one = nd->Cellular(source_cellular_index);
   NeuronCellular *cellular_two = nd->Cellular(target_cellular_index);

   // get source and target voxel indices
   int source_index = cellular_one->CenterVoxel()->DataIndex();
   int target_index = cellular_two->CenterVoxel()->DataIndex();

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

      // only go to target index
      if (voxel->DataIndex() == target_index) break;

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

   // get voxel index
   int voxel_index = cellular_two->CenterVoxel()->DataIndex();
   DijkstraVoxelNode *data = &voxel_data[voxel_index];

   // create a distribution of all steps on the dijkstra path
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

   int nsteps = (int)(path_steps.size());
   fwrite(&nsteps, sizeof(int), 1, output_fp);

   for (int is = 0; is < (int)path_steps.size(); ++is) {
      RNScalar step = path_steps[is];
      fwrite(&step, sizeof(RNScalar), 1, output_fp);
   }

   RNDistribution distribution = RNDistribution(path_steps);
   printf("Distance: %lf\n", distance);
   printf("NSteps: %d\n", distribution.NPoints());
   printf("Kurtosis: %lf\n", distribution.Kurtosis());
   printf("Maximum: %lf\n", distribution.Maximum());
   printf("Mean: %lf\n", distribution.Mean());
   printf("Median: %lf\n", distribution.Median());
   printf("Minimum: %lf\n", distribution.Minimum());
   printf("Skew: %lf\n", distribution.Skew());
   printf("StdDev: %lf\n", distribution.StdDev());
   printf("First Quintile: %lf\n", distribution.Quintile(1));
   printf("Second Quintile: %lf\n", distribution.Quintile(2));
   printf("Third Quintile: %lf\n", distribution.Quintile(3));
   printf("Fouth Quintile: %lf\n", distribution.Quintile(4));

   // remove voxel data
   delete[] voxel_data;

   // print statistics
   if (print_verbose) printf("Ran dijkstra algorithm in %.2f seconds\n", start_time.Elapsed());

   // return success
   return 1;
}



static int DijkstraPlot(char root_filename[4096])
{
   // read smoothing file
   char smoothing_filename[4096];
   sprintf(smoothing_filename, "%s/%s_random_forest_smoothing.binary", graphcut_directory, root_filename);

   FILE *smoothing_fp = fopen(smoothing_filename, "rb");
   if (!smoothing_fp) { fprintf(stderr, "Failed to read %s\n", smoothing_filename); return  0; }

   int nboundaries;
   fread(&nboundaries, sizeof(int), 1, smoothing_fp);
   rn_assertion(nboundaries == nd->NBoundaries());

   float *smoothing_results = new float[nboundaries];
   fread(smoothing_results, sizeof(float), nboundaries, smoothing_fp);

   // close file
   fclose(smoothing_fp);

   // read voxels
   nd->ReadVoxels();

   // set scaling factors
   for (int dim = 0; dim <= 2; ++dim) {
      scaling[dim] = nd->Resolution(RN_X) / nd->Resolution(dim);
   }

   // choosen randomly
   int boundary_one_index = -1;
   int boundary_two_index = -1;

   for (int ic = 1761; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);

      int same_neuron = -1;
      int diff_neuron = -1;

      for (int ib = 0; ib < cellular->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = cellular->Boundary(ib);
         NeuronSupervoxel *neighbor = boundary->OtherSupervoxel(cellular);
         if (neighbor->IsExtracellular()) continue;

         int smoothing_result = (int)(smoothing_results[boundary->DataIndex()]);

         if (neighbor->MajorityHumanLabel() == cellular->MajorityHumanLabel()) {
            if (!smoothing_result) continue;
            same_neuron = boundary->DataIndex();
         }
         else {
            if (smoothing_result) continue;
            if (boundary->Mean() < 0.70) continue;
            diff_neuron = boundary->DataIndex();
         }
      }
      if ((same_neuron != -1) && (diff_neuron != -1)) {
         boundary_one_index = same_neuron;
         boundary_two_index = diff_neuron;
         break;
      }
   }

   if (boundary_one_index == -1 || boundary_two_index == -1) rn_assertion(FALSE);

   // choose two random boundaries
   NeuronBoundary *boundary_one = nd->Boundary(boundary_one_index);
   NeuronBoundary *boundary_two = nd->Boundary(boundary_two_index);

   printf("%lf %lf\n", boundary_one->Maximum(), boundary_one->Mean());
   printf("%lf %lf\n", boundary_two->Maximum(), boundary_two->Mean());

   NeuronSupervoxel *supervoxel_one = boundary_one->SupervoxelOne();
   NeuronSupervoxel *supervoxel_two = boundary_one->SupervoxelTwo();
   NeuronSupervoxel *supervoxel_three = boundary_two->SupervoxelOne();
   if (supervoxel_three == supervoxel_one || supervoxel_three == supervoxel_two)
      supervoxel_three = boundary_two->SupervoxelTwo();
   printf("%d: (%d, %d, %d)\n", supervoxel_one->DataIndex(), supervoxel_one->CenterVoxel()->XCoordinate(), supervoxel_one->CenterVoxel()->YCoordinate(), supervoxel_one->CenterVoxel()->ZCoordinate());
   printf("%d: (%d, %d, %d)\n", supervoxel_two->DataIndex(), supervoxel_two->CenterVoxel()->XCoordinate(), supervoxel_two->CenterVoxel()->YCoordinate(), supervoxel_two->CenterVoxel()->ZCoordinate());
   printf("%d: (%d, %d, %d)\n", supervoxel_three->DataIndex(), supervoxel_three->CenterVoxel()->XCoordinate(), supervoxel_three->CenterVoxel()->YCoordinate(), supervoxel_three->CenterVoxel()->ZCoordinate());

   // get output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_dijkstra_distributions.results", results_directory, root_filename);

   // open file
   FILE *output_fp = fopen(output_filename, "wb");
   if (!output_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // call dijkstra 
   if (!RunDijkstraVoxels(boundary_one->SupervoxelOne()->DataIndex(), boundary_one->SupervoxelTwo()->DataIndex(), output_fp)) return 0;
   if (!RunDijkstraVoxels(boundary_two->SupervoxelOne()->DataIndex(), boundary_two->SupervoxelTwo()->DataIndex(), output_fp)) return 0;

   // close file
   fclose(output_fp);

   // free memory
   delete[] smoothing_results;

   // return success
   return 1;
}



static int MediansPlot(char root_filename[4096])
{
   // read all voxels
   nd->ReadVoxels();

   static const int NKs = 2;
   static const char *extensions[NKs] = { "facility", "kmedian" };
   static int ks[NKs] = { 25, 600 };
   for (int is = 0; is < NKs; ++is) {
      for (int ik = 0; ik < NKs; ++ik) {
         int K = ks[ik];

         // read boundary file
         char boundary_filename[4096];
         sprintf(boundary_filename, "%s/%s_%s_%04d.boundary", boundaries_directory, root_filename, extensions[is], K);

         // open file
         FILE *fp = fopen(boundary_filename, "rb");
         if (!fp) { fprintf(stderr, "Failed to read %s\n", boundary_filename); return 0; }

         // read the number of boundaries
         int nboundaries;
         fread(&nboundaries, sizeof(int), 1, fp);
         rn_assertion(nboundaries == nd->NBoundaries());

         // read all boundaries
         RNBoolean *boundary_merges = new RNBoolean[nd->NBoundaries()];
         if (fread(boundary_merges, sizeof(RNBoolean), nd->NBoundaries(), fp) != (unsigned int)nd->NBoundaries()) {
            fprintf(stderr, "Failed to read %s\n", boundary_filename);
            return 0;
         }

         // merge clusters
         int *cellular_clusters = new int[nd->NCellulars()];
         for (int ic = 0; ic < nd->NCellulars(); ++ic)
            cellular_clusters[ic] = ic;

         if (print_verbose) { printf("Merging boundaries for K = %d...\n  ", K); fflush(stdout); }
         for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
            if (print_verbose) RNProgressBar(ib, nd->NBoundaries());
            NeuronBoundary *boundary = nd->Boundary(ib);
            if (boundary_merges[ib]) {
               NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
               NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
               rn_assertion(supervoxel_one->IsCellular());
               rn_assertion(supervoxel_two->IsCellular());

               int cellular_one_label = cellular_clusters[supervoxel_one->DataIndex()];
               int cellular_two_label = cellular_clusters[supervoxel_two->DataIndex()];

               if (cellular_one_label == cellular_two_label) continue;

               // update all cellulars
               for (int ic = 0; ic < nd->NCellulars(); ++ic) {
                  if (cellular_clusters[ic] == cellular_two_label) {
                     cellular_clusters[ic] = cellular_one_label;
                  }
               }
            }
         }
         if (print_verbose) printf("\ndone.\n");

         // deflate integer array
         RNDeflateIntegerArray(cellular_clusters, nd->NCellulars());

         // create output grid
         char output_filename[4096];
         sprintf(output_filename, "%s/%s_%s_%04d", results_directory, root_filename, extensions[is], K);

         RNMeta meta = RNMeta("Uint16", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 1);
         R3Grid *grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
         for (int iv = 0; iv < nd->NVoxels(); ++iv) {
            NeuronVoxel *voxel = nd->Voxel(iv);
            if (voxel->IsExtracellular()) grid->SetGridValue(iv, 0);
            else grid->SetGridValue(iv, cellular_clusters[voxel->Supervoxel()->DataIndex()] + 1);
         }

         // write meta grid file
         if (!RNWriteNeuronMetaRawFile(output_filename, meta, grid)) return 0;

         // free memory
         delete[] boundary_merges;
         delete[] cellular_clusters;

         // close file
         fclose(fp);
      }
   }

   // return success
   return 1;
}



static int RandomizedPlot(char root_filename[4096])
{
   // get sum total 
   unsigned int *total_crossings = new unsigned int[nd->NBoundaries()];
   for (int ib = 0; ib < nd->NBoundaries(); ++ib)
      total_crossings[ib] = 0;

   float **ratios = new float *[200];
   for (int it = 0; it < 200; ++it)
      ratios[it] = new float[nd->NBoundaries()];

   

   // read all monte carlo files
   for (int it = 0; it < 200; ++it) {
      char input_filename[4096];
      sprintf(input_filename, "%s/%s_monte_carlo%03d.feature", monte_carlo_directory, root_filename, it);

      // open file
      FILE *input_fp = fopen(input_filename, "rb");
      if (!input_fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

      // read the number of boundaries
      int nboundaries;
      fread(&nboundaries, sizeof(int), 1, input_fp);
      rn_assertion(nboundaries == nd->NBoundaries());

      // read in all boundaries
      unsigned int *boundary_crossings = new unsigned int[nd->NBoundaries()];
      if (fread(boundary_crossings, sizeof(unsigned int), nd->NBoundaries(), input_fp) != (unsigned int)nd->NBoundaries()) {
         fprintf(stderr, "Failed to read %s\n", input_filename);
         return 0;
      }

      // update total crossing
      for (int ib = 0; ib < nd->NBoundaries(); ++ib)
         total_crossings[ib] += boundary_crossings[ib];

      // get current ratio
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         float ratio = total_crossings[ib] / (RNScalar)(it + 1);
         ratios[it][ib] = ratio / nd->Boundary(ib)->NAffinities();
      }

      // close file
      fclose(input_fp);

      // free memory
      delete[] boundary_crossings;
   }

   for (int it = 0; it < 200; ++it) {
      RNScalar difference = 0.0;
      // go through all boundaries
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         difference += sqrt((ratios[it][ib] - ratios[199][ib]) * (ratios[it][ib] - ratios[199][ib]));
      }
      printf("%3d: %lf\n", it, difference / nd->NBoundaries());
   }

   // free memory
   delete[] total_crossings;
   for (int it = 0; it < 200; ++it)
      delete[] ratios[it];
   delete[] ratios;

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
         else if (!strcmp(*argv, "-affinities")) affinities = 1;
         else if (!strcmp(*argv, "-dijkstra")) dijkstra = 1;
         else if (!strcmp(*argv, "-medians")) medians = 1;
         else if (!strcmp(*argv, "-randomized")) randomized = 1;
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

int
main(int argc, char **argv)
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

   ///////////////////////////////////
   //// RUN ALL OF THE STATISTICS ////
   ///////////////////////////////////

   if (affinities && !AffinityPlots(root_filename)) exit(-1);
   if (dijkstra && !DijkstraPlot(root_filename)) exit(-1);
   if (medians && !MediansPlot(root_filename)) exit(-1);
   if (randomized && !RandomizedPlot(root_filename)) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}