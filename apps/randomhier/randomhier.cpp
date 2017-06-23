// Source file to create graph cut initial positions



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static const char *extension = "random_forest_smoothing_global";
static int threshold = 50;
static int validation = 0;


// program variables

static NeuronData *nd;
static float *smoothing_terms = NULL;



// directory variables

static const char *hierarchical_directory = "algs_data/hierarchical";
static const char *graphcut_directory = "algs_data/graphcut";



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
// Hierarchical functions
////////////////////////////////////////////////////////////////////////

static int
HierarchicalClustering(const char root_filename[4096])
{
   // get merge filename
   char merge_filename[4096];
   sprintf(merge_filename, "%s/%s_%s.hier", hierarchical_directory, root_filename, extension);

   // open file
   FILE *merge_fp = fopen(merge_filename, "rb");
   if (!merge_fp) { fprintf(stderr, "Failed to read %s\n", merge_filename); return 0; }

   int nboundaries;
   fread(&nboundaries, sizeof(int), 1, merge_fp);

   int *boundary_merge_order = new int[nboundaries];
   fread(boundary_merge_order, sizeof(int), nboundaries, merge_fp);

   // close file
   fclose(merge_fp);

   // cellulars map to themselves before merging
   int ncellulars = nd->NCellulars();
   int *cellular_mapping = new int[ncellulars];
   for (int ic = 0; ic < ncellulars; ++ic)
      cellular_mapping[ic] = ic;

   RNTime merge_time;
   merge_time.Read();
   for (int ib = 0; ib < nboundaries; ++ib) {
      // get the boundary to merge
      int boundary_to_merge = boundary_merge_order[ib];
      NeuronBoundary *boundary = nd->Boundary(boundary_to_merge);

      // defined cut off points
      float smoothing_value = smoothing_terms[boundary->DataIndex()];
      if (smoothing_value < threshold / 100.0) break;

      // get the cellular indices
      int cellular_one_data_index = boundary->SupervoxelOne()->DataIndex();
      int cellular_two_data_index = boundary->SupervoxelTwo()->DataIndex();
      rn_assertion(cellular_one_data_index < nd->NCellulars());
      rn_assertion(cellular_two_data_index < nd->NCellulars());

      // see if a merge actually occurred
      int cellular_one_label = cellular_mapping[cellular_one_data_index];
      int cellular_two_label = cellular_mapping[cellular_two_data_index];
      if (cellular_one_label == cellular_two_label) continue;

      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (cellular_mapping[ic] == cellular_two_label)
            cellular_mapping[ic] = cellular_one_label;
      }
   }

   // deflate the cellular mapping
   RNDeflateIntegerArray(cellular_mapping, ncellulars);

   // get the voxel proposals
   int nvoxels = nd->NVoxels();
   int *voxel_proposals = new int[nvoxels];
   for (int iv = 0; iv < nvoxels; ++iv) {
      int voxel_supervoxel = nd->VoxelMapping(iv);
      if (voxel_supervoxel >= nd->NCellulars()) voxel_proposals[iv] = 0;
      else voxel_proposals[iv] = cellular_mapping[voxel_supervoxel] + 1;
   }

   RNScalar *test_metrics = new RNScalar[nd->NTestMetrics()];
   RNScalar *segmentation_metrics = new RNScalar[nd->NSegmentationMetrics()];

   // run the test metric
   if (!nd->TestMetric(voxel_proposals, test_metrics, FALSE)) return 0;
   if (!nd->SegmentationMetric(cellular_mapping, segmentation_metrics)) return 0;

   printf("%0.2lf,%0.6lf\n", threshold / 100.0, test_metrics[3]);

   // free memory
   delete[] voxel_proposals;
   delete[] test_metrics;
   delete[] segmentation_metrics;

   // free memory
   delete[] cellular_mapping;
   delete[] boundary_merge_order;

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
         else if (!strcmp(*argv, "-threshold")) { argv++; argc--; threshold = atoi(*argv); }
         else if (!strcmp(*argv, "-validation")) validation = 1;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input argument
   if (!input_filename) { fprintf(stderr, "Need to supply input filename\n"); return 0; }

   // Return OK status 
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

   // read in the neuron data
   if (!ReadData(input_filename)) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // allocate memory for array
   smoothing_terms = new float[nd->NBoundaries()];

   // get output filename
   char input_filename[4096];
   sprintf(input_filename, "%s/%s_%s.binary", graphcut_directory, root_filename, extension);

   FILE *binary_fp = fopen(input_filename, "rb");
   if (!binary_fp) { fprintf(stderr, "Failed to read %s\n", input_filename); return 0; }

   // write boundary terms
   int nboundaries = nd->NBoundaries();
   fread(&nboundaries, sizeof(int), 1, binary_fp);
   fread(smoothing_terms, sizeof(float), nd->NBoundaries(), binary_fp);

   fclose(binary_fp);

   if (validation) {
      for (threshold = 0; threshold <= 100; ++threshold) {
         if (!HierarchicalClustering(root_filename)) exit(-1);
      }
   }
   else {
      if (!HierarchicalClustering(root_filename)) exit(-1);
   }

   // free memory
   delete nd;
   delete[] smoothing_terms;

   // return success
   return 0;
}
