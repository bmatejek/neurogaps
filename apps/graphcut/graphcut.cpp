// Source file for the graph cut algorithm



// include files



#include "Neuron/Neuron.h"
#include "gco/GCoptimization.h"
#include <time.h>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_name = NULL;
static const char *unary_extension = NULL;
static const char *binary_extension = NULL;
static RNScalar binary_weight = 360.0;
static RNScalar unary_weight = 360.0;
static int write_meta = 0;



// program variables

static NeuronData *nd = NULL;
static float **data_terms = NULL;
static int *label_mapping_to_cellular_id = NULL;
static float **smoothing_terms = NULL;
static int *result_labels = NULL;
static const char *algs_directory = "algs_data/graphcut";
static const char *results_directory = "results/graphcut";
static const char *tmp_results_directory = "results/graphcut/tmp";
static const char *boundaries_directory = "boundaries";
static const char *output_directory = "output/graphcut";



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate neuron data\n");
      return 0;
   }

   nd->ReadFile(input_name, FALSE);

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }

   // return success
   return 1;
}



static int ReadBinaryFile(const char binary_filename[4096])
{
   FILE *binary_fp = fopen(binary_filename, "rb");
   if (!binary_fp) { fprintf(stderr, "Failed to read binary file %s\n", binary_filename); return 0; }

   // create array of smoothing terms
   int binary_nboundaries;
   fread(&binary_nboundaries, sizeof(int), 1, binary_fp);
   rn_assertion(binary_nboundaries == nd->NBoundaries());

   float *boundary_terms = new float[nd->NBoundaries()];
   if (fread(boundary_terms, sizeof(float), nd->NBoundaries(), binary_fp) != (unsigned int)nd->NBoundaries()) {
      fprintf(stderr, "Failed to read binary file %s\n", binary_filename); return 0;
   }

   smoothing_terms = new float *[nd->NCellulars()];
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      smoothing_terms[ic1] = new float[nd->NCellulars()];
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
	 smoothing_terms[ic1][ic2] = FLT_MAX;
      }
   }

   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get supervoxels
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      if (supervoxel_one->IsExtracellular()) continue;
      if (supervoxel_two->IsExtracellular()) continue;

      smoothing_terms[supervoxel_one->DataIndex()][supervoxel_two->DataIndex()] = boundary_terms[ib];
      smoothing_terms[supervoxel_two->DataIndex()][supervoxel_one->DataIndex()] = boundary_terms[ib];
   }

   // close file
   fclose(binary_fp);

   // return success
   return 1;
}



static int ReadUnaryFile(const char unary_filename[4096])
{
   // read in unary file
   FILE *unary_fp = fopen(unary_filename, "rb");
   if (!unary_fp) { fprintf(stderr, "Failed to read unary filename: %s\n", unary_filename); return 0; }

   // create array of data terms
   int unary_ncellulars;
   fread(&unary_ncellulars, sizeof(int), 1, unary_fp);
   rn_assertion(unary_ncellulars == nd->NCellulars());

   data_terms = new float *[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      data_terms[ic] = new float[nd->NCellulars()];
   }

   for (int i = 0; i < nd->NCellulars(); ++i) {
      fread(data_terms[i], sizeof(float), unary_ncellulars, unary_fp);
   }

   // close file
   fclose(unary_fp);

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Graphcut function calls functions
////////////////////////////////////////////////////////////////////////

static int RunGraphCut(char root_filename[4096])
{
   // run graph cut algorithm
   try {
      int num_labels = 0;
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         NeuronCellular *cellular = nd->Cellular(ic);
         if (cellular->IsOnBoundary()) num_labels++;
      }
      int num_sites = nd->NCellulars();

      // create mapping from label to cellular id
      label_mapping_to_cellular_id = new int[num_labels];
      int *cellular_id_mapping_to_label = new int[num_sites];
      int index = 0;
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         NeuronCellular *cellular = nd->Cellular(ic);
         if (cellular->IsOnBoundary()) {
            rn_assertion(cellular->DataIndex() == ic);
            label_mapping_to_cellular_id[index] = cellular->DataIndex();
            cellular_id_mapping_to_label[ic] = index;
            index++;
         }
         else {
            cellular_id_mapping_to_label[ic] = -1;
         }
      }
      printf("Setting unary data costs..."); fflush(stdout);
      GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_sites, num_labels);
      for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
         for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
            NeuronCellular *cellular_two = nd->Cellular(ic2);
            if (!cellular_two->IsOnBoundary()) continue;

            int label = cellular_id_mapping_to_label[ic2];
            int unary_term = (int)(unary_weight * data_terms[ic1][ic2] + 0.5);

            // set all data costs
            gc->setDataCost(ic1, label, unary_term);
         }
      }
      printf("done.\n");

      // set smoothing (binary) costs
      printf("Setting binary smoothing costs..."); fflush(stdout);
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         NeuronBoundary *boundary = nd->Boundary(ib);
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
         if (supervoxel_one->IsExtracellular()) continue;
         if (supervoxel_two->IsExtracellular()) continue;

         NeuronCellular *cellular_one = (NeuronCellular *)supervoxel_one;
         NeuronCellular *cellular_two = (NeuronCellular *)supervoxel_two;

         int index_one = cellular_one->DataIndex();
         int index_two = cellular_two->DataIndex();

         int smoothing_term = (int)(binary_weight * smoothing_terms[index_one][index_two] + 0.5);
         gc->setNeighbors(index_one, index_two, smoothing_term);
      }
      printf("done.\n");

      // set label costs
      printf("Setting label costs..."); fflush(stdout);
      for (int il1 = 0; il1 < num_labels; ++il1) {
         for (int il2 = 0; il2 < num_labels; ++il2) {
            if (il1 != il2) gc->setSmoothCost(il1, il2, 1);
            else gc->setSmoothCost(il1, il2, 0);
         }
      }
      printf("done.\n");

      // optimize
      printf("\nBefore optimization energy is %lf\n", gc->compute_energy());
      printf("  Data cost is %lf\n", gc->giveDataEnergy());
      printf("  Smoothing cost is %lf\n", gc->giveSmoothEnergy());
      printf("  Label cost is %lf\n", gc->giveLabelEnergy());
      RNTime start_time;
      start_time.Read();
      gc->expansion();
      printf("Completed in %0.2f seconds\n", start_time.Elapsed());
      printf("\nAfter optimization energy is %lf\n", gc->compute_energy());
      printf("  Data cost is %lf\n", gc->giveDataEnergy());
      printf("  Smoothing cost is %lf\n", gc->giveSmoothEnergy());
      printf("  Label cost is %lf\n", gc->giveLabelEnergy());

      // set labels
      if (print_verbose) { printf("Populating results matrix...\n  "); fflush(stdout); }
      result_labels = new int[nd->NCellulars()];
      for (int ic = 0; ic < nd->NCellulars(); ++ic) {
         if (print_verbose) RNProgressBar(ic, nd->NCellulars());
         result_labels[ic] = gc->whatLabel(ic);
      }
      if (print_verbose) printf("\ndone!\n");
   }
   catch (GCException e){
      e.Report();
      return 0;
   }

   // deflate the entries to speed up running time
   if (print_verbose) { printf("Deflating proposal array..."); fflush(stdout); }
   RNDeflateIntegerArray(result_labels, nd->NCellulars());
   if (print_verbose) printf("done!\n");

   // output results to boundary folder
   char boundary_filename[4096];
   sprintf(boundary_filename, "%s/%s_graphcut_%06d.boundary", boundaries_directory, root_filename, (int)(unary_weight + 0.5));

   // open file
   FILE *boundaries_fp = fopen(boundary_filename, "wb");
   if (!boundaries_fp) { fprintf(stderr, "Failed to write to %s\n", boundary_filename); return 0; }

   int nboundaries = nd->NBoundaries();
   fwrite(&nboundaries, sizeof(int), 1, boundaries_fp);

   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get supervoxels
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      // get data indices
      int supervoxel_one_data_index = supervoxel_one->DataIndex();
      int supervoxel_two_data_index = supervoxel_two->DataIndex();

      RNBoolean merged = FALSE;
      if (!supervoxel_one->IsExtracellular() && !supervoxel_two->IsExtracellular()) {
         // get labels for both cellulars
         int cellular_one_label = result_labels[supervoxel_one_data_index];
         int cellular_two_label = result_labels[supervoxel_two_data_index];

         merged = (cellular_one_label == cellular_two_label);
      }

      // write the result
      fwrite(&merged, sizeof(RNBoolean), 1, boundaries_fp);
   }

   // close file
   fclose(boundaries_fp);

   // return success
   return 1;
}



static int RunRandTest(const char root_filename[4096])
{
   RNTime start_time;
   start_time.Read();

   // create the voxel proposals matrix
   int *voxel_proposals = new int[nd->NVoxels()];
   if (print_verbose) { printf("Creating voxel proposals...\n  "); fflush(stdout); }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (print_verbose) RNProgressBar(iv, nd->NVoxels());
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->IsExtracellular()) voxel_proposals[iv] = 0;
      else {
         NeuronCellular *cellular = (NeuronCellular *)voxel->Supervoxel();
         int cellular_data_index = cellular->DataIndex();
         int label = result_labels[cellular_data_index];
         voxel_proposals[iv] = label + 1;
      }
   }
   if (print_verbose) printf("\ndone!\n");

   // get test metrics
   int ntest_metrics = nd->NTestMetrics();
   RNScalar *test_metrics = new RNScalar[ntest_metrics];

   if (!nd->TestMetric(voxel_proposals, test_metrics)) return 0;
   delete[] voxel_proposals;

   // get segmentation metrics
   int nsegmentation_metrics = nd->NSegmentationMetrics();
   RNScalar *segmentation_metrics = new RNScalar[nsegmentation_metrics];

   if (!nd->SegmentationMetric(result_labels, segmentation_metrics)) return 0;

   // get filename
   char rand_error_filename[4096];
   if (write_meta) sprintf(rand_error_filename, "%s/%s_%s_%s.rand", results_directory, root_filename, unary_extension, binary_extension);
   else sprintf(rand_error_filename, "%s/%s_%s_%s_%06d.rand", tmp_results_directory, root_filename, unary_extension, binary_extension, (int)(unary_weight + 0.5));

   // open file
   FILE *rand_fp = fopen(rand_error_filename, "wb");
   if (!rand_fp) { fprintf(stderr, "Failed to write to %s\n", rand_error_filename); return 0; }

   // write the metric results
   fwrite(&ntest_metrics, sizeof(int), 1, rand_fp);
   fwrite(&nsegmentation_metrics, sizeof(int), 1, rand_fp);
   fwrite(&unary_weight, sizeof(RNScalar), 1, rand_fp);
   fwrite(&binary_weight, sizeof(RNScalar), 1, rand_fp);
   fwrite(test_metrics, sizeof(RNScalar), ntest_metrics, rand_fp);
   fwrite(segmentation_metrics, sizeof(RNScalar), nsegmentation_metrics, rand_fp);

   // close file
   fclose(rand_fp);

   // free memory
   delete[] test_metrics;
   delete[] segmentation_metrics;

   // return success
   return 1;
}



static int WriteMetaRawFile(const char root_filename[4096])
{
   // get the meta/raw output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s", output_directory, root_filename);

   // output NMRMeta file
   RNMeta graphcut_meta = RNMeta("Int32", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 1);

   R3Grid *graphcut_raw = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);
            if (voxel->IsExtracellular()) graphcut_raw->SetGridValue(ix, iy, iz, 0);
            // cellulars have predictions
            else {
               NeuronCellular *cellular = (NeuronCellular *)voxel->Supervoxel();
               int label = result_labels[cellular->DataIndex()];
               graphcut_raw->SetGridValue(ix, iy, iz, label + 1);
            }
         }
      }
   }
   RNWriteNeuronMetaRawFile(output_filename, graphcut_meta, graphcut_raw);

   // free memory
   delete graphcut_raw;

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
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-binary_truth")) binary_extension = "truth";
         else if (!strcmp(*argv, "-unary_truth")) unary_extension = "truth";
         else if (!strcmp(*argv, "-binary_random_forest")) binary_extension = "random_forest_smoothing";
	 else if (!strcmp(*argv, "-binary_random_forest_global")) binary_extension = "random_forest_smoothing_global";
         else if (!strcmp(*argv, "-unary_random_forest")) unary_extension = "random_forest";
         else if (!strcmp(*argv, "-unary_conservative")) unary_extension = "conservative";
         else if (!strcmp(*argv, "-unary_ordermerge")) unary_extension = "ordermerge";
         else if (!strcmp(*argv, "-binary_weight")) { argv++; argc--; binary_weight = atof(*argv); }
         else if (!strcmp(*argv, "-unary_weight")) { argv++; argc--; unary_weight = atof(*argv); }
         else if (!strcmp(*argv, "-write_meta")) write_meta = 1;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // check filenames
   if (!input_name) { fprintf(stderr, "Must include input neuron data filename\n"); return 0; }
   if (!unary_extension) { fprintf(stderr, "Must include input unary extension\n"); return 0; }
   if (!binary_extension) { fprintf(stderr, "Must include input binary extension\n"); return 0; }

   // return OK status
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in neuron file
   if (!ReadData()) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strchr(root_filename, '.');
   *extp = '\0';

   // create binary and unary input files
   char binary_filename[4096];
   sprintf(binary_filename, "%s/%s_%s.binary", algs_directory, root_filename, binary_extension);
   char unary_filename[4096];
   sprintf(unary_filename, "%s/%s_%s.unary", algs_directory, root_filename, unary_extension);

   if (!ReadBinaryFile(binary_filename)) exit(-1);
   if (!ReadUnaryFile(unary_filename)) exit(-1);

   // run graphcut algorithm
   if (!RunGraphCut(root_filename)) exit(-1);

   // read in voxels
   nd->ReadVoxels();

   if (write_meta && !WriteMetaRawFile(root_filename)) exit(-1);

   // free memory
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      delete[] data_terms[ic];
      delete[] smoothing_terms[ic];
   }
   delete[] data_terms;
   delete[] smoothing_terms;

   // run testing algorithm
   if (!RunRandTest(root_filename)) exit(-1);

   // free memory
   delete[] result_labels;
   delete nd;

   // return success
   return 0;
}
