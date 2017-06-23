// Source file to create graph cut initial positions



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_name = NULL;
static const char *unary_extension = NULL;
static R3Grid *prediction_grid = NULL;



// program variables

static NeuronData *nd;



// directory variables

static const char *algs_directory = "algs_data/graphcut";
static const char *prediction_directory = "output/ordermerge";



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
   if (!nd->ReadFile(filename, FALSE)) {
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
         else if (!strcmp(*argv, "-truth")) unary_extension = "truth";
         else if (!strcmp(*argv, "-conservative")) unary_extension = "conservative";
         else if (!strcmp(*argv, "-ordermerge")) unary_extension = "ordermerge";
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
         else { fprintf(stderr, "Invalid program argument: %s", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input argument
   if (!input_name) { fprintf(stderr, "Need to supply input filename.\n"); return 0; }
   if (!unary_extension) { fprintf(stderr, "Need to provide metric.\n"); return 0; }

   // return OK status 
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
   if (!ReadData(input_name)) exit(-1);

   // get neuron root filename
   char root_filename[4096];
   strcpy(root_filename, nd->Filename());
   char *extp = strchr(root_filename, '.');
   *extp = '\0';

   // read prediction grid
   if (!strncmp(unary_extension, "ordermerge", 10)) {
      // read voxels
      nd->ReadVoxels();

      char meta_filename[4096];
      sprintf(meta_filename, "%s/%s", prediction_directory, root_filename);

      prediction_grid = RNReadNeuronMetaRawFile(meta_filename);
      if (!prediction_grid) exit(-1);
   }

   // create all memory
   float **data_terms = new float *[nd->NCellulars()];
   for (int i1 = 0; i1 < nd->NCellulars(); ++i1) {
      data_terms[i1] = new float[nd->NCellulars()];
      for (int i2 = 0; i2 < nd->NCellulars(); ++i2) {
         data_terms[i1][i2] = FLT_MIN;
      }
   }

   // create all data terms
   if (print_verbose) printf("Creating unary graph cut terms...\n  ");

   // go through every cellular boundary and determining the smoothing term
   for (int ic1 = 0; ic1 < nd->NCellulars(); ++ic1) {
      if (print_verbose) RNProgressBar(ic1, nd->NCellulars());
      NeuronCellular *cellular_one = nd->Cellular(ic1);
      for (int ic2 = 0; ic2 < nd->NCellulars(); ++ic2) {
         NeuronCellular *cellular_two = nd->Cellular(ic2);

         // make sure cellular two is on the boundary
         if (!cellular_two->IsOnBoundary()) continue;

         // get index of cellulars
         int cellular_one_data_index = cellular_one->DataIndex();
         int cellular_two_data_index = cellular_two->DataIndex();

         // update unary term
         if (!strncmp(unary_extension, "truth", 5)) {
            int cellular_one_human_label_index = cellular_one->MajorityHumanLabel()->DataIndex();
            int cellular_two_human_label_index = cellular_two->MajorityHumanLabel()->DataIndex();

            // update the data terms
            data_terms[cellular_one_data_index][cellular_two_data_index] = 1 - (cellular_one_human_label_index == cellular_two_human_label_index);
         }
         else if (!strncmp(unary_extension, "conservative", 12)) {
            if (ic1 == ic2) {
               data_terms[cellular_one_data_index][cellular_two_data_index] = 0.0;
            }
            else {
               if (cellular_one->IsOnBoundary() && cellular_two->IsOnBoundary())
                  data_terms[cellular_one_data_index][cellular_two_data_index] = 1.0 - (cellular_one->MajorityHumanLabel() == cellular_two->MajorityHumanLabel());
               else
                  data_terms[cellular_one_data_index][cellular_two_data_index] = 1.0;
            }
         }
         else if (!strncmp(unary_extension, "ordermerge", 10)) {
            int voxel_one_index = cellular_one->CenterVoxel()->DataIndex();
            int voxel_two_index = cellular_two->CenterVoxel()->DataIndex();

            int prediction_one_value = (int)(prediction_grid->GridValue(voxel_one_index) + 0.5);
            int prediction_two_value = (int)(prediction_grid->GridValue(voxel_two_index) + 0.5);

            data_terms[cellular_one_data_index][cellular_two_data_index] = 1 - (prediction_one_value == prediction_two_value);
         }
         else { fprintf(stderr, "Unrecognized extension: %s\n", unary_extension); exit(-1); }
      }
   }
   if (print_verbose) printf("\ndone!\n");

   // get unary filename
   char unary_filename[4096];
   sprintf(unary_filename, "%s/%s_%s.unary", algs_directory, root_filename, unary_extension);

   // open unary file for writing
   FILE *unary_fp = fopen(unary_filename, "wb");
   if (!unary_fp) { fprintf(stderr, "Failed to open %s\n", unary_filename); exit(-1); }

   // write the data terms
   int unary_ncellulars = nd->NCellulars();
   fwrite(&unary_ncellulars, sizeof(int), 1, unary_fp);

   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      fwrite(data_terms[ic], sizeof(float), unary_ncellulars, unary_fp);
   }

   // close file
   fclose(unary_fp);

   // free memory
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      delete[] data_terms[ic];
   }

   delete nd;

   // return success
   return 0;
}
