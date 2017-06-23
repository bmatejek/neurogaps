// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static RNScalar maximum_affinity = 0.80;
static RNScalar minimum_affinity = 0.20;
static int instance = -1;
static int merge_instances = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *tmp_algs_directory = "algs_data/montecarlo";
static const char *algs_directory = "algs_data/features";



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
         else if (!strcmp(*argv, "-instance")) { argv++; argc--; instance = atoi(*argv); }
         else if (!strcmp(*argv, "-merge")) { merge_instances = 1; }
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

   // make sure merge instances or instance is selected
   if (!merge_instances && instance == -1) {
      fprintf(stderr, "Need to either select instance or merge instance.\n");
      return 0;
   }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

static RNScalar ScaledAffinity(RNScalar affinity) {
   if (affinity > maximum_affinity) return 1.0;
   else if (affinity < minimum_affinity) return 0.0;
   else return (affinity - minimum_affinity) / (maximum_affinity - minimum_affinity);
}



static int RandomStep(NeuronVoxel *voxel)
{
   // consider all neighbors of voxel
   RNScalar affinity_sum = 0.0;
   for (int in = 0; in < voxel->NNeighbors(); ++in) {
      NeuronVoxel *neighbor = voxel->Neighbor(in);
      if (!neighbor) continue;

      affinity_sum += ScaledAffinity(voxel->AffinityToNeighbor(neighbor));
   }

   // determine the random direction
   RNScalar random_direction = affinity_sum * RNRandomScalar();
   affinity_sum = 0.0;
   for (int in = 0; in < voxel->NNeighbors(); ++in) {
      NeuronVoxel *neighbor = voxel->Neighbor(in);
      if (!neighbor) continue;

      affinity_sum += ScaledAffinity(voxel->AffinityToNeighbor(neighbor));
      if (affinity_sum > random_direction) return in;
   }

   return -1;
}



int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // get the root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   if (!merge_instances) {
      // read voxels
      nd->ReadVoxels();

      // start statistics
      RNTime start_time;
      start_time.Read();

      unsigned int *boundary_crossings = new unsigned int[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         boundary_crossings[ib] = 0;
      }

      // run several monte carlo simulations from random voxels
      int niterations = 2621440;
      if (print_verbose) { printf("Running monte carlo simulations...\n  "); fflush(stdout); }
      for (int it = 0; it < niterations; ++it) {
         if (print_verbose) RNProgressBar(it, niterations);
         // choose a random starting voxel
         int voxel_index = (int)(RNRandomScalar() * nd->NVoxels());
         NeuronVoxel *voxel = nd->Voxel(voxel_index);

         // skip extracellulars
         if (voxel->IsExtracellular()) {
            --it;
            continue;
         }

         // run for nsteps
         int nsteps = 1000;
         for (int is = 0; is < nsteps; ++is) {
            int direction = RandomStep(voxel);

            if (direction == -1) {
               --it;
               break;
            }

            // get the neighbor direction
            NeuronVoxel *neighbor = voxel->Neighbor(direction);
            rn_assertion(neighbor != NULL);

            // see if the supervoxels are different
            if (voxel->Supervoxel() != neighbor->Supervoxel()) {
               NeuronSupervoxel *supervoxel = voxel->Supervoxel();
               NeuronSupervoxel *neighbor_supervoxel = neighbor->Supervoxel();

               // find the boundary between these supervoxels
               for (int ib = 0; ib < supervoxel->NBoundaries(); ++ib) {
                  NeuronBoundary *boundary = supervoxel->Boundary(ib);
                  // increment this boundary
                  if (boundary->OtherSupervoxel(supervoxel) == neighbor_supervoxel) {
                     int boundary_data_index = boundary->DataIndex();
                     boundary_crossings[boundary_data_index]++;
                     break;
                  }
               }
            }

            // update the current voxel
            voxel = neighbor;
         }
      }


      // get output filename
      char raw_output_filename[4096];
      sprintf(raw_output_filename, "%s/%s_monte_carlo%03d.feature", tmp_algs_directory, root_filename, instance);

      // open file
      FILE *raw_fp = fopen(raw_output_filename, "wb");
      if (!raw_fp) { fprintf(stderr, "Failed to write to %s\n", raw_output_filename); exit(-1); }

      // write the number of boundaries
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, raw_fp);

      // write to all boundaries
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         fwrite(&(boundary_crossings[ib]), sizeof(unsigned int), 1, raw_fp);
      }

      // close file
      fclose(raw_fp);

      // free memory
      delete[] boundary_crossings;

      // print statistics
      printf("\nCompleted in %0.2f seconds\n", start_time.Elapsed());
   }
   else {
      /* TODO fix hard coded value */
      int ninstances = 200;

      // keep track of boundary cross occurrences
      unsigned long *boundary_crossings = new unsigned long[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib)
         boundary_crossings[ib] = 0;

      // gather all instances
      for (int ii = 0; ii < ninstances; ++ii) {
         // get input filename
         char raw_input_filename[4096];
         sprintf(raw_input_filename, "%s/%s_monte_carlo%03d.feature", tmp_algs_directory, root_filename, ii);

         // open file
         FILE *raw_fp = fopen(raw_input_filename, "rb");
         if (!raw_fp) { fprintf(stderr, "Failed to read %s\n", raw_input_filename); exit(-1); }

         // read the number of boudaries
         int ninput_boundaries;
         fread(&ninput_boundaries, sizeof(int), 1, raw_fp);

         unsigned int *input_crossings = new unsigned int[nd->NBoundaries()];
         fread(input_crossings, sizeof(unsigned int), nd->NBoundaries(), raw_fp);

         // close file
         fclose(raw_fp);

         for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
            boundary_crossings[ib] += input_crossings[ib];
         }

         delete[] input_crossings;
      }

      // get output filenames
      char ratio_output_filename[4096];
      sprintf(ratio_output_filename, "%s/%s_monte_carlo_ratio.feature", algs_directory, root_filename);
      char raw_output_filename[4096];
      sprintf(raw_output_filename, "%s/%s_monte_carlo.feature", algs_directory, root_filename);

      // open files
      FILE *ratio_fp = fopen(ratio_output_filename, "wb");
      if (!ratio_fp) { fprintf(stderr, "Failed to write to %s\n", ratio_output_filename); exit(-1); }

      FILE *raw_fp = fopen(raw_output_filename, "wb");
      if (!raw_fp) { fprintf(stderr, "Failed to write to %s\n", raw_output_filename); exit(-1); }

      // write the number of boundaries
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, ratio_fp);
      fwrite(&nboundaries, sizeof(int), 1, raw_fp);


      // write to all boundaries
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         float ratio = boundary_crossings[ib] / (RNScalar)nd->Boundary(ib)->NAffinities();
         fwrite(&ratio, sizeof(float), 1, ratio_fp);
         unsigned int raw = boundary_crossings[ib];
         fwrite(&raw, sizeof(unsigned int), 1, raw_fp);
         // check for overflow errors
         rn_assertion(raw == boundary_crossings[ib]);
      }

      // close files
      fclose(ratio_fp);
      fclose(raw_fp);
   }

   // free up memory
   delete nd;


   // return success
   return 0;
}
