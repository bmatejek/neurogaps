// Source file for the subsection algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
// proportion along all dimensions
static int xsplit = 0;
static int ysplit = 0;
static int zsplit = 0;



// global variables

static NeuronData *nd = NULL;
static char *input_filename = NULL;
static char root_filename[4096];



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
         else if (!strcmp(*argv, "-debug")) { print_debug = 1; print_verbose = 1; }
         // scale changes
         else if (!strcmp(*argv, "-x")) { argv++; argc--; xsplit = atoi(*argv); }
         else if (!strcmp(*argv, "-y")) { argv++; argc--; ysplit = atoi(*argv); }
         else if (!strcmp(*argv, "-z")) { argv++; argc--; zsplit = atoi(*argv); }
         // number of files to generate
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!input_filename) { fprintf(stderr, "Need to supply input filename.\n");  return 0; }

   // make sure each split integer is define
   if (!xsplit || !ysplit || !zsplit) { fprintf(stderr, "Must define the number of splits for each dimensions\n"); return 0; }

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

   // read in the nueron data
   if (!ReadData(input_filename)) exit(-1);

   // create input filename root
   strncpy(root_filename, nd->Filename(), 4096);
   char *endp = strrchr(root_filename, '.');
   if (endp) *endp = '\0';

   // make sure the dimension is split evenly
   if (nd->XResolution() % xsplit != 0) {
      fprintf(stderr, "Must split grid evenly (x dimension violated)\n");
      return 0;
   }
   if (nd->YResolution() % ysplit != 0) {
      fprintf(stderr, "Must split grid evenly (y dimension violated)\n");
      return 0;
   }
   if (nd->ZResolution() % zsplit != 0) {
      fprintf(stderr, "Must split grid evenly (z dimension violated)\n");
      return 0;
   }

   // size of each dimension
   int nx = nd->XResolution() / xsplit;
   int ny = nd->YResolution() / ysplit;
   int nz = nd->ZResolution() / zsplit;

   int nfiles = 1;
   // go through x dimension
   for (int ix = 0; ix < nd->XResolution(); ix += nx) {
      // go through y dimension
      for (int iy = 0; iy < nd->YResolution(); iy += ny) {
         // go through z dimension
         for (int iz = 0; iz < nd->ZResolution(); iz += nz, ++nfiles) {
            // start statistics
            RNTime subsection_time;
            subsection_time.Read();

            // get dimensions
            int xmin = ix;
            int ymin = iy;
            int zmin = iz;
            int xmax = ix + nx - 1;
            int ymax = iy + ny - 1;
            int zmax = iz + nz - 1;

            // create filename
            char output_root_filename[4096];
            sprintf(output_root_filename, "%s%d", root_filename, nfiles);

            // file widths
            int xwidth = xmax - xmin + 1;
            int ywidth = ymax - ymin + 1;
            int zwidth = zmax - zmin + 1;

            // create new R3Grids
            R3Grid *affinities_grids[3];
            R3Grid *human_labels_grid;
            R3Grid *image_grid;
            R3Grid *machine_labels_grid;

            // allocate memory
            for (int dim = 0; dim <= 2; ++dim) {
               affinities_grids[dim] = new R3Grid(xwidth, ywidth, zwidth);
            }
            human_labels_grid = new R3Grid(xwidth, ywidth, zwidth);
            image_grid = new R3Grid(xwidth, ywidth, zwidth);
            machine_labels_grid = new R3Grid(xwidth, ywidth, zwidth);

            int subi = 0;
            for (int oldi = xmin; oldi <= xmax; ++oldi, ++subi) {
               int subj = 0;
               for (int oldj = ymin; oldj <= ymax; ++oldj, ++subj) {
                  int subk = 0;
                  for (int oldk = zmin; oldk <= zmax; ++oldk, ++subk) {
                     NeuronVoxel *voxel = nd->Voxel(oldi, oldj, oldk);

                     // update grids
                     for (int dim = 0; dim <= 2; ++dim)
                        affinities_grids[dim]->SetGridValue(subi, subj, subk, voxel->Affinity(dim));
                     if (voxel->HumanLabel())
                        human_labels_grid->SetGridValue(subi, subj, subk, voxel->HumanLabel()->DataIndex() + 1);
                     else
                        human_labels_grid->SetGridValue(subi, subj, subk, 0);
                     image_grid->SetGridValue(subi, subj, subk, voxel->Image());
                     if (voxel->IsCellular())
                        machine_labels_grid->SetGridValue(subi, subj, subk, voxel->Supervoxel()->DataIndex() + 1);
                     else
                        machine_labels_grid->SetGridValue(subi, subj, subk, 0);
                  }
               }
            }

            // write meta filenames
            RNMeta affinity_meta = RNMeta("Float32", xwidth, ywidth, zwidth, 3);
            RNMeta human_label_meta = RNMeta("Uint16", xwidth, ywidth, zwidth, 1);
            RNMeta image_meta = RNMeta("Float32", xwidth, ywidth, zwidth, 1);
            RNMeta machine_label_meta = RNMeta("Int32", xwidth, ywidth, zwidth, 1);

            char affinity_filename[4096];
            sprintf(affinity_filename, "%s/affinities/%s_affinities", nd->FilePath(), output_root_filename);
            char human_label_filename[4096];
            sprintf(human_label_filename, "%s/human_labels/%s_human_labels", nd->FilePath(), output_root_filename);
            char image_filename[4096];
            sprintf(image_filename, "%s/images/%s_image", nd->FilePath(), output_root_filename);
            char machine_label_filename[4096];
            sprintf(machine_label_filename, "%s/machine_labels/%s_machine_labels", nd->FilePath(), output_root_filename);

            // write neuron meta/raw files
            if (!RNWriteNeuronMetaRawFile(affinity_filename, affinity_meta, affinities_grids)) { fprintf(stderr, "Failed to write %s\n", affinity_filename); exit(-1); }
            if (!RNWriteNeuronMetaRawFile(human_label_filename, human_label_meta, human_labels_grid)) { fprintf(stderr, "Failed to write %s\n", human_label_filename); exit(-1); }
            if (!RNWriteNeuronMetaRawFile(image_filename, image_meta, image_grid)) { fprintf(stderr, "Failed to write %s\n", image_filename); exit(-1); }
            if (!RNWriteNeuronMetaRawFile(machine_label_filename, machine_label_meta, machine_labels_grid)) { fprintf(stderr, "Failed to write %s\n", machine_label_filename); exit(-1); }

            char output_filename[4096];
            sprintf(output_filename, "%s/%s.txt", nd->FilePath(), output_root_filename);

            // open file
            FILE *fp = fopen(output_filename, "w");
            if (!fp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

            char *affinity_root = strrchr(affinity_filename, '/');
            affinity_root++;
            char *human_label_root = strrchr(human_label_filename, '/');
            human_label_root++;
            char *image_root = strrchr(image_filename, '/');
            image_root++;
            char *machine_label_root = strrchr(machine_label_filename, '/');
            machine_label_root++;

            fprintf(fp, "affinities/%s.meta\n", affinity_root);
            fprintf(fp, "human_labels/%s.meta\n", human_label_root);
            fprintf(fp, "images/%s.meta\n", image_root);
            fprintf(fp, "machine_labels/%s.meta\n", machine_label_root);
            fprintf(fp, "Scale(%f,%f,%f)\n", nd->XScale(), nd->YScale(), nd->ZScale());

            // close file
            fclose(fp);

            // free memory
            for (int dim = 0; dim <= 2; ++dim)
               delete affinities_grids[dim];
            delete human_labels_grid;
            delete image_grid;
            delete machine_labels_grid;

            // print statistics
            if (print_verbose) printf("Created subsection (%d, %d, %d) - (%d, %d, %d) in %.2f seconds\n", xmin, ymin, zmin, xmax, ymax, zmax, subsection_time.Elapsed());
         }
      }
   }

   // future steps
   printf("Need to deflateraw for machine labels and human labels.\n");

   // free memory
   delete nd;

   // return success
   return 0;
}