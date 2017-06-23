// Source file for the neuron statistics algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static RNScalar input_delta = -10;



// global variables

static NeuronData *nd = NULL;
static char *input_filename = NULL;



// directory structure

static const char *results_directory = "results/affinity";
static const char *visuals_directory = "visuals/affinity";



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

static int WatershedStatistics(void)
{
   unsigned int nvessels = 0;
   unsigned int nbackground = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel()) nvessels++;
      else nbackground++;
   }

   unsigned int nincorrectly_labeled_vessel = 0;
   unsigned int ncorrectly_labeled_vessel = 0;
   unsigned int nincorrectly_labeled_background = 0;
   unsigned int ncorrectly_labeled_background = 0;
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      NeuronCellular *cellular = nd->Cellular(ic);
      for (int iv = 0; iv < cellular->NVoxels(); ++iv) {
         NeuronVoxel *voxel = cellular->LocalVoxel(iv);
         if (!voxel->HumanLabel()) nincorrectly_labeled_vessel++;
         else ncorrectly_labeled_vessel++;
      }
   }
   for (int ie = 0; ie < nd->NExtracellulars(); ++ie) {
      NeuronExtracellular *extracellular = nd->Extracellular(ie);
      for (int iv = 0; iv < extracellular->NVoxels(); ++iv) {
         NeuronVoxel *voxel = extracellular->LocalVoxel(iv);
         if (voxel->HumanLabel()) nincorrectly_labeled_background++;
         else ncorrectly_labeled_background++;
      }
   }

   // print total number of incorrectly labels
   printf("Total voxels in vessel: %u\n", nvessels);
   printf("Total voxels in background: %u\n", nbackground);
   printf("Proportion of voxels in vessels: %lf\n\n", nvessels / (RNScalar)nd->NVoxels());

   printf("Correctly labeled vessel: %u\n", ncorrectly_labeled_vessel);
   printf("Incorrectly labeled vessel: %u\n", nincorrectly_labeled_vessel);
   printf("Proportion labeled vessels that are vessels: %lf\n\n", (ncorrectly_labeled_vessel) / (RNScalar)(ncorrectly_labeled_vessel + nincorrectly_labeled_vessel));

   printf("Correctly labeled background: %u\n", ncorrectly_labeled_background);
   printf("Incorrectly labeled background: %u\n", nincorrectly_labeled_background);
   printf("Proportion labeled background that are background: %lf\n\n", (ncorrectly_labeled_background) / (RNScalar)(ncorrectly_labeled_background + nincorrectly_labeled_background));

   printf("Correctly labeled: %u\n", ncorrectly_labeled_background + ncorrectly_labeled_vessel);
   printf("Proportion correctly labeled: %lf\n", (nd->NVoxels() - nincorrectly_labeled_background - nincorrectly_labeled_vessel) / (RNScalar)nd->NVoxels());

   // return success
   return 1;
}



struct AffinityStruct {
   RNScalar affinity;
   RNBoolean both_vessel;
};



int AffinityCompare(AffinityStruct a, AffinityStruct b)
{
   return a.affinity > b.affinity;
}



static int AffinityStatistics(const char root_filename[4096])
{
   // keep track of how many affinities occur between blood vessels
   unsigned long long nbetween_vessels = 0;

   // get all of the affinities
   std::vector<AffinityStruct> affinities = std::vector<AffinityStruct>();
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      for (int in = 0; in < voxel->NNeighbors(); in += 2) {
         NeuronVoxel *neighbor = voxel->Neighbor(in);
         if (!neighbor) continue;
         AffinityStruct affinity;
         affinity.affinity = voxel->AffinityToNeighbor(neighbor);

         // are they both blood vessels
         if (voxel->HumanLabel() && neighbor->HumanLabel()) {
            affinity.both_vessel = TRUE;
            nbetween_vessels++;
         }
         else affinity.both_vessel = FALSE;

         affinities.push_back(affinity);
      }
   }

   // sort the affinities
   sort(affinities.begin(), affinities.end(), &AffinityCompare);

   // get filename
   char pandr_filename[4096];
   sprintf(pandr_filename, "%s/%s_affinities_precision_recall.pyimage", results_directory, root_filename);

   // create the plots
   RNPyPointPlot *point_plot = new RNPyPointPlot();
   point_plot->SetTitle("Affinity Precision and Recall");
   point_plot->SetXLabel("Recall");
   point_plot->SetYLabel("Precision");

   unsigned long long ncorrect = 0;
   for (unsigned int ia = 0; ia < affinities.size(); ++ia) {
      AffinityStruct affinity = affinities[ia];
      if (affinity.both_vessel) ncorrect++;

      RNScalar recall = ncorrect / (RNScalar)nbetween_vessels;
      RNScalar precision = ncorrect / (RNScalar)(ia + 1);

      point_plot->InsertPoint(R2Point(recall, precision));
   }

   // create an image
   RNPyImage *image = new RNPyImage();
   image->SetOutputDirectory(visuals_directory);
   image->InsertPyPlot(point_plot);

   // write the image file
   image->WriteImageFile(pandr_filename);

   // return success
   return 1;
}



static int RecedeAffinity(const char root_filename[4096], RNScalar delta)
{
   // remove voxels from blood vessels if the highest affinity is less than delta
   RNBoolean *blood_vessel = new RNBoolean[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      blood_vessel[iv] = FALSE;
   }

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      for (int dim = 0; dim <= 2; ++dim) {
         if (voxel->Affinity(dim) > delta) blood_vessel[iv] = TRUE;
      }
   }

   unsigned int nvessels = 0;
   unsigned int nbackground = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->HumanLabel()) nvessels++;
      else nbackground++;
   }

   unsigned int nincorrectly_labeled_vessel = 0;
   unsigned int ncorrectly_labeled_vessel = 0;
   unsigned int nincorrectly_labeled_background = 0;
   unsigned int ncorrectly_labeled_background = 0;
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (blood_vessel[iv]) {
         if (voxel->HumanLabel()) ncorrectly_labeled_vessel++;
         else nincorrectly_labeled_vessel++;
      }
      else {
         if (voxel->HumanLabel()) nincorrectly_labeled_background++;
         else ncorrectly_labeled_background++;
      }
   }

   // print total number of incorrectly labels
   printf("Results for delta: %0.4f\n", delta);
   printf("  Total voxels in vessel: %u\n", nvessels);
   printf("  Total voxels in background: %u\n", nbackground);
   printf("  Proportion of voxels in vessels: %lf\n", nvessels / (RNScalar)nd->NVoxels());

   printf("  Correctly labeled vessel: %u\n", ncorrectly_labeled_vessel);
   printf("  Incorrectly labeled vessel: %u\n", nincorrectly_labeled_vessel);
   printf("  Proportion labeled vessels that are vessels: %lf\n", (ncorrectly_labeled_vessel) / (RNScalar)(ncorrectly_labeled_vessel + nincorrectly_labeled_vessel));

   printf("  Correctly labeled background: %u\n", ncorrectly_labeled_background);
   printf("  Incorrectly labeled background: %u\n", nincorrectly_labeled_background);
   printf("  Proportion labeled background that are background: %lf\n", (ncorrectly_labeled_background) / (RNScalar)(ncorrectly_labeled_background + nincorrectly_labeled_background));

   printf("  Correctly labeled: %u\n", ncorrectly_labeled_background + ncorrectly_labeled_vessel);
   printf("  Proportion correctly labeled: %lf\n", (nd->NVoxels() - nincorrectly_labeled_background - nincorrectly_labeled_vessel) / (RNScalar)nd->NVoxels());

   if (input_delta > 0.0) {
      // create output meta file
      RNMeta meta = RNMeta("Uint16", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 1);
      R3Grid *raw = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());

      // write to raw 
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (blood_vessel[iv]) raw->SetGridValue(iv, 1);
         else raw->SetGridValue(iv, 0);
      }

      // filename
      char filename[4096];
      sprintf(filename, "%s/%s_receded", results_directory, root_filename);

      RNWriteNeuronMetaRawFile(filename, meta, raw);

      // free memory
      delete raw;
   }

   // free memory
   delete[] blood_vessel;

   // return sucess
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
         else if (!strcmp(*argv, "-delta")) { argv++; argc--; input_delta = atof(*argv); }
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

   if (!AffinityStatistics(root_filename)) exit(-1);
   if (!WatershedStatistics()) exit(-1);
   /*if (input_delta < 0.0) {
      RNScalar delta = 0.90;
      for (int id = 0; id < 100; ++id) {
         if (!RecedeAffinity(root_filename, delta)) exit(-1);
         delta += 0.1 / 100;
      }
   }
   else if (!RecedeAffinity(root_filename, input_delta)) exit(-1);*/



   // free up memory
   delete nd;

   // return success
   return 0;
}