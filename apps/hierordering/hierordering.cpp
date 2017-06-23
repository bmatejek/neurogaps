// Source file for the simple file conversion algorithm



// include files

#include "Neuron/Neuron.h"
#include <algorithm>
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static const char *input_filename = NULL;
static const char *extension = NULL;



// global variables

NeuronData *nd = NULL;


// directory structure

static const char *hierarchical_directory = "algs_data/hierarchical";
static const char *results_directory = "results/hierarchical";
static const char *visuals_directory = "visuals/hierarchical";



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static NeuronData *ReadData(const char *filename)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   NeuronData *nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate memory for neuron data\n");
      return NULL;
   }

   // read in the file
   if (!nd->ReadFile(filename)) {
      fprintf(stderr, "Failed to read file\n");
      return NULL;
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

   // return the neuron structure
   return nd;
}



struct BoundaryStruct {
   BoundaryStruct(NeuronBoundary *boundary, RNScalar prediction, RNBoolean ground_truth) {
      this->boundary = boundary;
      this->prediction = prediction;
      this->ground_truth = ground_truth;
   }

   NeuronBoundary *boundary;
   RNScalar prediction;
   RNBoolean ground_truth;
};



int BoundaryCompare(BoundaryStruct one, BoundaryStruct two) {
   return one.prediction > two.prediction;
}



static int ConstructMergeOrder(void)
{
   // to add to the title for the plots
   char title_qualifier[4096];
   if (!strcmp(extension, "boundary_max")) strncpy(title_qualifier, "Boundary Max", 4096);
   else if (!strcmp(extension, "boundary_mean")) strncpy(title_qualifier, "Boundary Mean", 4096);
   else rn_assertion(FALSE);
   
   // go through all boundaries 
   std::vector<BoundaryStruct> merge_order = std::vector<BoundaryStruct>();
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get these supervoxels
      NeuronSupervoxel *supervoxel_one = (NeuronSupervoxel *)boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = (NeuronSupervoxel *)boundary->SupervoxelTwo();

      // skip boundaries with extracellulars
      if (supervoxel_one->IsExtracellular() || supervoxel_two->IsExtracellular()) continue;

      RNBoolean ground_truth = FALSE;
      if (supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel()) ground_truth = TRUE;

      if (!strcmp(extension, "boundary_max")) {
         BoundaryStruct boundary_pair = BoundaryStruct(boundary, boundary->Maximum(), ground_truth);
         merge_order.push_back(boundary_pair);
      }
      else if (!strcmp(extension, "boundary_mean")) {
         BoundaryStruct boundary_pair = BoundaryStruct(boundary, boundary->Mean(), ground_truth);
         merge_order.push_back(boundary_pair);
      }
      else { fprintf(stderr, "Unrecognized extension: %s\n", extension); return 0; }
   }

   // sort the result truth pairs
   std::sort(merge_order.begin(), merge_order.end(), &BoundaryCompare);

   // get the root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // save the hierachical merge order for smoothing
   char merge_filename[4096];
   sprintf(merge_filename, "%s/%s_%s.hier", hierarchical_directory, root_filename, extension);

   // open file
   FILE *merge_fp = fopen(merge_filename, "wb");
   if (!merge_fp) { fprintf(stderr, "Failed to write to %s\n", merge_filename); return 0; }

   int nboundaries = merge_order.size();
   fwrite(&nboundaries, sizeof(int), 1, merge_fp);
   for (int ib = 0; ib < nboundaries; ++ib) {
      NeuronBoundary *boundary = merge_order[ib].boundary;
      int boundary_data_index = boundary->DataIndex();
      fwrite(&boundary_data_index, sizeof(int), 1, merge_fp);
   }

   // close file
   fclose(merge_fp);

   // create image
   RNPyImage *proportion_image = new RNPyImage();
   proportion_image->SetOutputDirectory(visuals_directory);

   // create plots
   RNPyPointPlot *pdf_distribution = new RNPyPointPlot();
   RNPyPointPlot *cdf_distribution = new RNPyPointPlot();

   // create titles
   char pdf_title[4096];
   sprintf(pdf_title, "%s %s Probability Distribution Function", root_filename, title_qualifier);
   char cdf_title[4096];
   sprintf(cdf_title, "%s %s Cumulative Distribution Function", root_filename, title_qualifier);
   pdf_distribution->SetTitle(pdf_title);
   cdf_distribution->SetTitle(cdf_title);

   // set x and y labels
   pdf_distribution->SetXLabel("1.0 - Boundary Score");
   cdf_distribution->SetYLabel("1.0 - Boundary Score");
   pdf_distribution->SetXLabel("Proportion Correct");
   cdf_distribution->SetYLabel("Proportion Correct");

   // create the distributions
   int nbins = 1000;
   unsigned int *ncorrect = new unsigned int[nbins];
   unsigned int *nincorrect = new unsigned int[nbins];
   for (int ib = 0; ib < nbins; ++ib) {

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
         else if (!strcmp(*argv, "-boundary_max")) extension = "boundary_max";
         else if (!strcmp(*argv, "-boundary_mean")) extension = "boundary_mean";
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_filename) input_filename = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure an input file is supplied
   if (!input_filename) { fprintf(stderr, "Need to supply input filename\n"); return 0; }
   if (!extension) { fprintf(stderr, "Need to supply merge order extension (boundary_max, boundary_mean)\n"); return 0; }

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

   // read neuron data
   nd = ReadData(input_filename);
   if (!nd) exit(-1);

   // read in neuron files
   if (!ConstructMergeOrder()) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}