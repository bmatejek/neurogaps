// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd = NULL;
static std::vector<const char *> training_filenames = std::vector<const char *>();
static std::vector<const char *> testing_filenames = std::vector<const char *>();
static RNPyImage *image = NULL;



// best fit variables

static RNScalar alpha = RN_UNKNOWN;
static RNScalar beta = RN_UNKNOWN;
static RNScalar rsquared = RN_UNKNOWN;
static RNScalar *probabilities = NULL;
static const int nbins = 100;



// directory structure

static const char *algs_directory = "algs_data/affinity";
static const char *results_directory = "results/affinity";
static const char *visuals_directory = "visuals/affinity";
static const char *filters_directory = "filters";



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

static int ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
         else if (!strcmp(*argv, "-train")) {
            argv++; argc--;
            while ((argc > 0) && ((*argv)[0] != '-')) {
               training_filenames.push_back(*argv);
               argv++; argc--;
            }
            argv--; argc++;
         }
         else if (!strcmp(*argv, "-test")) {
            argv++; argc--;
            while ((argc > 0) && ((*argv)[0] != '-')) {
               testing_filenames.push_back(*argv);
               argv++; argc--;
            }
            argv--; argc++;
         }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0;
      }
      argv++; argc--;
   }

   // make sure there are training and testing files
   if (!training_filenames.size() || !testing_filenames.size()) {
      fprintf(stderr, "Need to supply training and testing files.\n");
      return 0;
   }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int UpdateBinProbabilities(unsigned int *nsame, unsigned int *ndifferent)
{
   // just checking
   rn_assertion(nd != NULL);
   rn_assertion(nsame != NULL);
   rn_assertion(ndifferent != NULL);

   // go through all affinities
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      // go through all of the neighbors
      for (int in = 0; in < voxel->NNeighbors(); in += 2) {
         NeuronVoxel *neighbor = voxel->Neighbor(in);
         if (!neighbor) continue;

         // get ground truth for this affinity 
         RNBoolean label = (voxel->HumanLabel() == neighbor->HumanLabel());
         if (label && !voxel->HumanLabel()) label = FALSE;

         // get the affinity between these voxels
         RNScalar affinity = voxel->AffinityToNeighbor(neighbor);

         // to which bin does this belong
         int bin = (int)(nbins * affinity + 0.5);

         // increment counters
         if (label) nsame[bin]++;
         else ndifferent[bin]++;
      }
   }

   // return success
   return 1;
}



int GenerateBestFitLine(void)
{
   // create bins for histogram
   unsigned int *nsame = new unsigned int[nbins];
   unsigned int *ndifferent = new unsigned int[nbins];
   for (int ib = 0; ib < nbins; ++ib) {
      nsame[ib] = 0;
      ndifferent[ib] = 0;
   }

   // read in all training files
   for (unsigned int it = 0; it < training_filenames.size(); ++it) {
      // read the data
      if (!ReadData(training_filenames[it])) exit(-1);

      // update the probability bins for this neuron file
      if (!UpdateBinProbabilities(nsame, ndifferent));

      // free up memory
      delete nd;
   }

   // create an output py plot
   RNPyPointPlot *point_plot = new RNPyPointPlot();
   point_plot->SetTitle("Affinity to Probability Map");
   point_plot->SetXLabel("Affinity");
   point_plot->SetYLabel("Probability");

   // set axis constraints
   point_plot->SetXAxisMin(0);
   point_plot->SetXAxisMax(1);

   point_plot->SetYAxisMin(0);
   point_plot->SetYAxisMax(1);

   // fit function to training set
   probabilities = new RNScalar[nbins];
   RNScalar *x = new RNScalar[nbins];
   RNScalar *y = new RNScalar[nbins];

   for (int ib = nbins - 1; ib >= 0; --ib) {
      // only consider places where training exists
      if (nsame[ib] + ndifferent[ib] != 0) {
         probabilities[ib] = nsame[ib] / (RNScalar)(nsame[ib] + ndifferent[ib]);
         // for monotonic function
         if (ib != nbins - 1 && probabilities[ib] > probabilities[ib + 1])
            probabilities[ib] = probabilities[ib + 1];
      }
      // forces monotocity and completeness
      else {
         if (ib == nbins - 1) probabilities[ib] = 1.0;
         else probabilities[ib] = probabilities[ib + 1];
      }

      // set x and y values and add to point plot
      x[ib] = ib / (RNScalar)nbins;
      y[ib] = probabilities[ib];
      point_plot->InsertPoint(R2Point(x[ib], y[ib]));
   }

   // add to the overall image
   image->InsertPyPlot(point_plot);

   // run the best fit approximation
   RNBestFitLine(x, y, nbins, alpha, beta, rsquared);

   // save the best fit line
   char output_filename[4096];
   sprintf(output_filename, "%s/affinity_map.function", algs_directory);

   // open file
   FILE *fp = fopen(output_filename, "w");
   if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // write file information
   fprintf(fp, "LINEAR_FUNCTION\n");
   fprintf(fp, "NINDEPENDENT_VARIABLES 1\n");
   fprintf(fp, "R-SQUARED: %lf\n", rsquared);
   fprintf(fp, "%lf\n", alpha);
   fprintf(fp, "%lf\n", beta);
   fprintf(fp, "\nPROBABILITIES\n");
   for (int ib = 0; ib < nbins; ++ib) {
      fprintf(fp, "%4d %lf\n", ib, probabilities[ib]);
   }

   // close file
   fclose(fp);

   // free up memory
   delete[] x;
   delete[] y;
   delete[] nsame;
   delete[] ndifferent;

   // print the results
   printf("Created best fit probability functions: \n");
   printf("  alpha:    %0.6lf\n", alpha);
   printf("  beta:     %0.6lf\n", beta);
   printf("  rsquared: %0.6lf\n", rsquared);

   // return success
   return 1;
}



int TestDataSet(void)
{
   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create an output py plot
   RNPyPointPlot *point_plot = new RNPyPointPlot();

   // create title
   char title[4096];
   sprintf(title, "Affinity to Probability Map for %s", root_filename);

   // set the point plots
   point_plot->SetTitle(title);
   point_plot->SetXLabel("Affinity");
   point_plot->SetYLabel("Proportion Same Neuron");

   // set axis constraints
   point_plot->SetXAxisMin(0);
   point_plot->SetXAxisMax(1);

   point_plot->SetYAxisMin(0);
   point_plot->SetYAxisMax(1);

   // create counters for same and different labels
   unsigned int *nsame = new unsigned int[nbins];
   unsigned int *ndifferent = new unsigned int[nbins];
   for (int ib = 0; ib < nbins; ++ib) {
      nsame[ib] = 0;
      ndifferent[ib] = 0;
   }

   // update the bin probabilities
   if (!UpdateBinProbabilities(nsame, ndifferent)) return 0;

   // see the average bin difference
   printf("Calculating differences in probabilities for %s\n", nd->Filename());
   for (int ib = nbins - 1; ib >= 0; --ib) {
      // ignore bins that do not have any data
      if (nsame[ib] + ndifferent[ib] == 0) continue;

      // find out the probability of this bin
      RNScalar probability = (RNScalar)nsame[ib] / (RNScalar)(nsame[ib] + ndifferent[ib]);

      // get the result from linear regression
      RNScalar x = ib / (RNScalar)nbins;
      RNScalar linear_prediction = alpha + beta * x;

      // print difference between prediction and probability
      printf("%4.2lf %6.2lf%% %6.2lf%%\n", ib / (RNScalar)nbins, 100 * (linear_prediction - probability), 100 * (probabilities[ib] - probability));

      // insert point to plot
      point_plot->InsertPoint(R2Point(x, probability));
   }

   // add to overall image
   image->InsertPyPlot(point_plot);

   // output probability filter
   RNMeta meta_filter = RNMeta("Float32", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 3);
   R3Grid *raw_filters[3];
   for (int dim = 0; dim <= 2; ++dim)
      raw_filters[dim] = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);

      // get affinities in each direction
      RNScalar xaffinity = voxel->XAffinity();
      int xbin = (int)(nbins * xaffinity + 0.5);
      RNScalar xprob = probabilities[xbin];
      RNScalar yaffinity = voxel->YAffinity();
      int ybin = (int)(nbins * yaffinity + 0.5);
      RNScalar yprob = probabilities[ybin];
      RNScalar zaffinity = voxel->ZAffinity();
      int zbin = (int)(nbins * zaffinity + 0.5);
      RNScalar zprob = probabilities[zbin];

      // set the probability filter
      raw_filters[RN_X]->SetGridValue(iv, xprob);
      raw_filters[RN_Y]->SetGridValue(iv, yprob);
      raw_filters[RN_Z]->SetGridValue(iv, zprob);
   }

   // get the output filename
   char meta_filename[4096];
   sprintf(meta_filename, "%s/%s_probabilities", filters_directory, root_filename);

   // save meta/raw file
   if (!RNWriteNeuronMetaRawFile(meta_filename, meta_filter, raw_filters)) { fprintf(stderr, "Failed to write meta/raw %s\n", meta_filename); return 0; }

   // free memory
   delete[] nsame;
   delete[] ndifferent;
   for (int dim = 0; dim <= 2; ++dim)
      delete raw_filters[dim];

   // return success
   return 1;
}



int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // create file for image
   image = new RNPyImage();
   image->SetOutputDirectory(visuals_directory);

   // generate the best fit line
   if (!GenerateBestFitLine()) exit(-1);

   // run through all tests
   for (unsigned int it = 0; it < testing_filenames.size(); ++it) {
      // read the data
      if (!ReadData(testing_filenames[it])) exit(-1);

      // test on this data set
      if (!TestDataSet()) exit(-1);

      // free up memory
      delete nd;
   }

   // free memory
   delete[] probabilities;

   // create output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/affinity_map.pyimage", results_directory);

   // write the image file
   image->WriteImageFile(output_filename);

   // return success
   return 0;
}
