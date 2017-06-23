// Source file for the simple file conversion algorithm



// include files 

#include "Neuron/Neuron.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;



// directory structure

static const char *features_directory = "algs_data/boundarydata";
static const char *filter_directory = "filters";



// useful string constant

static const int NDERIVATIVES = 2;
static const char *derivative_names[NDERIVATIVES] = {
   "first_derivative", "second_derivative"
};



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
   if (!nd->ReadFile(filename, TRUE, TRUE)) {
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
// Feature creation functions
////////////////////////////////////////////////////////////////////////

static int
CreateBoundaryFeatures(char root_filename[4096])
{
   // arrays to save all data
   RNScalar *skews = new RNScalar[nd->NBoundaries()];
   RNScalar *kurtosises = new RNScalar[nd->NBoundaries()];

   // number of bins for histogram
   const int nbins = 25;
   RNScalar *histograms[nbins];
   for (int in = 0; in < nbins; ++in) {
      histograms[in] = new RNScalar[nd->NBoundaries()];
   }

   RNScalar *fewer_than_40s = new RNScalar[nd->NBoundaries()];
   RNScalar *fewer_than_60s = new RNScalar[nd->NBoundaries()];
   RNScalar *fewer_than_80s = new RNScalar[nd->NBoundaries()];

   RNScalar *degree_differences = new RNScalar[nd->NBoundaries()];
   RNScalar *mutual_neighbors = new RNScalar[nd->NBoundaries()];

   // iterate through all boundaries
   if (print_verbose) { printf("Creating boundary statistics...\n  "); fflush(stdout); }
   for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
      if (print_verbose) RNProgressBar(ib, nd->NBoundaries());
      NeuronBoundary *boundary = nd->Boundary(ib);

      // get all of the affinities
      std::vector<RNScalar> affinities = std::vector<RNScalar>();
      for (int ia = 0; ia < boundary->NAffinities(); ++ia)
         affinities.push_back(boundary->Affinity(ia));

      // get mean and stddev of boundaries
      RNScalar mean = boundary->Mean();
      RNScalar stddev = boundary->StdDev();

      // find the skew and kurtosis of the boundaries
      RNScalar skew = 0.0;
      RNScalar kurtosis = 0.0;
      for (unsigned int ia = 0; ia < affinities.size(); ++ia) {
         skew += (affinities[ia] - mean) * (affinities[ia] - mean) * (affinities[ia] - mean);
         kurtosis += (affinities[ia] - mean) * (affinities[ia] - mean) * (affinities[ia] - mean) * (affinities[ia] - mean);
      }
      skew /= affinities.size();
      skew /= (stddev * stddev * stddev);
      kurtosis /= affinities.size();
      kurtosis /= (stddev * stddev * stddev * stddev) - 3;
      skews[ib] = skew;
      kurtosises[ib] = kurtosis;

      // create histogram with 25 bins
      RNScalar histogram[nbins];
      RNScalar hist_min = 0.15;
      RNScalar hist_max = 0.85;

      // initialize histogram to 0
      for (int ih = 0; ih < nbins; ++ih)
         histogram[ih] = 0;

      for (unsigned int ia = 0; ia < affinities.size(); ++ia) {
         RNScalar affinity = affinities[ia];

         // affinity is now between 0 and 1
         if (affinity > hist_max) affinity = 1.0;
         else if (affinity < hist_min) affinity = 0.0;
         else affinity = (affinity - hist_min) / (hist_max - hist_min);

         // get the bin for this affinity
         int bin = affinity * nbins;
         if (bin == nbins) bin--;

         // just checking
         rn_assertion((0 <= bin) && (bin < nbins));
         histogram[bin]++;
      }

      // normalize the bin values
      for (int ih = 0; ih < nbins; ++ih)
         histograms[ih][ib] = histogram[ih] / boundary->NAffinities();

      // calculate the number of affinities lower than 0.2, 0.5, 0.8
      int nfewer_than_40 = 0;
      int nfewer_than_60 = 0;
      int nfewer_than_80 = 0;
      for (unsigned int ia = 0; ia < affinities.size(); ++ia) {
         if (affinities[ia] < 0.20) nfewer_than_40++;
         if (affinities[ia] < 0.50) nfewer_than_60++;
         if (affinities[ia] < 0.80) nfewer_than_80++;
      }

      fewer_than_40s[ib] = nfewer_than_40 / (RNScalar)boundary->NAffinities();
      fewer_than_60s[ib] = nfewer_than_60 / (RNScalar)boundary->NAffinities();
      fewer_than_80s[ib] = nfewer_than_80 / (RNScalar)boundary->NAffinities();

      // find the number of common neighbors
      NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
      NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();

      // get the difference in degrees
      int degree_difference = abs(supervoxel_one->NBoundaries() - supervoxel_two->NBoundaries());
      int mutual_neighbor = 0;

      // go through all boundaries in supervoxel one
      for (int ib1 = 0; ib1 < supervoxel_one->NBoundaries(); ++ib1) {
         NeuronBoundary *boundary_one = supervoxel_one->Boundary(ib1);

         // get the neighbor supervoxel
         NeuronSupervoxel *neighbor = boundary_one->OtherSupervoxel(supervoxel_one);

         // see if this supervoxel neighbors supervoxel two
         for (int ib2 = 0; ib2 < neighbor->NBoundaries(); ++ib2) {
            NeuronBoundary *boundary_two = neighbor->Boundary(ib2);

            // see if this neighbor is supervoxel two
            if (boundary_two->OtherSupervoxel(neighbor) == supervoxel_two) {
               mutual_neighbor++;
            }
         }
      }

      // update degree difference and mutual neighbors
      degree_differences[ib] = degree_difference;
      mutual_neighbors[ib] = mutual_neighbor;

   }
   if (print_verbose) printf("\ndone.\n");

   // save files
   char skew_filename[4096];
   sprintf(skew_filename, "%s/%s_skew.feature", features_directory, root_filename);
   char kurtosis_filename[4096];
   sprintf(kurtosis_filename, "%s/%s_kurtosis.feature", features_directory, root_filename);
   char fewer_40_filename[4096];
   sprintf(fewer_40_filename, "%s/%s_fewer_40.feature", features_directory, root_filename);
   char fewer_60_filename[4096];
   sprintf(fewer_60_filename, "%s/%s_fewer_60.feature", features_directory, root_filename);
   char fewer_80_filename[4096];
   sprintf(fewer_80_filename, "%s/%s_fewer_80.feature", features_directory, root_filename);
   char degree_diff_filename[4096];
   sprintf(degree_diff_filename, "%s/%s_degree_differences.feature", features_directory, root_filename);
   char mutual_filename[4096];
   sprintf(mutual_filename, "%s/%s_mutual_neighbors.feature", features_directory, root_filename);

   // open files
   FILE *skew_fp = fopen(skew_filename, "wb");
   if (!skew_fp) { fprintf(stderr, "Failed to write %s\n", skew_filename); return 0; }
   FILE *kurtosis_fp = fopen(kurtosis_filename, "wb");
   if (!kurtosis_fp) { fprintf(stderr, "Failed to write %s\n", kurtosis_filename); return 0; }
   FILE *fewer_40_fp = fopen(fewer_40_filename, "wb");
   if (!fewer_40_fp) { fprintf(stderr, "Failed to write %s\n", fewer_40_filename); return 0; }
   FILE *fewer_60_fp = fopen(fewer_60_filename, "wb");
   if (!fewer_60_fp) { fprintf(stderr, "Failed to write %s\n", fewer_60_filename); return 0; }
   FILE *fewer_80_fp = fopen(fewer_80_filename, "wb");
   if (!fewer_80_fp) { fprintf(stderr, "Failed to write %s\n", fewer_80_filename); return 0; }
   FILE *degree_fp = fopen(degree_diff_filename, "wb");
   if (!degree_fp) { fprintf(stderr, "Failed to write %s\n", degree_diff_filename); return 0; }
   FILE *mutual_fp = fopen(mutual_filename, "wb");
   if (!mutual_fp) { fprintf(stderr, "Failed to write %s\n", mutual_filename); return 0; }

   int nboundaries = nd->NBoundaries();
   fwrite(&nboundaries, sizeof(int), 1, skew_fp);
   fwrite(skews, sizeof(RNScalar), nd->NBoundaries(), skew_fp);
   fwrite(&nboundaries, sizeof(int), 1, kurtosis_fp);
   fwrite(kurtosises, sizeof(RNScalar), nd->NBoundaries(), kurtosis_fp);
   fwrite(&nboundaries, sizeof(int), 1, fewer_40_fp);
   fwrite(fewer_than_40s, sizeof(RNScalar), nd->NBoundaries(), fewer_40_fp);
   fwrite(&nboundaries, sizeof(int), 1, fewer_60_fp);
   fwrite(fewer_than_60s, sizeof(RNScalar), nd->NBoundaries(), fewer_60_fp);
   fwrite(&nboundaries, sizeof(int), 1, fewer_80_fp);
   fwrite(fewer_than_80s, sizeof(RNScalar), nd->NBoundaries(), fewer_80_fp);
   fwrite(&nboundaries, sizeof(int), 1, degree_fp);
   fwrite(degree_differences, sizeof(RNScalar), nd->NBoundaries(), degree_fp);
   fwrite(&nboundaries, sizeof(int), 1, mutual_fp);
   fwrite(mutual_neighbors, sizeof(RNScalar), nd->NBoundaries(), mutual_fp);

   // close files
   fclose(skew_fp);
   fclose(kurtosis_fp);
   fclose(fewer_40_fp);
   fclose(fewer_60_fp);
   fclose(fewer_80_fp);
   fclose(degree_fp);
   fclose(mutual_fp);

   // save the histogram
   for (int ib = 0; ib < nbins; ++ib) {
      char histogram_filename[4096];
      sprintf(histogram_filename, "%s/%s_histogram_%02d.feature", features_directory, root_filename, ib);

      // open file
      FILE *hist_fp = fopen(histogram_filename, "wb");
      if (!hist_fp) { fprintf(stderr, "Failed to write %s\n", histogram_filename); return 0; }

      fwrite(&nboundaries, sizeof(int), 1, hist_fp);
      fwrite(histograms[ib], sizeof(RNScalar), nd->NBoundaries(), hist_fp);

      // close file
      fclose(hist_fp);
   }

   // free memory
   delete[] skews;
   delete[] kurtosises;
   delete[] fewer_than_40s;
   delete[] fewer_than_60s;
   delete[] fewer_than_80s;
   delete[] degree_differences;
   delete[] mutual_neighbors;
   for (int ib = 0; ib < nbins; ++ib)
      delete[] histograms[ib];

   // return success
   return 1;
}



static int CalculateDerivativeFeatures(char root_filename[4096])
{
   for (int id = 0; id < NDERIVATIVES; ++id) {
      R3Grid **derivative;

      char input_filename[4096];
      sprintf(input_filename, "%s/%s_%s", filter_directory, root_filename, derivative_names[id]);

      // read the derivative file
      derivative = RNReadNeuronMetaRawFile(input_filename, TRUE);
      if (!derivative) exit(-1);
      for (int dim = 0; dim <= 2; ++dim)
         if (!derivative[dim]) exit(-1);

      RNScalar *means = new RNScalar[nd->NBoundaries()];
      RNScalar *medians = new RNScalar[nd->NBoundaries()];
      RNScalar *maximums = new RNScalar[nd->NBoundaries()];
      RNScalar *minimums = new RNScalar[nd->NBoundaries()];
      RNScalar *stddevs = new RNScalar[nd->NBoundaries()];
      RNScalar *skews = new RNScalar[nd->NBoundaries()];
      RNScalar *kurtosises = new RNScalar[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         means[ib] = 0.0;
         medians[ib] = 0.0;
         maximums[ib] = 0.0;
         minimums[ib] = 0.0;
         stddevs[ib] = 0.0;
         skews[ib] = 0.0;
         kurtosises[ib] = 0.0;
      }


      // calculate boundary features
      if (print_verbose) { printf("Calculating results for %s\n  ", derivative_names[id]); fflush(stdout); }
      std::vector<RNScalar> *affinities = new std::vector<RNScalar>[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib)
         affinities[ib] = std::vector<RNScalar>();
      for (int iv = 0; iv < nd->NVoxels(); ++iv) {
         if (print_verbose) RNProgressBar(iv, nd->NVoxels());
         NeuronVoxel *voxel = nd->Voxel(iv);
         NeuronSupervoxel *supervoxel = voxel->Supervoxel();

         // go through all neighbors
         for (int in = 0; in < voxel->NNeighbors(); ++in) {
            NeuronVoxel *neighbor = voxel->Neighbor(in);
            if (!neighbor) continue;
            if (voxel->DataIndex() < neighbor->DataIndex()) continue;


            NeuronSupervoxel *neighbor_supervoxel = neighbor->Supervoxel();
            if (neighbor_supervoxel == supervoxel) continue;

            // find boundary
            int boundary_index = -1;
            for (int ib = 0; ib < supervoxel->NBoundaries(); ++ib) {
               NeuronBoundary *boundary = supervoxel->Boundary(ib);
               if (boundary->OtherSupervoxel(supervoxel) == neighbor_supervoxel) {
                  boundary_index = boundary->DataIndex();
                  break;
               }
            }

            rn_assertion(boundary_index != -1);
            RNScalar affinity;
            if (in == 0) affinity = derivative[RN_X]->GridValue(voxel->DataIndex());
            if (in == 1) affinity = derivative[RN_X]->GridValue(neighbor->DataIndex());
            if (in == 2) affinity = derivative[RN_Y]->GridValue(voxel->DataIndex());
            if (in == 3) affinity = derivative[RN_Y]->GridValue(neighbor->DataIndex());
            if (in == 4) affinity = derivative[RN_Z]->GridValue(voxel->DataIndex());
            if (in == 5) affinity = derivative[RN_Z]->GridValue(neighbor->DataIndex());

            affinities[boundary_index].push_back(affinity);
         }
      }
      if (print_verbose) printf("\ndone\n");

      if (print_verbose) { printf("Creating distributions for %s\n", derivative_names[id]); fflush(stdout); }
      for (int ib = 0; ib < nd->NBoundaries(); ++ib) {
         if (print_verbose) RNProgressBar(ib, nd->NBoundaries());
         // get the distribution for the affinities
         RNDistribution distribtion = RNDistribution(affinities[ib]);
         NeuronBoundary *boundary = nd->Boundary(ib);

         means[boundary->DataIndex()] = distribtion.Mean();
         medians[boundary->DataIndex()] = distribtion.Median();
         maximums[boundary->DataIndex()] = distribtion.Maximum();
         minimums[boundary->DataIndex()] = distribtion.Minimum();
         stddevs[boundary->DataIndex()] = distribtion.StdDev();
         skews[boundary->DataIndex()] = distribtion.Skew();
         kurtosises[boundary->DataIndex()] = distribtion.Kurtosis();
      }
      if (print_verbose) printf("\ndone.\n");

      // save files
      char mean_filename[4096];
      sprintf(mean_filename, "%s/%s_%s_mean.feature", features_directory, root_filename, derivative_names[id]);
      char median_filename[4096];
      sprintf(median_filename, "%s/%s_%s_median.feature", features_directory, root_filename, derivative_names[id]);
      char maximum_filename[4096];
      sprintf(maximum_filename, "%s/%s_%s_maximum.feature", features_directory, root_filename, derivative_names[id]);
      char minimum_filename[4096];
      sprintf(minimum_filename, "%s/%s_%s_minimum.feature", features_directory, root_filename, derivative_names[id]);
      char stddev_filename[4096];
      sprintf(stddev_filename, "%s/%s_%s_stddev.feature", features_directory, root_filename, derivative_names[id]);
      char skew_filename[4096];
      sprintf(skew_filename, "%s/%s_%s_skew.feature", features_directory, root_filename, derivative_names[id]);
      char kurtosis_filename[4096];
      sprintf(kurtosis_filename, "%s/%s_%s_kurtosis.feature", features_directory, root_filename, derivative_names[id]);

      // open files
      FILE *mean_fp = fopen(mean_filename, "wb");
      if (!mean_fp) { fprintf(stderr, "Failed to write %s\n", mean_filename); return 0; }
      FILE *median_fp = fopen(median_filename, "wb");
      if (!median_fp) { fprintf(stderr, "Failed to write %s\n", median_filename); return 0; }
      FILE *maximum_fp = fopen(maximum_filename, "wb");
      if (!maximum_fp) { fprintf(stderr, "Failed to write %s\n", maximum_filename); return 0; }
      FILE *minimum_fp = fopen(minimum_filename, "wb");
      if (!minimum_fp) { fprintf(stderr, "Failed to write %s\n", minimum_filename); return 0; }
      FILE *stddev_fp = fopen(stddev_filename, "wb");
      if (!stddev_fp) { fprintf(stderr, "Failed to write %s\n", stddev_filename); return 0; }
      FILE *skew_fp = fopen(skew_filename, "wb");
      if (!skew_fp) { fprintf(stderr, "Failed to write %s\n", skew_filename); return 0; }
      FILE *kurtosis_fp = fopen(kurtosis_filename, "wb");
      if (!kurtosis_fp) { fprintf(stderr, "Failed to write %s\n", kurtosis_filename); return 0; }

      // write boundaries
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, mean_fp);
      fwrite(&nboundaries, sizeof(int), 1, median_fp);
      fwrite(&nboundaries, sizeof(int), 1, maximum_fp);
      fwrite(&nboundaries, sizeof(int), 1, minimum_fp);
      fwrite(&nboundaries, sizeof(int), 1, stddev_fp);
      fwrite(&nboundaries, sizeof(int), 1, skew_fp);
      fwrite(&nboundaries, sizeof(int), 1, kurtosis_fp);

      // write results
      fwrite(means, sizeof(RNScalar), nd->NBoundaries(), mean_fp);
      fwrite(medians, sizeof(RNScalar), nd->NBoundaries(), median_fp);
      fwrite(maximums, sizeof(RNScalar), nd->NBoundaries(), maximum_fp);
      fwrite(minimums, sizeof(RNScalar), nd->NBoundaries(), minimum_fp);
      fwrite(stddevs, sizeof(RNScalar), nd->NBoundaries(), stddev_fp);
      fwrite(skews, sizeof(RNScalar), nd->NBoundaries(), skew_fp);
      fwrite(kurtosises, sizeof(RNScalar), nd->NBoundaries(), kurtosis_fp);

      // close file
      fclose(mean_fp);
      fclose(median_fp);
      fclose(maximum_fp);
      fclose(minimum_fp);
      fclose(stddev_fp);
      fclose(skew_fp);
      fclose(kurtosis_fp);

      // free memory
      for (int dim = 0; dim <= 2; ++dim) {
         delete derivative[dim];
      }
      delete[] derivative;
      delete[] means;
      delete[] medians;
      delete[] maximums;
      delete[] minimums;
      delete[] stddevs;
      delete[] skews;
      delete[] kurtosises;
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

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create boundary features
   if (!CreateBoundaryFeatures(root_filename)) exit(-1);
   if (!CalculateDerivativeFeatures(root_filename)) exit(-1);

   // free up memory
   delete nd;

   // return success
   return 0;
}
