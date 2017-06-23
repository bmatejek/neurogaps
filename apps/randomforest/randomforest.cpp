// Source file for the simple file conversion algorithm



// include files

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"
#include <algorithm>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static int validation = 0;
static int smoothing_terms = 0;
static int data_terms = 0;
static int boundary_terms = 0;
// feature options
static int boundary_features = 0;
static int shape_features = 0;
static int global_features = 0;
// string arguments
static const char *extension = NULL;
static char *random_forest_filename = NULL;
static char *testing_input_filename = NULL;



// global data sets and random forest

static NeuronData *testing_nd = NULL;
static RNDataset *testing_dataset = NULL;
static RNRandomForest *random_forest = NULL;



// directory structure
static const char *dataset_directory = "algs_data/dataset";
static const char *graphcut_directory = "algs_data/graphcut";
static const char *hierarchical_directory = "algs_data/hierarchical";

// output directories where the results are stored
static const char *results_directory = "results/random_forest";
static const char *tmp_results_directory = "results/random_forest/tmp";
static const char *visuals_directory = "visuals/random_forest";

// input directories where random forests are stored
static const char *random_forest_directory = "algs_data/random_forest";
static const char *tmp_random_forest_directory = "algs_data/random_forest/tmp";



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



////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////

// structure encapsulates random forest entry data
struct ForestResult
{
   ForestResult(void) {
      this->boundary = NULL;
      this->prediction = FLT_MAX;
      this->match = FALSE;
   }

   ForestResult(NeuronBoundary *boundary, RNScalar prediction, RNBoolean match) {
      this->boundary = boundary;
      this->prediction = prediction;
      this->match = match;
   }

   NeuronBoundary *boundary;
   RNScalar prediction;
   RNBoolean match;
};


// compare function for random forest results
int ForestResultMaxCompare(ForestResult one, ForestResult two)
{
   return one.prediction > two.prediction;
}



int ForestResultMinCompare(ForestResult one, ForestResult two)
{
   return one.prediction < two.prediction;
}



static int PrecisionRecall(NeuronData *nd, RNDataset *dataset, std::vector<float> &results)
{
   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // create title qualifier
   char title_qualifier[4096];
   if (!strcmp(extension, "smoothing")) sprintf(title_qualifier, "Smoothing");
   else if (!strcmp(extension, "data")) sprintf(title_qualifier, "Data");
   else if (!strcmp(extension, "boundary")) sprintf(title_qualifier, "Boundary");
   else if (!strcmp(extension, "shape")) sprintf(title_qualifier, "Shape");
   else if (!strcmp(extension, "data_boundary")) sprintf(title_qualifier, "Boundary Data");
   else if (!strcmp(extension, "global")) sprintf(title_qualifier, "Global");
   else if (!strcmp(extension, "smoothing_global")) sprintf(title_qualifier, "Smoothing Global");
   else { fprintf(stderr, "Unrecognized extension: %s\n", extension); return 0; }

   // count the correct number of true and false results
   unsigned int positive_counter = 0;
   unsigned int negative_counter = 0;

   // create precision and recall curves for this random forest
   std::vector<ForestResult> predicted_positive_forest_results = std::vector<ForestResult>(); // only consider positive predictions
   std::vector<ForestResult> predicted_negative_forest_results = std::vector<ForestResult>(); // only consider negative predictions
   std::vector<ForestResult> positive_forest_results = std::vector<ForestResult>();           // sort by positive probability

   for (unsigned int ie = 0; ie < dataset->NEntries(); ++ie) {
      float prediction = results[ie];
      unsigned int label = dataset->GetLabel(ie);
      unsigned int identification = dataset->GetIdentification(ie);
      RNBoolean match = ((prediction > 0.5) == label);

      // set NeuronBoundary for smoothing terms
      NeuronBoundary *boundary = NULL;
      if (smoothing_terms) boundary = nd->Boundary(identification);

      // increment based on the number of true and false occurrences
      if (label) positive_counter++;
      else negative_counter++;

      // predicted true
      if (prediction > 0.5) {
         if (smoothing_terms) predicted_positive_forest_results.push_back(ForestResult(boundary, prediction, match));
         else predicted_negative_forest_results.push_back(ForestResult(NULL, prediction, match));
      }
      // predicted false
      else {
         if (smoothing_terms) predicted_negative_forest_results.push_back(ForestResult(boundary, 1.0 - prediction, match));
         else predicted_positive_forest_results.push_back(ForestResult(NULL, 1.0 - prediction, match));
      }

      // add to full result vectors
      if (smoothing_terms) positive_forest_results.push_back(ForestResult(boundary, prediction, match));
      else positive_forest_results.push_back(ForestResult(NULL, 1.0 - prediction, match));
   }

   // sort the vectors
   std::sort(predicted_positive_forest_results.begin(), predicted_positive_forest_results.end(), &ForestResultMaxCompare);
   std::sort(predicted_negative_forest_results.begin(), predicted_negative_forest_results.end(), &ForestResultMaxCompare);
   std::sort(positive_forest_results.begin(), positive_forest_results.end(), &ForestResultMaxCompare);

   // create output images
   RNPyImage *precision_recall = new RNPyImage();
   precision_recall->SetOutputDirectory(visuals_directory);

   // create point plots
   RNPyPointPlot *positive_precision = new RNPyPointPlot();
   RNPyPointPlot *negative_precision = new RNPyPointPlot();
   RNPyPointPlot *proportion_plot = new RNPyPointPlot();

   // set titles
   char positive_title[4096];
   sprintf(positive_title, "%s %s Positive Precision - Recall", root_filename, title_qualifier);
   positive_precision->SetTitle(positive_title);

   char negative_title[4096];
   sprintf(negative_title, "%s %s Negative Precision - Recall", root_filename, title_qualifier);
   negative_precision->SetTitle(negative_title);

   char proportion_title[4096];
   sprintf(proportion_title, "%s %s Merge Score versus Proportion Positive", root_filename, title_qualifier);
   proportion_plot->SetTitle(proportion_title);

   // set X and Y axis labels
   positive_precision->SetXLabel("Recall");
   positive_precision->SetYLabel("Precision");
   negative_precision->SetXLabel("Recall");
   negative_precision->SetYLabel("Precision");
   proportion_plot->SetXLabel("Data Point Score");
   proportion_plot->SetYLabel("Proportion Correct");

   // force x and y axes
   positive_precision->SetXAxisMax(1.0);
   positive_precision->SetXAxisMin(0.0);
   positive_precision->SetYAxisMax(1.0);
   positive_precision->SetYAxisMin(0.0);
   negative_precision->SetXAxisMax(1.0);
   negative_precision->SetXAxisMin(0.0);
   negative_precision->SetYAxisMax(1.0);
   negative_precision->SetYAxisMin(0.0);
   proportion_plot->SetXAxisMax(1.0);
   proportion_plot->SetXAxisMin(0.0);
   proportion_plot->SetYAxisMax(1.0);
   proportion_plot->SetYAxisMin(0.0);

   unsigned int ntrue_positives = 0;
   unsigned int nfalse_positives = 0;
   for (unsigned int ie = 0; ie < predicted_positive_forest_results.size(); ++ie) {
      // get this data entry result
      ForestResult result = predicted_positive_forest_results[ie];
      if (result.match) ntrue_positives++;
      else nfalse_positives++;

      // calculate the precision and recall for this point
      RNScalar precision = ntrue_positives / (RNScalar)(ntrue_positives + nfalse_positives);
      RNScalar recall = ntrue_positives / (RNScalar)positive_counter;

      // add to the plot
      positive_precision->InsertPoint(R2Point(recall, precision));
   }
   unsigned int ntrue_negatives = 0;
   unsigned int nfalse_negatives = 0;
   for (unsigned int ie = 0; ie < predicted_negative_forest_results.size(); ++ie) {
      // get this data entry result
      ForestResult result = predicted_negative_forest_results[ie];
      if (result.match) ntrue_negatives++;
      else nfalse_positives++;

      /// calculate the precision and recall for this point
      RNScalar precision = ntrue_negatives / (RNScalar)(ntrue_negatives + nfalse_negatives);
      RNScalar recall = ntrue_negatives / (RNScalar)negative_counter;

      // add to the plot
      negative_precision->InsertPoint(R2Point(recall, precision));
   }

   // create proportion histogram
   int nbins = 100;
   unsigned int *ncorrect = new unsigned int[nbins];
   unsigned int *nincorrect = new unsigned int[nbins];
   for (int ib = 0; ib < nbins; ++ib) {
      ncorrect[ib] = 0;
      nincorrect[ib] = 0;
   }

   // calculate the number of correct positive entries
   for (unsigned int ie = 0; ie < positive_forest_results.size(); ++ie) {
      ForestResult result = positive_forest_results[ie];
      RNScalar prediction = result.prediction;
      int bin = (int)(nbins * prediction + 0.5);

      if (result.match) ncorrect[bin]++;
      else nincorrect[bin]++;
   }

   // add points to proportion plot
   for (int ib = 0; ib < nbins; ++ib) {
      if (ncorrect[ib] + nincorrect[ib] != 0) {
         RNScalar prediction = ib / (RNScalar)nbins;
         RNScalar proportion = ncorrect[ib] / (RNScalar)(ncorrect[ib] + nincorrect[ib]);
         proportion_plot->InsertPoint(R2Point(prediction, proportion));
      }
   }

   // add the plots to the images
   precision_recall->InsertPyPlot(positive_precision);
   precision_recall->InsertPyPlot(negative_precision);
   precision_recall->InsertPyPlot(proportion_plot);

   // create output filename
   char precision_recall_filename[4096];
   sprintf(precision_recall_filename, "%s/%s_%s_precision_recall.pyimage", results_directory, root_filename, extension);

   // write the image
   precision_recall->WriteImageFile(precision_recall_filename);

   // save merge order
   if (smoothing_terms) {
      // get file
      char hierarchical_merge_filename[4096];
      sprintf(hierarchical_merge_filename, "%s/%s_random_forest_%s.hier", hierarchical_directory, root_filename, extension);

      // open file
      FILE *hier_fp = fopen(hierarchical_merge_filename, "wb");
      if (!hier_fp) { fprintf(stderr, "Failed to write to %s\n", hierarchical_merge_filename); return 0; }

      // write the number of boundaries
      int nboundaries = positive_forest_results.size();
      fwrite(&nboundaries, sizeof(int), 1, hier_fp);

      // write all of the boundaries
      for (unsigned int ie = 0; ie < positive_forest_results.size(); ++ie) {
         ForestResult result = positive_forest_results[ie];
         NeuronBoundary *boundary = result.boundary;
         rn_assertion(boundary != NULL);

         // write this boundary to file
         int boundary_data_index = boundary->DataIndex();
         fwrite(&boundary_data_index, sizeof(int), 1, hier_fp);
      }

      // close file
      fclose(hier_fp);
   }

   // free memory
   delete positive_precision;
   delete negative_precision;
   delete proportion_plot;
   delete[] ncorrect;
   delete[] nincorrect;
   delete precision_recall;

   // return success
   return 1;
}



static int OptimizeThreshold(NeuronData *nd, RNDataset *dataset, std::vector<float> &results)
{
   // get the root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   std::vector<ForestResult> forest_results = std::vector<ForestResult>();

   // go through all boundaries
   for (unsigned int ie = 0; ie < dataset->NEntries(); ++ie) {
      // get the identification for this boundary
      int identification = dataset->GetIdentification(ie);

      // get the boundary and accompanying supervoxels
      NeuronBoundary *boundary = NULL;
      NeuronSupervoxel *supervoxel_one = NULL;
      NeuronSupervoxel *supervoxel_two = NULL;
      if (smoothing_terms) {
         boundary = nd->Boundary(identification);
         supervoxel_one = boundary->SupervoxelOne();
         supervoxel_two = boundary->SupervoxelTwo();
      }
      else {
         int ic1 = identification / nd->NCellulars();
         int ic2 = identification % nd->NCellulars();

         supervoxel_one = nd->Supervoxel(ic1);
         supervoxel_two = nd->Supervoxel(ic2);
      }
      rn_assertion(supervoxel_one->IsCellular());
      rn_assertion(supervoxel_two->IsCellular());

      // get the predicted label
      float result = results[ie];

      // create forest result and add to vector
      ForestResult forest_result;
      if (smoothing_terms) forest_result = ForestResult(boundary, result, supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel());
      else forest_result = ForestResult(boundary, 1.0 - result, supervoxel_one->MajorityHumanLabel() == supervoxel_two->MajorityHumanLabel());
      forest_results.push_back(forest_result);
   }

   // sort the results
   sort(forest_results.begin(), forest_results.end(), ForestResultMinCompare);

   // see the results of cutting here
   int *left_negatives = new int[dataset->NEntries()];
   int *right_positives = new int[dataset->NEntries()];
   for (unsigned int ie = 0; ie < dataset->NEntries(); ++ie) {
      left_negatives[ie] = 0;
      right_positives[ie] = 0;
   }

   // count the number of negatives boundaries with smaller smoothing term values
   for (int ie = 1; ie < (int)dataset->NEntries(); ++ie) {
      if (!forest_results[ie - 1].match) left_negatives[ie] = left_negatives[ie - 1] + 1;
      else left_negatives[ie] = left_negatives[ie - 1];
   }

   // count the number of positive boundaries with larger smoothing term values
   for (int ie = (int)dataset->NEntries() - 2; ie >= 0; --ie) {
      if (forest_results[ie + 1].match) right_positives[ie] = right_positives[ie + 1] + 1;
      else right_positives[ie] = right_positives[ie + 1];
   }

   RNPyImage *image = new RNPyImage();
   image->SetOutputDirectory(visuals_directory);

   char title[4096];
   sprintf(title, "%s Performance With Thresholding", root_filename);

   // create plot
   RNPyPointPlot *threshold = new RNPyPointPlot();
   threshold->SetXLabel("Threshold");
   threshold->SetYLabel("Proportion Correct");
   threshold->SetTitle(title);

   threshold->SetLegend("Proportion Correct Guesses");
   threshold->SetExtremaType(MAX_EXTREMA);

   // force x and y axis bounds
   threshold->SetXAxisMax(1.0);
   threshold->SetXAxisMin(0.0);
   threshold->SetYAxisMax(1.0);
   threshold->SetYAxisMin(0.0);

   for (int ie = 0; ie < (int)dataset->NEntries(); ++ie) {
      int ncorrect = right_positives[ie] + left_negatives[ie];
      if (forest_results[ie].match) ncorrect++;
      threshold->InsertPoint(R2Point(forest_results[ie].prediction, ncorrect / (RNScalar)forest_results.size()));
   }

   image->InsertPyPlot(threshold);

   // get output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s_threshold_%s.pyimage", results_directory, root_filename, extension);

   image->WriteImageFile(output_filename);

   // free memory
   delete[] left_negatives;
   delete[] right_positives;
   delete threshold;
   delete image;

   // return success
   return 1;
}



static int WriteGraphcutFiles(NeuronData *nd, RNDataset *dataset)
{
   // get the vector results
   std::vector<float> results = random_forest->PredictBinaryProbabilities(*dataset);

   // get root filenames
   char root_filename[4096];
   strcpy(root_filename, nd->Filename());
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // ids correspond to boundaries
   if (smoothing_terms) {
      // allocate memory for array
      float *smoothing_terms = new float[nd->NBoundaries()];
      for (int ib = 0; ib < nd->NBoundaries(); ++ib)
         // default is no smoothing penalty
         smoothing_terms[ib] = 0.0;

      // go through all dataset entries
      for (unsigned int in = 0; in < dataset->NEntries(); ++in) {
         // get boundary identification
         int identification = dataset->GetIdentification(in);

         // get the boundary
         NeuronBoundary *boundary = nd->Boundary(identification);
         NeuronSupervoxel *supervoxel_one = boundary->SupervoxelOne();
         NeuronSupervoxel *supervoxel_two = boundary->SupervoxelTwo();
         rn_assertion(supervoxel_one->IsCellular() && supervoxel_two->IsCellular());

         // set the results
         smoothing_terms[identification] = results[in];
      }

      // get output filename
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_random_forest_%s.binary", graphcut_directory, root_filename, extension);

      FILE *binary_fp = fopen(output_filename, "wb");
      if (!binary_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

      // write boundary terms
      int nboundaries = nd->NBoundaries();
      fwrite(&nboundaries, sizeof(int), 1, binary_fp);
      fwrite(smoothing_terms, sizeof(float), nd->NBoundaries(), binary_fp);

      fclose(binary_fp);

      // free memory
      delete[] smoothing_terms;
   }
   else if (data_terms) {
      // allocate data term array
      int ncellulars = nd->NCellulars();
      float **data_terms = new float *[ncellulars];
      for (int ic1 = 0; ic1 < ncellulars; ++ic1) {
         data_terms[ic1] = new float[ncellulars];
         for (int ic2 = 0; ic2 < ncellulars; ++ic2) {
            data_terms[ic1][ic2] = FLT_MIN;
         }
      }

      // go through every entry in training dataset
      for (unsigned int in = 0; in < dataset->NEntries(); ++in) {
         // corresponds to boundary
         int identification = dataset->GetIdentification(in);

         int ic1 = identification / ncellulars;
         int ic2 = identification % ncellulars;
         rn_assertion(ic1 < nd->NCellulars());
         rn_assertion(ic2 < nd->NCellulars());

         data_terms[ic1][ic2] = 1.0 - results[in];
         data_terms[ic2][ic1] = 1.0 - results[in];
      }

      // get output filenames
      char output_filename[4096];
      sprintf(output_filename, "%s/%s_random_forest.unary", graphcut_directory, root_filename);

      // open files
      FILE *data_fp = fopen(output_filename, "wb");
      if (!data_fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

      fwrite(&ncellulars, sizeof(int), 1, data_fp);
      for (int ic = 0; ic < ncellulars; ++ic) {
         fwrite(data_terms[ic], sizeof(float), ncellulars, data_fp);
      }

      // close files
      fclose(data_fp);

      // free memory
      for (int ic = 0; ic < ncellulars; ++ic)
         delete[] data_terms[ic];
      delete[] data_terms;
   }
   else { fprintf(stderr, "Neither data or smoothing term selected.\n"); return 0; }

   // create precision and recall terms
   //if (!PrecisionRecall(nd, dataset, results)) return 0;
   //if (!OptimizeThreshold(nd, dataset, results)) return 0;

   // return success
   return 1;
}



static int ReadDataset(void)
{
   // get data root filename
   char testing_root_filename[4096];
   strncpy(testing_root_filename, testing_nd->Filename(), 4096);
   char *testing_extp = strrchr(testing_root_filename, '.');
   *testing_extp = '\0';

   // get dataset filenames
   char testing_dataset_filename[4096];
   sprintf(testing_dataset_filename, "%s/%s_%s.mldb", dataset_directory, testing_root_filename, extension);

   // train random forest classifier
   char input_filename[4096];
   if (validation) sprintf(input_filename, "%s/%s_%s.xml", tmp_random_forest_directory, random_forest_filename, extension);
   else sprintf(input_filename, "%s/%s_%s.xml", random_forest_directory, random_forest_filename, extension);

   // load the random forest
   random_forest = new RNRandomForest();
   random_forest->LoadClassifier(input_filename);

   // read in the testing dataset
   testing_dataset = new RNDataset();
   testing_dataset->ReadFile(testing_dataset_filename);
   if (!validation && !WriteGraphcutFiles(testing_nd, testing_dataset)) exit(-1);

   // print statistics
   printf("Testing dataset results for %s terms:\n", extension);
   unsigned int **confusion_matrix = random_forest->PrintPerformance(*testing_dataset);

   if (validation) {
      // save the confusion matrix
      char results_filename[4096];
      sprintf(results_filename, "%s/%s_%s_%s.confusion", tmp_results_directory, testing_root_filename, random_forest_filename, extension);

      // open file
      FILE *results_fp = fopen(results_filename, "wb");
      if (!results_fp) { fprintf(stderr, "Failed to open results filename %s\n", results_filename); return 0; }

      // write the dimension of the confusion matrix
      int ndimensions = testing_dataset->NClasses();
      fwrite(&ndimensions, sizeof(int), 1, results_fp);

      // write the confusion matrix
      for (int i = 0; i < ndimensions; ++i) {
         for (int j = 0; j < ndimensions; ++j) {
            fwrite(&(confusion_matrix[i][j]), sizeof(int), 1, results_fp);
         }
      }

      // close file
      fclose(results_fp);
   }

   // free memory
   for (unsigned int ic = 0; ic < testing_dataset->NClasses(); ++ic)
      delete[] confusion_matrix[ic];
   delete[] confusion_matrix;
   delete testing_dataset;
   delete random_forest;

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
         else if (!strcmp(*argv, "-forest")) { argv++; argc--; random_forest_filename = *argv; }
         else if (!strcmp(*argv, "-test")) { argv++; argc--; testing_input_filename = *argv; }
         else if (!strcmp(*argv, "-validation")) { validation = 1; }
         else if (!strcmp(*argv, "-smoothing_term")) { smoothing_terms = 1; }
         else if (!strcmp(*argv, "-data_term")) { data_terms = 1; }
         else if (!strcmp(*argv, "-boundary_term")) { boundary_terms = 1; }
         // feature options
         else if (!strcmp(*argv, "-boundary_features")) { boundary_features = 1; }
         else if (!strcmp(*argv, "-shape_features")) { shape_features = 1; }
         else if (!strcmp(*argv, "-global_features")) { global_features = 1; }
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
	 if (!testing_input_filename) testing_input_filename = *argv;
	 else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!random_forest_filename) {
      fprintf(stderr, "Need to supply random forest filename.\n");
      return 0;
   }
   if (!testing_input_filename) {
      fprintf(stderr, "Need to supply testing input filename.\n");
      return 0;
   }

   // make sure there are consistent labeling options
   if (data_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      extension = "data";
   }
   else if (boundary_terms) {
      if (boundary_features || shape_features || global_features) {
         fprintf(stderr, "Cannot select boundary, shape, or global features with data terms\n");
         return 0;
      }

      data_terms = 1;
      extension = "data_boundary";
   }
   else if (smoothing_terms) {
      if (global_features) extension = "smoothing_global";
      else extension = "smoothing";

      // turn on all feature options besides global
      boundary_features = 1;
      shape_features = 1;
   }
   else if (boundary_features || shape_features || global_features) {
      // make sure only one feature is selected
      if (boundary_features + shape_features + global_features > 1) {
         fprintf(stderr, "Can only select one among boundary, shape, and global features.\n");
         return 0;
      }

      // get the correct extension
      if (boundary_features) extension = "boundary";
      if (shape_features) extension = "shape";
      if (global_features) extension = "global";

      // turn on smoothing terms to go to the correct functions
      smoothing_terms = 1;
   }
   else {
      fprintf(stderr, "Must select either data or smoothing terms.\n");
      return 0;
   }

   // additional input verification
   char *extp = strrchr(random_forest_filename, '.');
   if (extp) { fprintf(stderr, "Root filename should not have an extension.\n"); return 0; }
   char *prep = strrchr(random_forest_filename, '/');
   if (prep) { fprintf(stderr, "Root filename should not include directory structure.\n"); return 0; }

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

   // read in neuron files
   testing_nd = ReadData(testing_input_filename);
   if (!testing_nd) exit(-1);

   if (!ReadDataset()) exit(-1);

   // free up memory
   delete testing_nd;

   // return success
   return 0;
}
