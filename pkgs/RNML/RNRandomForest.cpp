// Source file for the machine learning random forest classifier class



//////////////////////////////////////////////////////////////
// Include files
//////////////////////////////////////////////////////////////

#include "RNML.h"



/* class type definitions */

RN_CLASS_TYPE_DEFINITIONS(RNRandomForest);



/////////////////////////////////////////////////////////////
// Constructors/destructors
/////////////////////////////////////////////////////////////

RNRandomForest::
RNRandomForest(void) :
random_forest(NULL)
{
   sprintf(classifier_type_name, "random_forest");
}



RNRandomForest::
RNRandomForest(RNDataset &data, int max_depth = 5, int min_sample_count = 1, RNBoolean rebalance = FALSE, int nactive_vars = 0, int ntrees = 1000, float forest_accuracy = 0.01) :
random_forest(NULL)
{
   sprintf(classifier_type_name, "random_forest");

   // assert that the data set is labeled
   rn_assertion(data.IsLabeled());

   // define training and classification matricies
   cv::Mat *training_data = data.OpenCVDataMatrix();
   cv::Mat *training_classifications = data.OpenCVClassiciationsMatrix();

   // define all attributes as numerical
   cv::Mat *var_type = data.OpenCVAttributesMatrix();

   // determine if rebalancing is needed
   float *priors = NULL;
   if (rebalance) {
      // calculate the number of occurrences
      int *occurrences = new int[data.NClasses()];
      for (unsigned int ic = 0; ic < data.NClasses(); ++ic)
         occurrences[ic] = 0;
      for (unsigned int ie = 0; ie < data.NEntries(); ++ie)
         occurrences[data.GetLabel(ie)]++;

      // create priors array
      priors = new float[data.NClasses()];
      for (unsigned int ic = 0; ic < data.NClasses(); ++ic)
         priors[ic] = occurrences[ic] / (float)data.NEntries();

      // free memory
      delete[] occurrences;
   }

   // use parameters input from 
   float regression_accuracy = 0.0;  // not needed for classification
   bool use_surrogates = false;      // no missing data entries allowed
   int max_categories = 10;          // not used for 2-class classification problems
   bool calc_var_importance = true;  // for internal use
   // terminate learning by forest accuracy
   CvRTParams params = CvRTParams(max_depth, min_sample_count, regression_accuracy, use_surrogates, max_categories, priors, calc_var_importance, nactive_vars, ntrees, forest_accuracy, CV_TERMCRIT_ITER);

   // create
   random_forest = new CvRTrees();

   // train the random forest
   printf("Training the random forest classifier...\n"); fflush(stdout);
   RNTime train_time;
   train_time.Read();
   random_forest->train(*training_data, CV_ROW_SAMPLE, *training_classifications, cv::Mat(), cv::Mat(), *var_type, cv::Mat(), params);

   // output the performance on the training dataset
   printf("  Tree Count: %d\n", random_forest->get_tree_count());
   printf("done in %0.2f seconds!\n", train_time.Elapsed());

   // delete excess memory
   delete training_data;
   delete training_classifications;
   delete var_type;
   if (priors) delete[] priors;
}



RNRandomForest::
~RNRandomForest(void)
{
   delete random_forest;
}



////////////////////////////////////////////////////////////
// Classification functions
////////////////////////////////////////////////////////////

unsigned int RNRandomForest::
Classify(const std::vector<float> &attributes) const
{
   // create the testing feature attributes
   cv::Mat *testing_attributes = new cv::Mat(1, attributes.size(), CV_32FC1);
   for (unsigned int ia = 0; ia < attributes.size(); ++ia) {
      testing_attributes->at<float>(0, ia) = attributes[ia];
   }

   // calculate the result
   float result = random_forest->predict(*testing_attributes, cv::Mat());

   // return the best integer guess
   return (unsigned int)(0.5 + result);
}



float RNRandomForest::
PredictBinaryProbability(const std::vector<float> &attributes) const
{
   // create the testing matrix
   cv::Mat *testing_attributes = new cv::Mat(1, attributes.size(), CV_32FC1);
   for (unsigned int ia = 0; ia < attributes.size(); ++ia) {
      testing_attributes->at<float>(0, ia) = attributes[ia];
   }

   // determine the probability for this feature vector
   float probability = random_forest->predict_prob(*testing_attributes, cv::Mat());

   delete testing_attributes;

   // return the probability
   return probability;
}



std::vector<unsigned int> RNRandomForest::
Classify(RNDataset& testing_data) const
{
   std::vector<unsigned int> results = std::vector<unsigned int>();

   for (unsigned int in = 0; in < testing_data.NEntries(); ++in) {
      std::vector<float> attributes = testing_data.GetFeature(in);
      results.push_back(Classify(attributes));
   }

   // return the vector of results
   return results;
}



std::vector<float> RNRandomForest::
PredictBinaryProbabilities(RNDataset& testing_data) const
{
   std::vector<float> results = std::vector<float>();

   printf("Calculating binary probabilities...\n  "); fflush(stdout);
   for (unsigned int in = 0; in < testing_data.NEntries(); ++in) {
      RNProgressBar(in, testing_data.NEntries());
      std::vector<float> attributes = testing_data.GetFeature(in);
      results.push_back(PredictBinaryProbability(attributes));
   }
   printf("\ndone.\n");

   // return the vector of results
   return results;
}



////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////

unsigned int **RNRandomForest::
PrintPerformance(RNDataset& data) const
{
   // make sure the data is labeled
   rn_assertion(data.IsLabeled());

   // define training and classification matricies
   cv::Mat *testing_data = data.OpenCVDataMatrix();

   unsigned int **confusion_matrix = new unsigned int *[data.NClasses()];
   for (unsigned int i = 0; i < data.NClasses(); ++i) {
      confusion_matrix[i] = new unsigned int[data.NClasses()];
      for (unsigned int j = 0; j < data.NClasses(); ++j) {
         confusion_matrix[i][j] = 0;
      }
   }

   for (unsigned int ie = 0; ie < data.NEntries(); ++ie) {
      cv::Mat testing_sample = testing_data->row(ie);

      float result = random_forest->predict(testing_sample, cv::Mat());
      int categorical_result = (int)(0.5 + result);
      int label = data.GetLabel(ie);
      confusion_matrix[label][categorical_result]++;
   }
   printf("------------------------------------------------------\n");
   printf("Confusion Matrix:\n");
   printf("                  Predicted Class\n");
   printf("Actual Class");
   for (unsigned int i = 0; i < data.NClasses(); ++i) {
      printf("%10d", i);
   }
   printf("\n");
   for (unsigned int i = 0; i < data.NClasses(); ++i) {
      printf("%12d", i);
      for (unsigned int j = 0; j < data.NClasses(); ++j) {
         printf("%10u", confusion_matrix[i][j]);
      }
      printf("\n");
   }

   unsigned int total_correct = 0;
   unsigned int *total_per_class = new unsigned int[data.NClasses()];
   for (unsigned int i = 0; i < data.NClasses(); ++i) {
      total_correct += confusion_matrix[i][i];
      total_per_class[i] = 0;
      for (unsigned int j = 0; j < data.NClasses(); ++j) {
         total_per_class[i] += confusion_matrix[i][j];
      }
   }
   printf("------------------------------------------------------\n");
   printf("Correctly Labeled: %lf%%\n", 100 * (RNScalar)total_correct / data.NEntries());
   for (unsigned int i = 0; i < data.NClasses(); ++i)
      printf("Class %d Correctly Labeled: %lf%%\n", i, 100 * (RNScalar)confusion_matrix[i][i] / total_per_class[i]);
   printf("------------------------------------------------------\n");

   delete testing_data;
   delete[] total_per_class;

   return confusion_matrix;
}



void RNRandomForest::
SaveClassifier(const char output_filename[4096]) const
{
   // save the random forest
   random_forest->save(output_filename);
}



void RNRandomForest::
LoadClassifier(const char input_filename[4096])
{
   // load the random forest
   random_forest = new CvRTrees();
   random_forest->load(input_filename);
}
