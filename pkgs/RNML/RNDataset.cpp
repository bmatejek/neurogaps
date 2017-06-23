// Source file for the machine learning dataset class



///////////////////////////////////////////////////////////////
//// Include files
///////////////////////////////////////////////////////////////

#include "RNML.h"



///////////////////////////////////////////////////////////////
// Constructor/destructor functions 
///////////////////////////////////////////////////////////////

RNDataset::
RNDataset(void) :
features(),
labels()
{
}



RNDataset::
RNDataset(char filename[4096]) :
features(),
labels(),
names(),
identifications()
{
   ReadFile(filename);
}



RNDataset::
RNDataset(std::vector<std::vector<float> > &features, std::vector<unsigned int> &labels, std::vector<std::string> &names, std::vector<unsigned int> &identifications, unsigned int nclasses) :
nclasses(nclasses),
features(features),
labels(labels),
names(names),
identifications(identifications)
{
}



RNDataset::
~RNDataset(void)
{
}



///////////////////////////////////////////////////////////////
// OpenCV functions 
///////////////////////////////////////////////////////////////

cv::Mat *RNDataset::
OpenCVDataMatrix(void) const
{
   // create data matrix
   cv::Mat *data = new cv::Mat(NEntries(), NAttributes(), CV_32FC1);

   for (unsigned int in = 0; in < NEntries(); ++in) {
      for (unsigned int ia = 0; ia < NAttributes(); ++ia) {
         data->at<float>(in, ia) = features[in][ia];
      }
   }

   // return data matrix
   return data;
}



cv::Mat *RNDataset::
OpenCVClassiciationsMatrix(void) const
{
   // create classifications matrix
   cv::Mat *classifications = new cv::Mat(NEntries(), 1, CV_32FC1);

   for (unsigned int in = 0; in < NEntries(); ++in) {
      classifications->at<float>(in, 0) = labels[in];
   }

   // return classifications matrix
   return classifications;
}



cv::Mat *RNDataset::
OpenCVAttributesMatrix(void) const
{
   // create the attributes matrix (all are numerical)
   cv::Mat *attributes = new cv::Mat(NAttributes() + 1, 1, CV_8U);
   attributes->setTo(cv::Scalar(CV_VAR_NUMERICAL));

   // it is a classification problem so there are a discrete number of classes
   attributes->at<uchar>(NAttributes(), 0) = CV_VAR_CATEGORICAL;

   // return the attributes matrix
   return attributes;
}



///////////////////////////////////////////////////////////////
// Property functions
///////////////////////////////////////////////////////////////

// create entropy struct
struct Attribute {
   Attribute(float attribute_value, int label) :
      attribute_value(attribute_value),
      label(label)
   {}

   float attribute_value;
   int label;
};



int AttributeSort(Attribute one, Attribute two)
{
   return one.attribute_value < two.attribute_value;
}



static RNScalar Entropy(int npositives, int nnegatives)
{
   RNScalar positive_proportion = npositives / (RNScalar)(npositives + nnegatives);
   RNScalar negative_proportion = nnegatives / (RNScalar)(npositives + nnegatives);

   // return the entropy
   return -1 * positive_proportion * log2(positive_proportion) - negative_proportion * log2(negative_proportion);
}



void RNDataset::
CalculateEntropy(void) const
{
   // make sure the data is labeled
   rn_assertion(IsLabeled());

   // create array for entropies
   RNScalar *information_gains = new RNScalar[names.size()];
   RNScalar *proportion_correct = new RNScalar[names.size()];
   RNScalar *split_values = new RNScalar[names.size()];

   // get number of positive and negative results
   int nparent_positives = 0;
   int nparent_negatives = 0;
   for (unsigned int ie = 0; ie < NEntries(); ++ie) {
      if (labels[ie]) nparent_positives++;
      else nparent_negatives++;
   }

   // get parent entropy
   RNScalar parent_entropy = Entropy(nparent_positives, nparent_negatives);

   // go through every attribute
   for (unsigned int ia = 0; ia < names.size(); ++ia) {
      // create vector for all attribue values for this feature
      std::vector<Attribute> attributes = std::vector<Attribute>();

      // go through every entry in the set
      for (unsigned int ie = 0; ie < NEntries(); ++ie) {
         attributes.push_back(Attribute(features[ie][ia], labels[ie]));
      }

      // sort all of the attributes
      sort(attributes.begin(), attributes.end(), AttributeSort);

      // find the number of right and left ones and zeros
      int *left_negatives = new int[attributes.size()];
      int *right_positives = new int[attributes.size()];
      for (int ie = 0; ie < (int)attributes.size(); ++ie) {
         left_negatives[ie] = 0;
         right_positives[ie] = 0;
      }

      // go from from left to right to get left negatives array
      for (int ie = 1; ie < (int)attributes.size(); ++ie) {
         if (attributes[ie - 1].label) left_negatives[ie] = left_negatives[ie - 1];
         else left_negatives[ie] = left_negatives[ie - 1] + 1;
      }

      // go from right to left to get right positives array
      for (int ie = (int)attributes.size() - 2; ie >= 0; --ie) {
         if (attributes[ie + 1].label) right_positives[ie] = right_positives[ie + 1] + 1;
         else right_positives[ie] = right_positives[ie + 1];
      }

      RNScalar best_information_gain = 0.0;

      // calculate the entropy at each location
      for (int ie = 0; ie < (int)attributes.size(); ++ie) {
         int right_positive = right_positives[ie];
         int left_negative = left_negatives[ie];
         int left_positive = (ie - left_negative);
         int right_negative = ((attributes.size() - 1) - ie - right_positive);

         // find current label
         int current_label = labels[ie];

         // give benefit of doubt for splitting here
         if (current_label) {
            if (left_positive > right_positive) left_positive++;
            else right_positive++;
         }
         else {
            if (left_negative > right_negative) left_negative++;
            else right_negative++;
         }

         // calculate the left and right entropy
         RNScalar left_entropy = Entropy(left_positive, left_negative);
         RNScalar right_entropy = Entropy(right_positive, right_negative);

         RNScalar left_proportion = (left_positive + left_negative) / (RNScalar)(attributes.size());
         RNScalar right_proportion = (right_positive + right_negative) / (RNScalar)(attributes.size());

         RNScalar average_child_entropy = left_proportion * left_entropy + right_proportion * right_entropy;
         RNScalar information_gain = parent_entropy - average_child_entropy;

         // used for discrete variables
         if ((ie != (int)attributes.size() - 1) && (attributes[ie].attribute_value == attributes[ie + 1].attribute_value)) continue;

         if (information_gain > best_information_gain) {
            best_information_gain = information_gain;

            RNScalar split_right = (right_negative + left_positive) / (RNScalar)attributes.size();
            RNScalar split_left = (left_negative + right_positive) / (RNScalar)attributes.size();

            if (split_right > split_left)
               proportion_correct[ia] = split_right;
            else
               proportion_correct[ia] = split_left;

            split_values[ia] = attributes[ie].attribute_value;
         }
      }

      information_gains[ia] = best_information_gain;
   }

   // print out header
   printf("        Attribute Name        |   Information Gain   |  Proportion Correct  |   Split Values   \n");
   printf("------------------------------+----------------------+----------------------+------------------\n");
   // print out all results
   for (unsigned int in = 0; in < names.size(); ++in) {
      char attribute_name[4096];
      strncpy(attribute_name, names[in].c_str(), 4096);
      char *extp = strrchr(attribute_name, '\n');
      if (extp) *extp = '\0';

      printf("%-30s&%19.6lf   &%20.6lf  &%15.6lf\\\\\n", attribute_name, information_gains[in], proportion_correct[in], split_values[in]);
   }

   // free memory
   delete[] information_gains;
   delete[] proportion_correct;
   delete[] split_values;
}




///////////////////////////////////////////////////////////////
// Input/output functions
///////////////////////////////////////////////////////////////

int RNDataset::
ReadFile(char filename[4096])
{
   // get file ending
   char *extp = strrchr(filename, '.');
   if (!strncmp(extp, ".txt", 4)) {
      if (!ReadASCIIFile(filename)) return 0;
   }
   else if (!strncmp(extp, ".mldb", 5)) {
      if (!ReadMLDBFile(filename)) return 0;
   }
   else { fprintf(stderr, "Failed to recognize extension %s\n", extp); return 0; }

   // return success
   return 1;
}



int RNDataset::
WriteFile(char filename[4096]) const
{
   // get file ending
   char *extp = strrchr(filename, '.');
   if (!strncmp(extp, ".txt", 4)) {
      if (!WriteASCIIFile(filename)) return 0;
   }
   else if (!strncmp(extp, ".mldb", 5)) {
      if (!WriteMLDBFile(filename)) return 0;
   }
   else { fprintf(stderr, "Failed to recognize extension %s\n", extp); return 0; }

   // return success
   return 1;
}



int RNDataset::
ReadASCIIFile(char ascii_filename[4096])
{
   // clear all instance variables
   features.clear();
   labels.clear();
   names.clear();
   identifications.clear();

   // open file
   FILE *fp = fopen(ascii_filename, "r");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", ascii_filename); return 0; }

   // read in the number of data points and features
   int ndata_points;
   int nfeatures;
   int is_labeled;
   fscanf(fp, "%d %d %d %d\n", &ndata_points, &nfeatures, &nclasses, &is_labeled);

   // read in all of the names
   for (int in = 0; in < nfeatures; ++in) {
      char feature_name[128];
      fscanf(fp, "%s\n", &(feature_name[0]));
      names.push_back(std::string(feature_name));
   }

   // read in all of the features for each data point
   for (int id = 0; id < ndata_points; ++id) {
      // create new feature vector
      features.push_back(std::vector<float>());
      for (int ia = 0; ia < nfeatures; ++ia) {
         float feature;
         fscanf(fp, "%f ", &feature);
         features[id].push_back(feature);
      }
      if (is_labeled) {
         unsigned int label;
         fscanf(fp, "%u ", &label);
         labels.push_back(label);
      }
      unsigned int identification;
      fscanf(fp, "%u\n", &identification);
      identifications.push_back(identification);
   }

   // close file 
   fclose(fp);

   // return success
   return 1;
}



int RNDataset::
ReadMLDBFile(char mldb_filename[4096])
{
   // clear all instance variables
   features.clear();
   labels.clear();
   names.clear();
   identifications.clear();

   // open file
   FILE *fp = fopen(mldb_filename, "rb");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", mldb_filename); return 0; }

   int ndata_points;
   int nfeatures;
   RNBoolean is_labeled;
   fread(&ndata_points, sizeof(int), 1, fp);
   fread(&nfeatures, sizeof(int), 1, fp);
   fread(&nclasses, sizeof(int), 1, fp);
   fread(&is_labeled, sizeof(RNBoolean), 1, fp);

   // read in the feature names
   for (int ia = 0; ia < nfeatures; ++ia) {
      char feature_name[128];
      fread(&(feature_name[0]), sizeof(char), 128, fp);
      names.push_back(std::string(feature_name));
   }

   // read in the data points
   printf("Reading data points...\n  "); fflush(stdout);
   for (int in = 0; in < ndata_points; ++in) {
      RNProgressBar(in, ndata_points);
      features.push_back(std::vector<float>());

      // read attributes
      float *attribute_values = new float[nfeatures];
      fread(attribute_values, sizeof(float), nfeatures, fp);
      for (int ia = 0; ia < nfeatures; ++ia) {
         features[in].push_back(attribute_values[ia]);
      }
      delete[] attribute_values;

      if (is_labeled) {
         unsigned int label;
         fread(&label, sizeof(unsigned int), 1, fp);
         labels.push_back(label);
      }
      unsigned int identification;
      fread(&identification, sizeof(unsigned int), 1, fp);
      identifications.push_back(identification);
   }
   printf("\ndone.\n");

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNDataset::
WriteASCIIFile(char ascii_filename[4096]) const
{
   // open file
   FILE *fp = fopen(ascii_filename, "w");
   if (!fp) { fprintf(stderr, "Failed to write %s\n", ascii_filename); return 0; }

   // print out the number of data points and features
   int ndata_points = features.size();
   rn_assertion(ndata_points != 0);
   int nfeatures = features[0].size();
   fprintf(fp, "%d %d %d %d\n", ndata_points, nfeatures, nclasses, IsLabeled());

   // print out the names of all of the features
   for (int in = 0; in < nfeatures; ++in) {
      fprintf(fp, "%s\n", names[in].c_str());
   }

   // print out all of the features for each data point
   for (int id = 0; id < ndata_points; ++id) {
      for (int ia = 0; ia < nfeatures; ++ia) {
         fprintf(fp, "%f ", features[id][ia]);
      }
      if (IsLabeled()) {
         fprintf(fp, "%u ", labels[id]);
      }
      fprintf(fp, "%u\n", identifications[id]);
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNDataset::
WriteMLDBFile(char mldb_filename[4096]) const
{
   // open file
   FILE *fp = fopen(mldb_filename, "wb");
   if (!fp) { fprintf(stderr, "Failed to write %s\n", mldb_filename); return 0; }

   int ndata_points = NEntries();
   int nfeatures = NAttributes();
   RNBoolean is_labeled = IsLabeled();
   fwrite(&ndata_points, sizeof(int), 1, fp);
   fwrite(&nfeatures, sizeof(int), 1, fp);
   fwrite(&nclasses, sizeof(int), 1, fp);
   fwrite(&is_labeled, sizeof(RNBoolean), 1, fp);

   // read in the feature names
   for (int ia = 0; ia < nfeatures; ++ia) {
      char feature_name[128];
      sprintf(feature_name, "%s\n", names[ia].c_str());
      fwrite(&(feature_name[0]), sizeof(char), 128, fp);
   }

   // read in the data points
   for (int in = 0; in < ndata_points; ++in) {
      for (int ia = 0; ia < nfeatures; ++ia) {
         float attribute_value = features[in][ia];
         fwrite(&attribute_value, sizeof(float), 1, fp);
      }
      if (is_labeled) {
         unsigned int label = labels[in];
         fwrite(&label, sizeof(unsigned int), 1, fp);
      }
      unsigned int identification = identifications[in];
      fwrite(&identification, sizeof(unsigned int), 1, fp);
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



///////////////////////////////////////////////////////////////
// createdata functions
///////////////////////////////////////////////////////////////

void RNDataset::
InsertDatapoint(const std::vector<float> &attributes, unsigned int label, unsigned int identification)
{
   features.push_back(attributes);
   labels.push_back(label);
   identifications.push_back(identification);
}



void RNDataset::
SetNames(std::vector<std::string> &names)
{
   this->names = std::vector<std::string>(names);
}



void RNDataset::
SetNClasses(unsigned int nclasses)
{
   this->nclasses = nclasses;
}


