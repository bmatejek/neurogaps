// Include file for the machine learning random forest classifier class


#ifndef __RN_RANDOM_FOREST_H__
#define __RN_RANDOM_FOREST_H__


//////////////////////////////////////////////////////////////////////
//// CLASS DEFINITION
//////////////////////////////////////////////////////////////////////

class RNRandomForest : public RNClassifier {
public:

   //////////////////////////////////
   //// CONSTRUCTORS/DESTRUCTORS ////
   //////////////////////////////////

   // constructor/destructors
   RNRandomForest(void);
   RNRandomForest(RNDataset &data, int max_depth, int min_sample_count, RNBoolean rebalance, int nactive_vars, int ntrees, float forest_accuracy);
   virtual ~RNRandomForest(void);


   //////////////////////////////////
   //// CLASSIFICATION FUNCTIONS ////
   //////////////////////////////////

   virtual unsigned int Classify(const std::vector<float> &attributes) const;
   virtual float PredictBinaryProbability(const std::vector<float> &attributes) const;

   virtual std::vector<unsigned int> Classify(RNDataset& testing_data) const;
   virtual std::vector<float> PredictBinaryProbabilities(RNDataset &testing_data) const;


   ////////////////////////////////
   //// INPUT/OUTPUT FUNCTIONS ////
   ////////////////////////////////

   virtual unsigned int **PrintPerformance(RNDataset& testing_dataset) const;
   virtual void LoadClassifier(const char input_filename[4096]);
   virtual void SaveClassifier(const char output_filename[4096]) const;


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////
public:
   // class type definitions
   RN_CLASS_TYPE_DECLARATIONS(RNRandomForest);


private:
   char classifier_type_name[128];
   CvRTrees *random_forest;
};

#endif
