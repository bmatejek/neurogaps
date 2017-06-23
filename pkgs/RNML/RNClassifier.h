// Include file for the machine learning classifier class


#ifndef __RN_CLASSIFIER_H__
#define __RN_CLASSIFIER_H__


//////////////////////////////////////////////////////////////////////
//// CLASS DEFINITION
//////////////////////////////////////////////////////////////////////

class RNClassifier {
  public:

   /////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTORS ////
   /////////////////////////////////
   
   // constructors/destructors
   RNClassifier(void);
   virtual ~RNClassifier(void);
   
   
   //////////////////////////////////
   //// CLASSIFICATION FUNCTIONS ////
   //////////////////////////////////

   // save the classifier
   virtual unsigned int Classify(const std::vector<float>& attributes) const;
   virtual float PredictBinaryProbability(const std::vector<float>& attributes) const;
   
   virtual std::vector<unsigned int> Classify(RNDataset& testing_data) const;
   virtual std::vector<float> PredictBinaryProbabilities(RNDataset& testing_data) const;

   
   ////////////////////////////////
   //// INPUT/OUTPUT FUNCTIONS ////
   ////////////////////////////////

   // performance characteristics
   virtual unsigned int **PrintPerformance(RNDataset& dataset) const;
   virtual void LoadClassifier(const char input_filename[4096]);
   virtual void SaveClassifier(const char output_filename[4096]) const;


public:
   // class type declarations   
   RN_CLASS_TYPE_DECLARATIONS(RNClassifier);
};

#endif
