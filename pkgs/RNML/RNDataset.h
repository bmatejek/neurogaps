// Include file for the machine learning dataset class

#ifndef __RN_DATASET_H__
#define __RN_DATASET_H__



////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////

class RNDataset {
public:
   //////////////////////////////////
   //// CONSTRUCTORS/DESTRUCTORS ////
   //////////////////////////////////

   // constructors/destructors
   RNDataset(void);
   RNDataset(char filename[4096]);
   RNDataset(std::vector<std::vector<float> > &features, std::vector<unsigned int> &labels, std::vector<std::string> &names, std::vector<unsigned int> &identifications, unsigned int nclasses);
   ~RNDataset(void);


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   // property functions
   void CalculateEntropy(void) const;
   unsigned int NAttributes(void) const;
   unsigned int NClasses(void) const;
   unsigned int NEntries(void) const;
   int IsLabeled(void) const;


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // access functions
   unsigned int GetLabel(unsigned int index) const;
   const std::vector<unsigned int>& GetLabels(void) const;

   const std::vector<float>& GetFeature(unsigned int index) const;
   const std::vector<std::vector<float> >& GetFeatures(void) const;

   const std::string& GetName(unsigned int index) const;
   const std::vector<std::string>& GetNames(void) const;

   unsigned int GetIdentification(unsigned int index) const;
   const std::vector<unsigned int>& GetIdentifications(void) const;


   //////////////////////////
   //// OPENCV FUNCTIONS ////
   //////////////////////////

   cv::Mat *OpenCVDataMatrix(void) const;
   cv::Mat *OpenCVClassiciationsMatrix(void) const;
   cv::Mat *OpenCVAttributesMatrix(void) const;


   ///////////////////////
   //// I/O FUNCTIONS ////
   ///////////////////////
public:
   int ReadFile(char filename[4096]);
   int WriteFile(char filename[4096]) const;

private:
   int ReadASCIIFile(char ascii_filename[4096]);
   int ReadMLDBFile(char mldb_filename[4096]);
   int WriteASCIIFile(char ascii_filename[4096]) const;
   int WriteMLDBFile(char mldb_filename[4096]) const;


   //////////////////////////////
   //// CREATEDATA FUNCTIONS ////
   //////////////////////////////
public:
   void InsertDatapoint(const std::vector<float> &attributes, unsigned int label, unsigned int identification);
   void SetNames(std::vector<std::string> &names);
   void SetNClasses(unsigned int nclasses);

private:
   // instance variables
   unsigned int nclasses;
   std::vector<std::vector<float> > features;
   std::vector<unsigned int> labels;
   std::vector<std::string> names;
   std::vector<unsigned int> identifications;
};



// inline functions



inline unsigned int RNDataset::
NAttributes(void) const
{
   // return the number of attributes
   return names.size();
}



inline unsigned int RNDataset::
NClasses(void) const
{
   // return the number of classes
   return nclasses;
}



inline unsigned int RNDataset::
NEntries(void) const
{
   // return the number of entries in the dataset
   return features.size();
}



inline int RNDataset::
IsLabeled(void) const
{
   // return if the dataset is labeled
   return labels.size() == features.size();
}



inline unsigned int RNDataset::
GetLabel(unsigned int index) const
{
   // return the label for this index
   return labels[index];
}



inline const std::vector<unsigned int>& RNDataset::
GetLabels(void) const
{
   // return labels
   return labels;
}



inline const std::vector<float>& RNDataset::
GetFeature(unsigned int index) const
{
   // return the feature for this index
   return features[index];
}



inline const std::vector<std::vector<float> >& RNDataset::
GetFeatures(void) const
{
   // return the features vector
   return features;
}



inline const std::string& RNDataset::
GetName(unsigned int index) const
{
   // return the name of this feature
   return names[index];
}



inline const std::vector<std::string>& RNDataset::
GetNames(void) const
{
   // return the names vector
   return names;
}



inline unsigned int RNDataset::
GetIdentification(unsigned int index) const
{
   // return the number of identifying features
   return identifications[index];
}



inline const std::vector<unsigned int>& RNDataset::
GetIdentifications(void) const
{
   return identifications;
}



#endif
