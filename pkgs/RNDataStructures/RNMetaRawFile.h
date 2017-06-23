// Include file for the meta raw file

#ifndef __RN_META_RAW_FILE_H__
#define __RN_META_RAW_FILE_H__


/////////////////////////////////////////////////////////////////////
// META DEFINITION
/////////////////////////////////////////////////////////////////////

class RNMeta {
public:
   // constructors/destructors
   RNMeta(void);
   RNMeta(const char *input_data_type, int xres, int yres, int zres, int ngrids);
   ~RNMeta(void);

   // access functions
   const char *DataType(void) const;
   int XResolution(void) const;
   int YResolution(void) const;
   int ZResolution(void) const;
   int NGrids(void) const;

private:
   char data_type[4096];
   int xres; 
   int yres;
   int zres;
   int ngrids;
};



/////////////////////////////////////////////////////////////////////
// META/RAW INPUT/OUTPUT DEFINITION
/////////////////////////////////////////////////////////////////////

// input functions
R3Grid *RNReadNeuronMetaRawFile(const char *root_filename);
R3Grid **RNReadNeuronMetaRawFile(const char *filename, RNBoolean multidimensional);
RNMeta RNReadNeuronMetaFile(const char *meta_filename);
int RNReadNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grid);
int RNReadNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grids[3]);

// output functions
int RNWriteNeuronMetaRawFile(const char *filename, const RNMeta& meta, R3Grid *grid);
int RNWriteNeuronMetaRawFile(const char *filename, const RNMeta& meta, R3Grid *grids[3]);
int RNWriteNeuronMetaFile(const char *meta_filename, const RNMeta& meta);
int RNWriteNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grid);
int RNWriteNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grids[3]);

#endif