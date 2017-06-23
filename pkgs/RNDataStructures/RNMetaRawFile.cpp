// Source file for meta raw file 



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RNDataStructures/RNDataStructures.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

RNMeta::
RNMeta(void) :
xres(-1),
yres(-1),
zres(-1),
ngrids(-1)
{
   data_type[0] = '\0';
}



RNMeta::
RNMeta(const char *input_data_type, int xres, int yres, int zres, int ngrids) :
xres(xres),
yres(yres),
zres(zres),
ngrids(ngrids)
{
   strncpy(data_type, input_data_type, 4096);
}



RNMeta::
~RNMeta(void)
{
}



const char *RNMeta::
DataType(void) const
{
   // return data type
   return data_type;
}



int RNMeta::
XResolution(void) const
{
   // return xresolution
   return xres;
}



int RNMeta::
YResolution(void) const
{
   // return y resolution
   return yres;
}



int RNMeta::
ZResolution(void) const
{
   // return z resolution
   return zres;
}



int RNMeta::
NGrids(void) const
{
   // return number of grids
   return ngrids;
}



////////////////////////////////////////////////////////////////////////
// Int/RNScalar input/output functions
////////////////////////////////////////////////////////////////////////

static int
RNReadNeuronRawValue(FILE *fp, int format, RNScalar &dest)
{
   // convert file pointer to RNScalar
   switch (format) {
   case 0: {
      // Uint8
      unsigned char value;
      fread(&value, sizeof(unsigned char), 1, fp);
      dest = (RNScalar)value;
      break; }
   case 1: {
      // Uint16
      unsigned short value;
      fread(&value, sizeof(unsigned short), 1, fp);
      dest = (RNScalar)value;
      break; }
   case 2: {
      // Uint32
      unsigned int value;
      fread(&value, sizeof(unsigned int), 1, fp);
      dest = (RNScalar)value;
      break; }
   case 3: {
      // Int32
      int value;
      fread(&value, sizeof(int), 1, fp);
      dest = (RNScalar)value;
      break; }
   case 4: {
      // Float32
      float value;
      fread(&value, sizeof(float), 1, fp);
      dest = (RNScalar)value;
      break; }
   case 5: {
      // RNScalar
      fread(&dest, sizeof(RNScalar), 1, fp);
      break; }
   default: {
      // unrecognized format
      fprintf(stderr, "Unrecognized format %d\n", format);
      return 0; }
   }

   // return success
   return 1;
}



static int
RNWriteNeuronRawValue(FILE *fp, int format, RNScalar source)
{
   // write int to file pointer
   switch (format) {
   case 0: {
      // Uint8
      unsigned char value = (unsigned char)(source + 0.5);
      fwrite(&value, sizeof(unsigned char), 1, fp);
      break; }
   case 1: {
      // Uint16
      unsigned short value = (unsigned short)(source + 0.5);
      fwrite(&value, sizeof(unsigned short), 1, fp);
      break; }
   case 2: {
      // Uint32
      unsigned int value = (unsigned int)(source + 0.5);
      fwrite(&value, sizeof(unsigned int), 1, fp);
      break; }
   case 3: {
      // Int32
      int value = (int)(source + 0.5);
      fwrite(&value, sizeof(int), 1, fp);
      break; }
   case 4: {
      // Float32
      float value = (float)source;
      fwrite(&value, sizeof(float), 1, fp);
      break; }
   case 5: {
      // RNScalar
      fwrite(&source, sizeof(RNScalar), 1, fp);
      break; }
   default: {
      // unrecognized format
      fprintf(stderr, "Unrecognized format %d\n", format);
      return 0; }
   }

   // return success
   return 1;
}



/////////////////////////////////////////////////////////////////////
// Meta input/output definition
/////////////////////////////////////////////////////////////////////

R3Grid *RNReadNeuronMetaRawFile(const char *filename)
{
   // make root filename if it isn't already
   char root_filename[4096];
   strncpy(root_filename, filename, 4096);
   char *extp = strrchr(root_filename, '.');
   if (extp) *extp = '\0';

   // get meta filename
   char meta_filename[4096];
   sprintf(meta_filename, "%s.meta", root_filename);

   // get raw filename
   char raw_filename[4096];
   sprintf(raw_filename, "%s.raw", root_filename);

   // read meta file
   RNMeta meta = RNReadNeuronMetaFile(meta_filename);
   rn_assertion((meta.XResolution() > 0) && (meta.YResolution() > 0) && (meta.ZResolution() > 0) && (meta.NGrids() == 1));

   // create grids 
   R3Grid *grid = new R3Grid(meta.XResolution(), meta.YResolution(), meta.ZResolution());

   // read raw file
   if (!RNReadNeuronRawFile(raw_filename, meta, grid)) { fprintf(stderr, "Failed to read raw file\n"); return NULL; }

   // return success
   return grid;
}



R3Grid **RNReadNeuronMetaRawFile(const char *filename, RNBoolean multidimensional)
{
   // make root filename if it isn't already
   char root_filename[4096];
   strncpy(root_filename, filename, 4096);
   char *extp = strrchr(root_filename, '.');
   if (extp) *extp = '\0';

   // get meta filename
   char meta_filename[4096];
   sprintf(meta_filename, "%s.meta", root_filename);

   // get raw filename
   char raw_filename[4096];
   sprintf(raw_filename, "%s.raw", root_filename);

   // read meta file
   RNMeta meta = RNReadNeuronMetaFile(meta_filename);
   rn_assertion((meta.XResolution() > 0) && (meta.YResolution() > 0) && (meta.ZResolution() > 0) && (meta.NGrids() == 3));

   // create grids 
   R3Grid **grids = new R3Grid *[meta.NGrids()];
   for (int i = 0; i < meta.NGrids(); ++i) {
      grids[i] = new R3Grid(meta.XResolution(), meta.YResolution(), meta.ZResolution());
   }

   // read raw file
   if (!RNReadNeuronRawFile(raw_filename, meta, grids)) { fprintf(stderr, "Failed to read raw file\n"); return NULL; }

   // return success
   return grids;
}



RNMeta RNReadNeuronMetaFile(const char *meta_filename) {
   // open file
   FILE *fp = fopen(meta_filename, "r");
   if (!fp) { fprintf(stderr, "Unable to open meta file: %s\n", meta_filename); return RNMeta(); }
   // read meta file
   char buffer[4096];
   if (!fgets(buffer, 4096, fp)) { fprintf(stderr, "Unable to read meta file %s\n", meta_filename); return RNMeta(); }

   // parse data type
   char *data_type = strtok(buffer, "(");
   if (!data_type) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return RNMeta(); }

   // parse dimensions
   char *bufferp = strtok(NULL, ",");
   if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return RNMeta(); }
   int xres = atoi(bufferp);
   bufferp = strtok(NULL, ",");
   if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return RNMeta(); }
   int yres = atoi(bufferp);
   bufferp = strtok(NULL, ",)");
   if (!bufferp) { fprintf(stderr, "Error parsing meta file %s\n", meta_filename); return RNMeta(); }
   int zres = atoi(bufferp);
   bufferp = strtok(NULL, ")");
   int ngrids = (bufferp) ? atoi(bufferp) : 1;

   // close file
   fclose(fp);

   // return meta data type
   RNMeta meta = RNMeta(data_type, xres, yres, zres, ngrids);

   return meta;
}



int RNReadNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grid)
{
   rn_assertion((grid != NULL) && (meta.NGrids() == 1));

   // open file
   FILE *fp = fopen(raw_filename, "rb");
   if (!fp) { fprintf(stderr, "Unable to open raw file %s\n", raw_filename); return 0; }

   // determine data format
   int format = -1;
   if (!strcmp(meta.DataType(), "Uint8")) format = 0;
   else if (!strcmp(meta.DataType(), "Uint16")) format = 1;
   else if (!strcmp(meta.DataType(), "Uint32")) format = 2;
   else if (!strcmp(meta.DataType(), "Int32")) format = 3;
   else if (!strcmp(meta.DataType(), "Float32")) format = 4;
   else if (!strcmp(meta.DataType(), "RNScalar")) format = 5;

   // read file
   for (int iz = 0; iz < meta.ZResolution(); ++iz) {
      for (int iy = 0; iy < meta.YResolution(); ++iy) {
         for (int ix = 0; ix < meta.XResolution(); ++ix) {
            RNScalar value = 0;
            if (!RNReadNeuronRawValue(fp, format, value)) return 0;
            grid->SetGridValue(ix, iy, iz, value);
         }
      }
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNReadNeuronRawFile(const char *raw_filename, const RNMeta& meta, R3Grid *grids[3])
{
   rn_assertion((meta.NGrids() == 3) && (grids[RN_X] != NULL) && (grids[RN_Y] != NULL) && (grids[RN_Z] != NULL));

   // open file
   FILE *fp = fopen(raw_filename, "rb");
   if (!fp) { fprintf(stderr, "Unable to open raw file %s\n", raw_filename); return 0; }

   // determine data format
   int format = -1;
   if (!strcmp(meta.DataType(), "Uint8")) format = 0;
   else if (!strcmp(meta.DataType(), "Uint16")) format = 1;
   else if (!strcmp(meta.DataType(), "Uint32")) format = 2;
   else if (!strcmp(meta.DataType(), "Int32")) format = 3;
   else if (!strcmp(meta.DataType(), "Float32")) format = 4;
   else if (!strcmp(meta.DataType(), "RNScalar")) format = 5;

   // read file
   for (int ig = 0; ig < meta.NGrids(); ++ig) {
      rn_assertion(grids[ig]);
      for (int iz = 0; iz < meta.ZResolution(); ++iz) {
         for (int iy = 0; iy < meta.YResolution(); ++iy) {
            for (int ix = 0; ix < meta.XResolution(); ++ix) {
               RNScalar value = 0;
               if (!RNReadNeuronRawValue(fp, format, value)) return 0;
               grids[ig]->SetGridValue(ix, iy, iz, value);
            }
         }
      }
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNWriteNeuronMetaRawFile(const char *filename, const RNMeta& meta, R3Grid *grid)
{
   rn_assertion((meta.NGrids() == 1) && (grid != NULL));

   // write meta file
   char meta_filename[4096];
   sprintf(meta_filename, "%s.meta", filename);
   if (!RNWriteNeuronMetaFile(meta_filename, meta)) return 0;

   // write raw file
   char raw_filename[4096];
   sprintf(raw_filename, "%s.raw", filename);
   if (!RNWriteNeuronRawFile(raw_filename, meta, grid));

   // return success
   return 1;
}



int RNWriteNeuronMetaRawFile(const char *filename, const RNMeta& meta, R3Grid *grids[3])
{
   rn_assertion((meta.NGrids() == 3) && (grids[RN_X] != NULL) && (grids[RN_Y] != NULL) && (grids[RN_Z] != NULL));

   // write meta file
   char meta_filename[4096];
   sprintf(meta_filename, "%s.meta", filename);
   if (!RNWriteNeuronMetaFile(meta_filename, meta)) return 0;

   // write raw file
   char raw_filename[4096];
   sprintf(raw_filename, "%s.raw", filename);
   if (!RNWriteNeuronRawFile(raw_filename, meta, grids));

   // return success
   return 1;
}



int RNWriteNeuronMetaFile(const char *meta_filename, const RNMeta& meta)
{
   // open meta file
   FILE *fp = fopen(meta_filename, "w");
   if (!fp) { fprintf(stderr, "Failed to open %s to write\n", meta_filename); return 0; }

   if (meta.NGrids() == 1)
      fprintf(fp, "%s(%d,%d,%d)", meta.DataType(), meta.XResolution(), meta.YResolution(), meta.ZResolution());
   else
      fprintf(fp, "%s(%d,%d,%d,%d)", meta.DataType(), meta.XResolution(), meta.YResolution(), meta.ZResolution(), meta.NGrids());

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNWriteNeuronRawFile(const char *raw_fileame, const RNMeta& meta, R3Grid *grid)
{
   rn_assertion((meta.NGrids() == 1) && (grid != NULL));
   // open file
   FILE *fp = fopen(raw_fileame, "wb");
   if (!fp) { fprintf(stderr, "Failed to open %s to write\n", raw_fileame); return 0; }

   // determine data format
   int format = -1;
   if (!strcmp(meta.DataType(), "Uint8")) format = 0;
   else if (!strcmp(meta.DataType(), "Uint16")) format = 1;
   else if (!strcmp(meta.DataType(), "Uint32")) format = 2;
   else if (!strcmp(meta.DataType(), "Int32")) format = 3;
   else if (!strcmp(meta.DataType(), "Float32")) format = 4;
   else if (!strcmp(meta.DataType(), "RNScalar")) format = 5;

   for (int iz = 0; iz < meta.ZResolution(); ++iz) {
      for (int iy = 0; iy < meta.YResolution(); ++iy) {
         for (int ix = 0; ix < meta.XResolution(); ++ix) {
            RNScalar value = grid->GridValue(ix, iy, iz);
            if (!RNWriteNeuronRawValue(fp, format, value)) return 0;
         }
      }
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNWriteNeuronRawFile(const char *raw_fileame, const RNMeta& meta, R3Grid *grids[3])
{
   rn_assertion((meta.NGrids() == 3) && (grids[RN_X] != NULL) && (grids[RN_Y] != NULL) && (grids[RN_Z] != NULL));

   // open file
   FILE *fp = fopen(raw_fileame, "wb");
   if (!fp) { fprintf(stderr, "Failed to open %s to write\n", raw_fileame); return 0; }

   // determine data format
   int format = -1;
   if (!strcmp(meta.DataType(), "Uint8")) format = 0;
   else if (!strcmp(meta.DataType(), "Uint16")) format = 1;
   else if (!strcmp(meta.DataType(), "Uint32")) format = 2;
   else if (!strcmp(meta.DataType(), "Int32")) format = 3;
   else if (!strcmp(meta.DataType(), "Float32")) format = 4;
   else if (!strcmp(meta.DataType(), "RNScalar")) format = 5;

   for (int dim = 0; dim < meta.NGrids(); ++dim) {
      rn_assertion(grids[dim] != NULL);
      for (int iz = 0; iz < meta.ZResolution(); ++iz) {
         for (int iy = 0; iy < meta.YResolution(); ++iy) {
            for (int ix = 0; ix < meta.XResolution(); ++ix) {
               RNScalar value = grids[dim]->GridValue(ix, iy, iz);
               if (!RNWriteNeuronRawValue(fp, format, value)) return 0;
            }
         }
      }
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}
