// Source file for RNSkeleton class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RNDataStructures/RNDataStructures.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

RNVessel::
RNVessel(void) :
nodes(),
spans(),
skeleton_index(-1),
type(NVESSEL_TYPES)
{
}



RNVessel::
~RNVessel(void)
{
}



RNSkeleton::
RNSkeleton(void) :
vessels()
{
   resolution[RN_X] = -1;
   resolution[RN_Y] = -1;
   resolution[RN_Z] = -1;
}



RNSkeleton::
~RNSkeleton(void)
{
}



////////////////////////////////////////////////////////////////////////
// Conversion functions
////////////////////////////////////////////////////////////////////////

R3Grid *RNSkeleton::
CreateR3Grid(int zcrop) const
{
   R3Grid *uncropped_grid = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
   if (!uncropped_grid) { return NULL; }

   // initialize to 0
   for (int iv = 0; iv < uncropped_grid->NEntries(); ++iv) {
      uncropped_grid->SetGridValue(iv, 0);
   }

   // rasterize all lines
   for (int iv = 0; iv < NVessels(); ++iv) {
      const RNVessel& vessel = Vessel(iv);
      for (int is = 0; is < vessel.NSpans(); ++is) {
         const R3Span& span = vessel.Span(is);
         uncropped_grid->RasterizeGridSpan(span.Start(), span.End(), 1.0);
      }
   }

   // post process
   for (int iv = 0; iv < uncropped_grid->NEntries(); ++iv) {
      int rasterized_value = (int)(uncropped_grid->GridValue(iv) + 0.5);
      if (rasterized_value > 0) uncropped_grid->SetGridValue(iv, 1);
      else uncropped_grid->SetGridValue(iv, 0);
   }


   // return the grid
   if (zcrop != 0) {
      R3Grid *cropped_grid = new R3Grid(resolution[RN_X], resolution[RN_Y], resolution[RN_Z] - zcrop);
      if (!cropped_grid) { return NULL; }

      for (int iz = 0; iz < cropped_grid->ZResolution(); ++iz) {
         for (int iy = 0; iy < cropped_grid->YResolution(); ++iy) {
            for (int ix = 0; ix < cropped_grid->XResolution(); ++ix) {
               cropped_grid->SetGridValue(ix, iy, iz, uncropped_grid->GridValue(ix, iy, iz + zcrop));
            }
         }
      }

      // free memory
      delete uncropped_grid;

      // return the cropped grid
      return cropped_grid;
   }

   // return the uncroppsed grid
   return uncropped_grid;
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int RNVessel::
ReadVessel(FILE *fp)
{
   // read this vessel from fp
   rn_assertion(fp != NULL);

   int skeleton_index;
   if (fscanf(fp, "VESSEL: %d\n", &skeleton_index) != 1) return 0;
   this->skeleton_index = skeleton_index;

   int type;
   if (fscanf(fp, "TYPE: %d\n", &type) != 1) return 0;
   rn_assertion((0 <= type) && (type < NVESSEL_TYPES));
   this->type = (enum VESSEL_TYPE) type;

   int nnodes;
   if (fscanf(fp, "NODES: %d\n", &nnodes) != 1) return 0;

   // read in all of the nodes
   for (int in = 0; in < nnodes; ++in) {
      int node_index;
      RNScalar xcoordinate;
      RNScalar ycoordinate;
      RNScalar zcoordinate;

      if (fscanf(fp, "%d (%lf, %lf, %lf)\n", &node_index, &xcoordinate, &ycoordinate, &zcoordinate) != 4) return 0;
      rn_assertion(node_index == in + 1);

      R3Point node = R3Point(xcoordinate, ycoordinate, zcoordinate);
      nodes.push_back(node);
      if (in > 0) {
         R3Span span = R3Span(nodes[in - 1], nodes[in]);
         spans.push_back(span);
      }
   }

   // return success
   return 1;
}



int RNVessel::
WriteVessel(FILE *fp) const
{
   // write this vessel to fp
   rn_assertion(fp != NULL);

   // write the skeleton index
   fprintf(fp, "VESSEL: %d\n", skeleton_index);

   // write the vessel type
   fprintf(fp, "TYPE: %d\n", type);

   // write the number of nodes
   fprintf(fp, "NODES: %lu\n", nodes.size());

   // write all of the nodes
   for (unsigned int in = 0; in < nodes.size(); ++in) {
      fprintf(fp, "%d (%lf, %lf, %lf)\n", in, nodes[in].X(), nodes[in].Y(), nodes[in].Z());
   }

   // return success
   return 1;
}



int RNSkeleton::
ReadSkelFile(const char skel_filename[4096])
{
   // open file
   FILE *fp = fopen(skel_filename, "r");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", skel_filename); return 0; }

   fscanf(fp, "SKELETON\n");
   if (fscanf(fp, "(%d, %d, %d)\n", &(resolution[RN_X]), &(resolution[RN_Y]), &(resolution[RN_Z])) != 3) { fprintf(stderr, "Failed to read %s\n", skel_filename); return 0; }
   int nvessels;
   if (fscanf(fp, "NVESSELS: %d\n", &nvessels) != 1) { fprintf(stderr, "Failed to read %s\n", skel_filename); return 0; }
   for (int iv = 0; iv < nvessels; ++iv) {
      RNVessel vessel = RNVessel();
      if (!vessel.ReadVessel(fp)) { fprintf(stderr, "Failed to read vessel %d\n", iv); return 0; }
      vessels.push_back(vessel);
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}



int RNSkeleton::
WriteSkelFile(const char skel_filename[4096]) const
{
   // open file
   FILE *fp = fopen(skel_filename, "w");
   if (!fp) { fprintf(stderr, "Failed to write %s\n", skel_filename); return 0; }

   fprintf(fp, "SKELETON\n");
   fprintf(fp, "(%d, %d, %d)\n", resolution[RN_X], resolution[RN_Y], resolution[RN_Z]);
   fprintf(fp, "NVESSELS: %lu\n", vessels.size());

   for (unsigned int iv = 0; iv < vessels.size(); ++iv) {
      const RNVessel& vessel = vessels[iv];
      if (!vessel.WriteVessel(fp)) { fprintf(stderr, "Failed to write vessel %d\n", iv); return 0; }
   }

   // close file
   fclose(fp);

   // return success
   return 1;
}