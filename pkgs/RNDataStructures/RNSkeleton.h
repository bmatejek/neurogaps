// Include file for the RNSkeleton class

#ifndef __RN_SKELETON_H__
#define __RN_SKELETON_H__

#include <vector>



/* useful constants */

enum VESSEL_TYPE {
   ARTERY = 1,
   ARTERIOLE = 2,
   VEIN = 3,
   VENULE = 4,
   CAPILLARY = 5,
   NVESSEL_TYPES = 6
};



static const char *vessel_types[NVESSEL_TYPES] = {
   "", "Artery", "Arteriole", "Vein", "Venule", "Capillary"
};



/////////////////////////////////////////////////////////////////////
// VESSEL CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNVessel {
public:
   // constructors/destructors
   RNVessel(void);
   ~RNVessel(void);

   // access functions
   int NNodes(void) const;
   const R3Point& Node(int node_index) const;
   int NSpans(void) const;
   const R3Span& Span(int span_index) const;
   int SkeletonIndex(void) const;
   const char *VesselType(void) const;

protected:
   // I/O functions
   int ReadVessel(FILE *fp);
   int WriteVessel(FILE *fp) const;

private:
   // instance variables
   friend class RNSkeleton;
   std::vector<R3Point> nodes;
   std::vector<R3Span> spans;
   int skeleton_index;
   enum VESSEL_TYPE type;
};



/* inline vessel functions */

inline int RNVessel::
NNodes(void) const
{
   // return the number of nodes in this vessel
   return nodes.size();
}



inline const R3Point& RNVessel::
Node(int node_index) const
{
   rn_assertion((0 <= node_index) && (node_index < (int)nodes.size()));
   // return the vessel 
   return nodes[node_index];
}



inline int RNVessel::
NSpans(void) const
{
   // return the number of lines in this vessel
   return spans.size();
}



inline const R3Span& RNVessel::
Span(int span_index) const
{
   rn_assertion((0 <= span_index) && (span_index < (int)spans.size()));
   // return the line
   return spans[span_index];
}



inline int RNVessel::
SkeletonIndex(void) const
{
   // return index within skeleton
   return skeleton_index;
}



inline const char *RNVessel::
VesselType(void) const
{
   // return the name of this vessel
   return vessel_types[type];
}



/////////////////////////////////////////////////////////////////////
// SKELETON CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNSkeleton {
public:
   // consructors/destructors
   RNSkeleton(void);
   ~RNSkeleton(void);

   // access functions
   int NVessels(void) const;
   const RNVessel& Vessel(int vessel_index) const;

   // conversion functions
   R3Grid *CreateR3Grid(int zcrop = 0) const;

   // I/O functions
   int ReadSkelFile(const char skel_filename[4096]);
   int WriteSkelFile(const char skel_filename[4096]) const;

private:
   // instance variables
   std::vector<RNVessel> vessels;
   int resolution[3];
};



/* inline functions */

inline int RNSkeleton::
NVessels(void) const
{
   // return the number of vessels
   return vessels.size();
}



inline const RNVessel& RNSkeleton::
Vessel(int vessel_index) const
{
   rn_assertion((0 <= vessel_index) && (vessel_index < (int)vessels.size()));
   // return this vessel
   return vessels[vessel_index];
}



#endif