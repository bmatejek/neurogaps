// Include file for the neuron cellular class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronCellular : public NeuronSupervoxel {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronCellular(void);
   virtual ~NeuronCellular(void);


   //////////////////////////
   //// ACCESS FUNCTIONS ////
   //////////////////////////

   // supervoxel access functions
   virtual int SupervoxelIndex(void) const;


   ////////////////////////////
   //// PROPERTY FUNCTIONS ////
   ////////////////////////////

   virtual RNBoolean IsCellular(void) const;
   virtual RNBoolean IsExtracellular(void) const;


   ///////////////////////////
   //// DISPLAY FUNCTIONS ////
   ///////////////////////////
public:
   // draw functions
   virtual void Draw(void) const;

   // print functions
   virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;


   ////////////////////////////////////////////////////////////////////////
   // INTERNAL STUFF BELOW HERE
   ////////////////////////////////////////////////////////////////////////
public:
   // class type definitions
   RN_CLASS_TYPE_DECLARATIONS(NeuronCellular);


protected:
   friend class NeuronData;
};



/////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/* property functions */

inline RNBoolean NeuronCellular::
IsCellular(void) const
{
   // return TRUE
   return TRUE;
}



inline RNBoolean NeuronCellular::
IsExtracellular(void) const
{
   // return FALSE
   return FALSE;
}