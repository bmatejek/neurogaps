// Include file for the neuron human label class



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class NeuronHumanLabel : public NeuronReconstruction {
public:
   //////////////////////////////////////////
   //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
   //////////////////////////////////////////

   // constructor/destructor functions
   NeuronHumanLabel(void);
   virtual ~NeuronHumanLabel(void);


   ////////////////////////////////
   //// MANIPULATION FUNCTIONS ////
   ////////////////////////////////
protected:
   // supervoxel manipulation functions
   virtual void InsertSupervoxel(NeuronSupervoxel *supervoxel);
   virtual void RemoveSupervoxel(NeuronSupervoxel *supervoxel);


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
   RN_CLASS_TYPE_DECLARATIONS(NeuronHumanLabel);

protected:
   friend class NeuronData;
};