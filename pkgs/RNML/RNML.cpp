// Source file for the machine learning module



// include files

#include "RNML.h"



// private variables

static int rnml_active_count = 0;



int RNMLInit(void)
{
   // check whether already initialized
   if ((rnml_active_count++) > 0) return TRUE;

   // initialize dependencies
   if (!R3InitShapes()) return FALSE;

   // return OK status
   return TRUE;
}



void RNMLStop(void)
{
   // check whether already initialized
   if ((--rnml_active_count) > 0) return;

   // stop dependencies
   R3StopShapes();
}
