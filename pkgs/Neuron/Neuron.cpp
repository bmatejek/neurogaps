// Source file for the Neuron module



// include files

#include "Neuron.h"



// private variables

static int neuron_active_count = 0;



int NeuronInit(void)
{
   // check whether already initialized
   if ((neuron_active_count++) > 0) return TRUE;

   // initialize dependencies
   if (!R3InitShapes()) return FALSE;

   // return OK status
   return TRUE;
}



void NeuronStop(void)
{
   // check whether already initialized
   if ((--neuron_active_count) > 0) return;

   // stop dependencies
   R3StopShapes();
}