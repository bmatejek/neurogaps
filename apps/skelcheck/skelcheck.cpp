// Source file for the subsection algorithm



// include files 

#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static char *human_label_filename = NULL;
static char *skeleton_filename = NULL;


////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
   // parse arguments
   argc--; argv++;
   while (argc > 0) {
      if ((*argv)[0] == '-') {
         if (!strcmp(*argv, "-v")) print_verbose = 1;
         else if (!strcmp(*argv, "-human_labels")) { argv++; argc--; human_label_filename = *argv; }
         else if (!strcmp(*argv, "-skeleton")) { argv++; argc--; skeleton_filename = *argv; }
         else if (!strcmp(*argv, "-debug")) { print_debug = 1; print_verbose = 1; }
         // number of files to generate
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0;
      }
      argv++; argc--;
   }

   // make sure there is an input filename
   if (!human_label_filename) { fprintf(stderr, "Need to supply human label filename.\n");  return 0; }
   if (!skeleton_filename) { fprintf(stderr, "Need to supply skeleton filename.\n"); return 0; }

   // return OK status 
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // get the root filenames
   char skeleton_root_filename[4096];
   strncpy(skeleton_root_filename, skeleton_filename, 4096);
   char *skeleton_extp = strrchr(skeleton_root_filename, '.');
   *skeleton_extp = '\0';

   char human_label_root_filename[4096];
   strncpy(human_label_root_filename, human_label_filename, 4096);
   char *human_label_extp = strrchr(human_label_root_filename, '.');
   *human_label_extp = '\0';
   
   R3Grid *skeleton = RNReadNeuronMetaRawFile(skeleton_root_filename);
   R3Grid *human_label = RNReadNeuronMetaRawFile(human_label_root_filename);

   // just checking
   rn_assertion(skeleton->NEntries() == human_label->NEntries());

   // see how many voxels match
   unsigned int ncorrect = 0; 
   unsigned int nincorrect = 0;
   for (int ie = 0; ie < skeleton->NEntries(); ++ie) {
      int skeleton_value = (int)(skeleton->GridValue(ie) + 0.5);
      if (skeleton_value) {
         int human_label_value = (int)(human_label->GridValue(ie) + 0.5);
         if (human_label_value) ncorrect++;
         else nincorrect++;
      }
   }

   printf("Proportion of skeleton in human labels: %lf\n", ncorrect / (RNScalar)(ncorrect + nincorrect));

   // free memory
   delete skeleton;
   delete human_label;

   // return success
   return 0;
}