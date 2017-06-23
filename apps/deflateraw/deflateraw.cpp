// Source file for the deflate raw truth/supervoxel algorithm



// include files

#include "RNDataStructures/RNDataStructures.h"



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static const char *input_filename = NULL;



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
            else if (!strcmp(*argv, "-debug")) { print_verbose = 1;  print_debug = 1; }
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        else {
           if (!input_filename) input_filename = *argv;
            else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
        }
        argv++; argc--;
    }

    // make sure there is an input filename
    if (!input_filename) {
        fprintf(stderr, "Need to supply input filename.\n");
        return 0;
    }

    // return OK status 
    return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program 
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    // parse program arguments
    if (!ParseArgs(argc, argv)) exit(-1);

    char root_filename[4096];
    strncpy(root_filename, input_filename, 4096);
    char *extp = strrchr(root_filename, '.');
    *extp = '\0';

    RNMeta meta = RNReadNeuronMetaFile(input_filename);

    // set grid array
    R3Grid *grid = RNReadNeuronMetaRawFile(root_filename);
    int *grid_values = new int[grid->NEntries()];
    for (int ie = 0; ie < grid->NEntries(); ++ie) {
       grid_values[ie] = (int)(grid->GridValue(ie) + 0.5);
    }

    RNDeflateIntegerArray(grid_values, grid->NEntries());

    // update grid values
    for (int ie = 0; ie < grid->NEntries(); ++ie) {
       grid->SetGridValue(ie, grid_values[ie]);
    }

    // get raw filename
    char raw_filename[4096];
    sprintf(raw_filename, "%s.raw", root_filename);

    // write meta/raw file
    if (!RNWriteNeuronRawFile(raw_filename, meta, grid)) exit(-1);

    // print outcome
    printf("Successfully deflated %s\n", input_filename);

    // return success
    return 0;
}
