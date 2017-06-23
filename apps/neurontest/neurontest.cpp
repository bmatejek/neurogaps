// Source file for the .neuron test algorithm



// include files 

#include "Neuron/Neuron.h"
#include <vector>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;



// global variables

static NeuronData *nd;
static char *input_filename = NULL;
static R3Grid *affinity_grid[3] = { NULL, NULL, NULL };
static R3Grid *image_grid = NULL;
static R3Grid *human_labels_grid = NULL;
static R3Grid *machine_labels_grid = NULL;



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(const char *filename)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate memory for neuron data\n");
      return 0;
   }

   // read in the file
   if (!nd->ReadFile(filename, TRUE, TRUE)) {
      fprintf(stderr, "Failed to read file\n");
      return 0;
   }

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }
   if (print_debug) {
      printf("  Bounding Box: (%0.2f, %0.2f, %0.2f) to (%0.2f, %0.2f, %0.2f)\n", nd->GridBox().XMin(), nd->GridBox().YMin(), nd->GridBox().ZMin(), nd->GridBox().XMax(), nd->GridBox().YMax(), nd->GridBox().ZMax());
      printf("  Voxels: %d\n", nd->NVoxels());
      printf("  Cellulars: %d\n", nd->NCellulars());
      printf("  Extracellulars: %d\n", nd->NExtracellulars());
      printf("  Boundaries: %d\n", nd->NBoundaries());
      printf("  Human Labels: %d\n", nd->NHumanLabels());
      printf("  Predictions: %d\n", nd->NPredictions());
   }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Test functions
////////////////////////////////////////////////////////////////////////

void TestVoxelIndices(void)
{
   RNTime start_time;
   start_time.Read();

   // index to indices functions and vice versa
   printf("Testing all index and indices functions...");
   fflush(stdout);

   int iv = 0;
   for (int iz = 0; iz < nd->ZResolution(); ++iz) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int ix = 0; ix < nd->XResolution(); ++ix, ++iv) {
            int coordinates[3] = { ix, iy, iz };

            // first function
            int index = nd->IndicesToIndex(ix, iy, iz);
            rn_assertion(index == iv);

            // second function
            index = nd->IndicesToIndex(coordinates);
            rn_assertion(index == iv);

            // third function
            nd->IndicesToIndex(ix, iy, iz, index);
            rn_assertion(index == iv);

            // fourth function
            nd->IndicesToIndex(coordinates, index);
            rn_assertion(index == iv);

            // fifth function
            int testx, testy, testz;
            nd->IndexToIndices(iv, testx, testy, testz);
            rn_assertion((testx == ix) && (testy == iy) && (testz == iz));

            // sixth function
            int test_coordinates[3];
            nd->IndexToIndices(iv, test_coordinates);
            rn_assertion((test_coordinates[RN_X] == ix) && (test_coordinates[RN_Y] == iy) && (test_coordinates[RN_Z] == iz));

            // test NeuronData::Voxel()
            NeuronVoxel *voxel_one = nd->Voxel(ix, iy, iz);
            NeuronVoxel *voxel_two = nd->Voxel(iv);
            rn_assertion(voxel_one == voxel_two);
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestVolumeLabels(void)
{
   RNTime start_time;
   start_time.Read();

   // test voxel extracellular/supervoxel label
   printf("Testing voxel cellular/extracellular labels...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);

            int voxel_index = nd->IndicesToIndex(ix, iy, iz);
            int voxel_mapping = nd->VoxelMapping(voxel_index);
            int volume_index = voxel->Supervoxel()->SupervoxelIndex();
            rn_assertion(volume_index == voxel_mapping);

            RNScalar supervoxel_value = machine_labels_grid->GridValue(ix, iy, iz);
            rn_assertion(supervoxel_value >= 0);
            int supervoxel_index = (int)(supervoxel_value + 0.5) - 1;

            rn_assertion(supervoxel_index != -1 || voxel->IsExtracellular());
            rn_assertion(supervoxel_index == -1 || voxel->IsCellular());

            if (supervoxel_index >= 0) {
               int cellular_data_index = voxel->Supervoxel()->DataIndex();
               rn_assertion(nd->Cellular(cellular_data_index) == voxel->Supervoxel());
            }
            else {
               int extracellular_data_index = voxel->Supervoxel()->DataIndex();
               rn_assertion(nd->Extracellular(extracellular_data_index) == voxel->Supervoxel());
            }
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestBoundaryFunctions(void)
{
   RNTime start_time;
   start_time.Read();

   // test on boundary functions
   printf("Testing on boundary voxel functions...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);
            if (ix == 0) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            if (iy == 0) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            if (iz == 0) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            if (ix == nd->XResolution() - 1) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            if (iy == nd->YResolution() - 1) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            if (iz == nd->ZResolution() - 1) {
               rn_assertion(voxel->IsOnBoundary());
               continue;
            }
            rn_assertion(!voxel->IsOnBoundary());
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestNeighborFunctions(void)
{
   RNTime start_time;
   start_time.Read();

   // test previous/next/neighbor voxels
   printf("Testing neighbor/next/previous functions...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);

            if (ix != 0) {
               NeuronVoxel *neighbor = nd->Voxel(ix - 1, iy, iz);
               NeuronVoxel *neighbor2 = voxel->Neighbor(0);
               NeuronVoxel *neighbor3 = voxel->PreviousVoxel(RN_X);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
            if (ix != nd->XResolution() - 1) {
               NeuronVoxel *neighbor = nd->Voxel(ix + 1, iy, iz);
               NeuronVoxel *neighbor2 = voxel->Neighbor(1);
               NeuronVoxel *neighbor3 = voxel->NextVoxel(RN_X);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
            if (iy != 0) {
               NeuronVoxel *neighbor = nd->Voxel(ix, iy - 1, iz);
               NeuronVoxel *neighbor2 = voxel->Neighbor(2);
               NeuronVoxel *neighbor3 = voxel->PreviousVoxel(RN_Y);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
            if (iy != nd->YResolution() - 1) {
               NeuronVoxel *neighbor = nd->Voxel(ix, iy + 1, iz);
               NeuronVoxel *neighbor2 = voxel->Neighbor(3);
               NeuronVoxel *neighbor3 = voxel->NextVoxel(RN_Y);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
            if (iz != 0) {
               NeuronVoxel *neighbor = nd->Voxel(ix, iy, iz - 1);
               NeuronVoxel *neighbor2 = voxel->Neighbor(4);
               NeuronVoxel *neighbor3 = voxel->PreviousVoxel(RN_Z);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
            if (iz != nd->ZResolution() - 1) {
               NeuronVoxel *neighbor = nd->Voxel(ix, iy, iz + 1);
               NeuronVoxel *neighbor2 = voxel->Neighbor(5);
               NeuronVoxel *neighbor3 = voxel->NextVoxel(RN_Z);
               rn_assertion((neighbor == neighbor2) && (neighbor2 == neighbor3));
            }
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestNeighborAffinities(void)
{
   RNTime start_time;
   start_time.Read();

   // test all neighbor voxels
   printf("Testing neighbor affinity functions...");
   fflush(stdout);

   for (int dim = 0; dim <= 2; ++dim) {
      for (int ix = 0; ix < nd->XResolution(); ++ix) {
         for (int iy = 0; iy < nd->YResolution(); ++iy) {
            for (int iz = 0; iz < nd->ZResolution(); ++iz) {
               NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);
               int prevx = ix - 1;
               int prevy = iy - 1;
               int prevz = iz - 1;
               int nextx = ix + 1;
               int nexty = iy + 1;
               int nextz = iz + 1;

               if (ix != 0) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(prevx, iy, iz));
                  rn_assertion(affinity == affinity_grid[RN_X]->GridValue(prevx, iy, iz));
               }
               if (iy != 0) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(ix, prevy, iz));
                  rn_assertion(affinity == affinity_grid[RN_Y]->GridValue(ix, prevy, iz));
               }
               if (iz != 0) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(ix, iy, prevz));
                  rn_assertion(affinity == affinity_grid[RN_Z]->GridValue(ix, iy, prevz));
               }
               if (ix != nd->XResolution() - 1) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(nextx, iy, iz));
                  rn_assertion(affinity == affinity_grid[RN_X]->GridValue(ix, iy, iz));
               }
               if (iy != nd->YResolution() - 1) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(ix, nexty, iz));
                  rn_assertion(affinity == affinity_grid[RN_Y]->GridValue(ix, iy, iz));
               }
               if (iz != nd->ZResolution() - 1) {
                  RNScalar affinity = voxel->AffinityToNeighbor(nd->Voxel(ix, iy, nextz));
                  rn_assertion(affinity == affinity_grid[RN_Z]->GridValue(ix, iy, iz));
               }
            }
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestImageTruthValues(void)
{
   RNTime start_time;
   start_time.Read();

   //  make sure the image values are correct
   printf("Testing voxel human label values...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);
            int human_label_index = (int)(human_labels_grid->GridValue(ix, iy, iz) + 0.5) - 1;
            if (voxel->HumanLabel()) { rn_assertion(voxel->HumanLabel()->DataIndex() == human_label_index); }
            else { rn_assertion(human_label_index == -1); }
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestImageValues(void)
{
   RNTime start_time;
   start_time.Read();

   //  make sure the image values are correct
   printf("Testing voxel image values...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);

            rn_assertion(voxel->Image() == image_grid->GridValue(ix, iy, iz));
         }
      }
   }
}

void TestDepthFirstSearch(void)
{
   RNTime start_time;
   start_time.Read();

   // test extracellular depth first search
   printf("Testing extracellular depth first search...");
   fflush(stdout);

   for (int ix = 0; ix < nd->XResolution(); ++ix) {
      for (int iy = 0; iy < nd->YResolution(); ++iy) {
         for (int iz = 0; iz < nd->ZResolution(); ++iz) {
            NeuronVoxel *voxel = nd->Voxel(ix, iy, iz);
            if (voxel->IsCellular()) continue;
            NeuronExtracellular *extracellular = (NeuronExtracellular *)voxel->Supervoxel();
            for (int in = 0; in < voxel->NNeighbors(); ++in) {
               NeuronVoxel *neighbor = voxel->Neighbor(in);
               if (!neighbor) continue;
               if (neighbor->IsCellular()) continue;

               NeuronExtracellular *neighbor_extracellular = (NeuronExtracellular *)neighbor->Supervoxel();
               rn_assertion(extracellular == neighbor_extracellular);
            }
         }
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestGlobalVoxel(void)
{
   RNTime start_time;
   start_time.Read();

   // make sure it is possible to access each unique voxel
   printf("Testing one to one relationship with all voxels (test GlobalVoxel)...");
   fflush(stdout);

   RNBoolean *reached = new RNBoolean[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      reached[iv] = FALSE;
   }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      rn_assertion(voxel->DataIndex() == iv);
      reached[voxel->DataIndex()] = TRUE;
   }
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      rn_assertion(reached[iv]);
   }
   delete[] reached;
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestVolumeIndices(void)
{
   RNTime start_time;
   start_time.Read();

   // test volume/data index
   printf("Testing volume/data index relations...");
   fflush(stdout);

   for (int is = 0; is < nd->NSupervoxels(); ++is) {
      NeuronSupervoxel *supervoxel = nd->Supervoxel(is);
      if (supervoxel->IsCellular()) { rn_assertion(supervoxel->DataIndex() == supervoxel->SupervoxelIndex()); }
      else { rn_assertion(supervoxel->DataIndex() + nd->NCellulars() == supervoxel->SupervoxelIndex()); }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestUniqueVolume(void)
{
   RNTime start_time;
   start_time.Read();

   // test index uniqueness for one to one
   printf("Testing unique volume and data indices as well as size...");
   fflush(stdout);

   RNBoolean *cellular = new RNBoolean[nd->NCellulars()];
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      cellular[ic] = FALSE;
   }
   RNBoolean *extracellular = new RNBoolean[nd->NExtracellulars()];
   for (int ie = 0; ie < nd->NExtracellulars(); ++ie) {
      extracellular[ie] = FALSE;
   }
   RNBoolean *supervoxel = new RNBoolean[nd->NSupervoxels()];
   for (int is = 0; is < nd->NSupervoxels(); ++is) {
      supervoxel[is] = FALSE;
   }
   for (int is = 0; is < nd->NSupervoxels(); ++is) {
      // make sure the volume is nonempty
      rn_assertion(nd->Supervoxel(is)->NVoxels());

      if (nd->Supervoxel(is)->IsCellular()) {
         int cellular_index = nd->Supervoxel(is)->DataIndex();
         rn_assertion((0 <= cellular_index) && (cellular_index < nd->NSupervoxels()));
         cellular[cellular_index] = TRUE;
      }
      else {
         int extracellular_index = nd->Supervoxel(is)->DataIndex();
         rn_assertion((0 <= extracellular_index) && (extracellular_index < nd->NExtracellulars()));
         extracellular[extracellular_index] = TRUE;
      }
      int supervoxel_index = nd->Supervoxel(is)->SupervoxelIndex();
      rn_assertion((0 <= supervoxel_index) && (supervoxel_index < nd->NSupervoxels()));
      supervoxel[supervoxel_index] = TRUE;
   }
   for (int ic = 0; ic < nd->NCellulars(); ++ic) {
      rn_assertion(cellular[ic]);
   }
   for (int ie = 0; ie < nd->NExtracellulars(); ++ie) {
      rn_assertion(extracellular[ie]);
   }
   for (int is = 0; is < nd->NSupervoxels(); ++is) {
      rn_assertion(supervoxel[is]);
   }
   delete[] cellular;
   delete[] extracellular;
   delete[] supervoxel;
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestUniqueTruthSegment(void)
{
   RNTime start_time;
   start_time.Read();

   // test index uniqueness for one to one
   printf("Testing unique truth and segment indices as well as size...");
   fflush(stdout);

   RNBoolean *human_label = new RNBoolean[nd->NHumanLabels()];
   for (int ih = 0; ih < nd->NHumanLabels(); ++ih) {
      human_label[ih] = FALSE;
   }
   RNBoolean *prediction = new RNBoolean[nd->NPredictions()];
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      prediction[ip] = FALSE;
   }
   for (int ih = 0; ih < nd->NHumanLabels(); ++ih) {
      // make sure the volume is nonempty
      rn_assertion(nd->HumanLabel(ih)->NVoxels());
      int human_label_index = nd->HumanLabel(ih)->DataIndex();
      rn_assertion((0 <= human_label_index) && (human_label_index < nd->NHumanLabels()));
      human_label[human_label_index] = TRUE;
   }
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      // make sure the volume is nonempty
      rn_assertion(nd->Prediction(ip)->NVoxels());
      int prediction_index = nd->Prediction(ip)->DataIndex();
      rn_assertion((0 <= prediction_index) && (prediction_index < nd->NPredictions()));
      prediction[prediction_index] = TRUE;
   }
   for (int ih = 0; ih < nd->NHumanLabels(); ++ih) {
      rn_assertion(human_label[ih]);
   }
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      rn_assertion(prediction[ip]);
   }
   delete[] human_label;
   delete[] prediction;
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}

void TestVolumeUniqueness(void)
{
   RNTime start_time;
   start_time.Read();

   // test uniqueness of supervoxels in segments/truths
   printf("Testing uniqueness of cellulars/extracellulars in predictions/human labels...");
   fflush(stdout);

   for (int ih = 0; ih < nd->NHumanLabels(); ++ih) {
      NeuronHumanLabel *human_label = nd->HumanLabel(ih);
      std::vector<NeuronSupervoxel *> supervoxels = std::vector<NeuronSupervoxel *>();
      for (int is = 0; is < human_label->NSupervoxels(); ++is) {
         NeuronSupervoxel *supervoxel = human_label->Supervoxel(is);
         for (unsigned int ia = 0; ia < supervoxels.size(); ++ia) {
            rn_assertion(supervoxel != supervoxels[ia]);
         }
         supervoxels.push_back(supervoxel);
      }
   }
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      NeuronPrediction *prediction = nd->Prediction(ip);
      std::vector<NeuronSupervoxel *> supervoxels = std::vector<NeuronSupervoxel *>();
      for (int is = 0; is < prediction->NSupervoxels(); ++is) {
         NeuronSupervoxel *supervoxel = prediction->Supervoxel(is);
         for (unsigned int ia = 0; ia < supervoxels.size(); ++ia) {
            rn_assertion(supervoxel != supervoxels[ia]);
         }
         supervoxels.push_back(supervoxel);
      }
   }
   printf("done in %.2f seconds.\n", start_time.Elapsed());
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int ParseArgs(int argc, char **argv)
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

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // get all of the raw grids
   RNTime grid_time;
   grid_time.Read();
   for (int dim = 0; dim <= 2; ++dim)
      affinity_grid[dim] = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!nd->ReadAffinities(affinity_grid)) exit(-1);
   human_labels_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!nd->ReadHumanLabels(human_labels_grid)) exit(-1);
   image_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!nd->ReadImages(image_grid)) exit(-1);
   machine_labels_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!nd->ReadMachineLabels(machine_labels_grid)) exit(-1);
   printf("Read in grids in %.2f seconds\n\n", grid_time.Elapsed());


   //////////////////////////////
   //// RUN ALL OF THE TESTS ////
   //////////////////////////////

   TestVoxelIndices();
   TestVolumeLabels();
   TestBoundaryFunctions();
   TestNeighborFunctions();
   TestNeighborAffinities();
   TestImageTruthValues();
   TestDepthFirstSearch();
   TestGlobalVoxel();
   TestVolumeIndices();
   TestUniqueVolume();
   TestUniqueTruthSegment();
   TestVolumeUniqueness();

   // free up memory
   delete nd;

   // return success
   return 0;
}
