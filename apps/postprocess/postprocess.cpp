// Source file for the simple file conversion 



// include files 

#include "Neuron/Neuron.h"
#include "RNML/RNML.h"
#include <vector>
#include <string>



// program arguments

static int print_verbose = 0;
static int print_debug = 0;
static char *input_filename = NULL;
static int predict_all = 0;


// global variables

static NeuronData *nd = NULL;
static float *merge_terms = NULL;
static float *distances = NULL;
static int *prev_node = NULL;


// directory structure

static const char *algs_directory = "algs_data/graphcut";
static const char *postprocess_directory = "output/postprocess";
static const char *distances_directory = "distances";
static const char *dataset_directory = "algs_data/dataset";



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
   if (!nd->ReadFile(filename, TRUE)) {
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



static int ReadGraphcutFile(char root_filename[4096])
{
   merge_terms = new float[nd->NPredictionBoundaries()];
   if (!merge_terms) return 0;

   // read graphcut file
   char graphcut_filename[4096];
   if (predict_all) sprintf(graphcut_filename, "%s/%s_random_forest_post_all.binary", algs_directory, root_filename);
   else sprintf(graphcut_filename, "%s/%s_random_forest_post.binary", algs_directory, root_filename);

   // open file
   FILE *fp = fopen(graphcut_filename, "rb");
   if (!fp) { fprintf(stderr, "Failed to read %s\n", graphcut_filename); return 0; }

   // read the number of boundaries
   int nboundaries;
   fread(&nboundaries, sizeof(int), 1, fp);
   rn_assertion(nboundaries == nd->NPredictionBoundaries());

   if (fread(merge_terms, sizeof(float), nd->NPredictionBoundaries(), fp) != (unsigned int)nd->NPredictionBoundaries()) {
      fprintf(stderr, "Failed to read %s\n", graphcut_filename);
      return 0;
   }

   // close file
   fclose(fp);

   // return success
   return 1;
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
         else if (!strcmp(*argv, "-all")) { predict_all = 1; }
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

struct PathOutput {
   PathOutput(int next_neighbor, int nsteps, RNScalar distance_to_boundary, RNScalar total_distance) :
      next_neighbor(next_neighbor),
      nsteps(nsteps),
      distance_to_boundary(distance_to_boundary),
      total_distance(total_distance)
   {}

   // instance variables
   int next_neighbor;
   int nsteps;
   RNScalar distance_to_boundary;
   RNScalar total_distance;
};


static PathOutput UpdatePath(int voxel_index)
{
   NeuronVoxel *starting_voxel = nd->Voxel(voxel_index);
   int path_index = voxel_index;


   int nsteps = 0;
   RNScalar distance_to_boundary = 0.0;
   RNScalar total_distance = distances[voxel_index];

   while (path_index != -1) {
      NeuronVoxel *previous_voxel = nd->Voxel(path_index);
      path_index = prev_node[path_index];
      NeuronVoxel *current_voxel = nd->Voxel(path_index);

      nsteps++;
      distance_to_boundary += previous_voxel->AffinityToNeighbor(current_voxel);
      if (!current_voxel->Prediction()) return PathOutput(-1, nsteps, distance_to_boundary, total_distance);
      if (starting_voxel->Prediction() != current_voxel->Prediction()) return PathOutput(current_voxel->Prediction()->DataIndex(), nsteps, distance_to_boundary, total_distance);
   }

   // return success
   return PathOutput(-1, nsteps, distance_to_boundary, total_distance);
}



int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // file conversion
   if (!ReadData(input_filename)) exit(-1);

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in graphcut prediction
   char graphcut_filename[4096];
   sprintf(graphcut_filename, "output/graphcut/%s", root_filename);

   R3Grid *graphcut_grid = RNReadNeuronMetaRawFile(graphcut_filename);
   if (!graphcut_grid) exit(-1);

   char ordermerge_filename[4096];
   sprintf(ordermerge_filename, "output/ordermerge/%s", root_filename);

   R3Grid *ordermerge_grid = RNReadNeuronMetaRawFile(ordermerge_filename);
   if (!ordermerge_grid) exit(-1);

   // create predictions
   nd->CreatePredictions(ordermerge_grid);

   int *prediction_mapping = new int[nd->NPredictions()];
   for (int ip = 0; ip < nd->NPredictions(); ++ip)
      prediction_mapping[ip] = ip;

   // remove entire encapsulation
   int nencapsulated = 0;
   for (int ip1 = 0; ip1 < nd->NPredictions(); ++ip1) {
      NeuronPrediction *prediction = nd->Prediction(ip1);
      if (prediction->NBoundaries() != 1) continue;
      nencapsulated++;
      // get the neighbor
      NeuronPredictionBoundary *boundary = prediction->Boundary(0);
      NeuronPrediction *neighbor = boundary->OtherPrediction(prediction);

      int prediction_one_label = prediction_mapping[ip1];
      int prediction_two_label = prediction_mapping[neighbor->DataIndex()];
      for (int ip2 = 0; ip2 < nd->NPredictions(); ++ip2) {
         if (prediction_mapping[ip2] == prediction_one_label) {
            prediction_mapping[ip2] = prediction_two_label;
         }
      }
   }
   printf("%d\n", nencapsulated);

   /*std::vector<std::string> names = std::vector<std::string>();
   std::vector<unsigned int> labels = std::vector<unsigned int>();
   std::vector<std::vector<float> > attributes = std::vector<std::vector<float> >();
   std::vector<unsigned int> identifications = std::vector<unsigned int>();

   names.push_back("Prediction Score");
   names.push_back("NVoxels");
   names.push_back("NCellulars");*/

   // find non boundary predictions
   int nmerged = 0;
   int ncorrect = 0;
   int nnon_boundaries = 0;
   for (int ip = 0; ip < nd->NPredictions(); ++ip) {
      NeuronPrediction *prediction = nd->Prediction(ip);
      if (prediction->IsOnBoundary()) continue;
      nnon_boundaries++;
      // find all supervoxels in this prediction
      int original_graphcut_label = -1;
      RNBoolean same_label = TRUE;
      for (int is = 0; is < prediction->NSupervoxels(); ++is) {
         NeuronSupervoxel *supervoxel = prediction->Supervoxel(is);
         NeuronVoxel *voxel = supervoxel->CenterVoxel();

         // see if they all have the same graphcut labeling
         int graphcut_label = (int)(graphcut_grid->GridValue(voxel->DataIndex()) + 0.5);

         if (is == 0) original_graphcut_label = graphcut_label;
         else if (graphcut_label != original_graphcut_label) same_label = FALSE;
      }

      // find the neighbor with this graphcut label
      if (!same_label) continue;

      int *prediction_counter = new int[nd->NPredictions()];
      for (int ipp = 0; ipp < nd->NPredictions(); ++ipp)
         prediction_counter[ipp] = 0;

      for (int is = 0; is < prediction->NSupervoxels(); ++is) {
         NeuronSupervoxel *supervoxel = prediction->Supervoxel(is);
         for (int ib = 0; ib < supervoxel->NBoundaries(); ++ib) {
            NeuronBoundary *boundary = supervoxel->Boundary(ib);

            // get neighbor
            NeuronSupervoxel *neighbor = boundary->OtherSupervoxel(supervoxel);
            if (neighbor->IsExtracellular()) continue;
            NeuronVoxel *neighbor_voxel = neighbor->CenterVoxel();
            if (neighbor_voxel->Prediction() == prediction) continue;

            // get neighbor prediction value
            int graphcut_prediction_value = (int)(graphcut_grid->GridValue(neighbor_voxel->DataIndex()) + 0.5);

            // ignore if this does not continue the graphcut
            if (graphcut_prediction_value != original_graphcut_label) continue;

            // get the prediction value for here
            int prediction_value = neighbor_voxel->Prediction()->DataIndex();
            prediction_counter[prediction_value]++;
         }
      }

      int total_boundaries = 0;
      int best_prediction_index = -1;
      int best_prediction_score = 0;
      for (int ipp = 0; ipp < nd->NPredictions(); ++ipp) {
         if (prediction_counter[ipp] > best_prediction_score) {
            best_prediction_score = prediction_counter[ipp];
            best_prediction_index = ipp;
         }
         total_boundaries += prediction_counter[ipp];
      }
      if (best_prediction_index == -1) continue;
      int label = (prediction->MajorityHumanLabel() == nd->Prediction(best_prediction_index)->MajorityHumanLabel());

      if (best_prediction_score != total_boundaries) continue;

      if (label) ncorrect++;

      //if (!label) continue;

      //attributes.push_back(std::vector<float>());
      //labels.push_back(label);
      //identifications.push_back(ip);

      // add features
      //attributes[nentries].push_back(best_prediction_score);
      //attributes[nentries].push_back(prediction->NVoxels());
      //attributes[nentries].push_back(prediction->NSupervoxels());

      // get the prediction values
      int prediction_one_value = prediction_mapping[ip];
      int prediction_two_value = prediction_mapping[best_prediction_index];

      for (int ipp = 0; ipp < nd->NPredictions(); ++ipp) {
         if (prediction_mapping[ipp] == prediction_one_value)
            prediction_mapping[ipp] = prediction_two_value;
      }
      ++nmerged;
   }

   printf("NCorrect: %d/%d (%0.2lf)\n", ncorrect, nmerged, ncorrect / (RNScalar)nmerged);
   printf("%d/%d\n", nmerged, nnon_boundaries);

   //RNDataset dataset = RNDataset(attributes, labels, names, identifications, 2);
   //dataset.CalculateEntropy();

   int *voxel_proposals = new int[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->Prediction()) voxel_proposals[iv] = prediction_mapping[voxel->Prediction()->DataIndex()];
      else voxel_proposals[iv] = 0;
      rn_assertion((0 <= voxel_proposals[iv]) && (voxel_proposals[iv] < nd->NPredictions()));
   }

   RNScalar test_metric[12];
   nd->TestMetric(voxel_proposals, test_metric, TRUE);

   // save post processed grid
   R3Grid *postprocessed = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      postprocessed->SetGridValue(iv, voxel_proposals[iv]);
   }

   RNMeta meta = RNMeta("Uint16", nd->XResolution(), nd->YResolution(), nd->ZResolution(), 1);

   // get output filename
   char output_filename[4096];
   sprintf(output_filename, "%s/%s", postprocess_directory, root_filename);
   if (!RNWriteNeuronMetaRawFile(output_filename, meta, postprocessed)) exit(-1);

   // free memory
   delete[] voxel_proposals;
   delete[] prediction_mapping;
   delete ordermerge_grid;
   delete graphcut_grid;
   delete postprocessed;
   delete nd;

   // return success
   return 0;
}
