// Source file for neuron visualizer



// include files

#include "Neuron/Neuron.h"
#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"



// GLUT defines

#define ENTER 13
#define ESCAPE 27
#define SPACEBAR 32



// program arguments

static char *input_name = NULL;
static int print_verbose = 0;
static int print_debug = 0;
static const char *extension = NULL;



// program variables

static NeuronData *nd = NULL;
static float *distances = NULL;
static int *prev_node = NULL;
static R3Viewer *viewer = NULL;
static R3Point selected_position;
static R3Box world_box;
static int selected_voxel[3];



// voxel grids

static R3Grid *affinity_grid[3] = { NULL, NULL, NULL };
static R3Grid *human_labels_grid = NULL;
static R3Grid *image_grid = NULL;
static R3Grid *machine_labels_grid = NULL;
static R3Grid *predictions_grid = NULL;



// projection variables

static R2Grid *affinity_selected_slice[3] = { NULL, NULL, NULL };
static R2Grid *human_labels_selected_slice = NULL;
static R2Grid *image_selected_slice = NULL;
static R2Grid *machine_labels_selected_slice = NULL;
static R2Grid *predictions_selected_slice = NULL;
static RNInterval affinity_range[3];
static RNInterval human_labels_range;
static RNInterval image_range;
static RNInterval machine_labels_range;
static RNInterval predictions_range;
// always update the grid and the interval
static R2Grid *selected_slice = NULL;
static R2Point selected_slice_position(RN_UNKNOWN, RN_UNKNOWN);
static RNInterval selected_slice_range(0, 0);
static RNScalar projection_scale = FLT_MAX;



// GLUT variables

static int GLUTwindow = 0;
static int GLUTwindow_height = 700;
static int GLUTwindow_width = 1200;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;



// color arrays

static RNScalar background_color[] = { 1, 1, 1 };


// projection color variables

static int color_type = 0; // 0=gray, 1=red-green-blue



// display projection variables

static int show_projection_affinity = 0;
static int show_projection_image = 1;
static int show_projection_human_labels = 0;
static int show_projection_machine_labels = 0;
static int show_projection_predictions = 0;
static int projection_dim = RN_Z;
static int affinity_dim = RN_X;
static int selected_slice_index = 0;



// display dimenson variables

static int show_bbox = 1;
static int show_slice = 0;
static int projection = 0;
static int show_prediction_results = 1;
static int show_dimension_affinity = 0;
static int show_dimension_image = 1;
static int show_dimension_human_labels = 0;
static int show_dimension_machine_labels = 0;
static int show_dimension_predictions = 0;




// size variables

static RNScalar voxel_size = 5.0;
static int *voxel_hash = NULL;



// path variables

std::vector<int> dijkstra_path = std::vector<int>();
std::vector<RNRgb> path_colors = std::vector<RNRgb>();



// directory structure

static const char *output_directory = "output";
static const char *distances_directory = "distances";



////////////////////////////////////////////////////////////////////////
// Path algorithms
////////////////////////////////////////////////////////////////////////

static int UpdatePath(void)
{
   // reset vectors
   dijkstra_path.clear();
   path_colors.clear();

   // start from selected voxel
   int iv = nd->IndicesToIndex(selected_voxel);

   int human_label_index = (int)(human_labels_grid->GridValue(iv) + 0.5);

   int path_index = iv;
   while (path_index != -1) {
      // get this human label
      int path_human_label_index = (int)(human_labels_grid->GridValue(path_index) + 0.5);

      dijkstra_path.push_back(path_index);

      if (human_label_index == 0)
         path_colors.push_back(RNblack_rgb);
      else if (human_label_index == path_human_label_index)
         path_colors.push_back(RNgreen_rgb);
      else
         path_colors.push_back(RNred_rgb);

      path_index = prev_node[path_index];
   }

   // return success
   return 1;
}





////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static R2Point ConvertGridToWorld(R2Point grid_point)
{
   // projection offsets
   int window_xdim = (projection_dim + 1) % 3;
   int window_ydim = (projection_dim + 2) % 3;
   int xoffset = (GLUTwindow_width - nd->Scale(window_xdim) * projection_scale * nd->Resolution(window_xdim)) / 2;
   int yoffset = (GLUTwindow_height - nd->Scale(window_ydim) * projection_scale * nd->Resolution(window_ydim)) / 2;

   // convert from grid to world coordinates
   return R2Point(grid_point.X() * projection_scale * nd->Scale(window_xdim) + xoffset, grid_point.Y() * projection_scale * nd->Scale(window_ydim) + yoffset);
}



static R2Point ConvertWorldToGrid(R2Point world_point)
{
   // projection offsets
   int window_xdim = (projection_dim + 1) % 3;
   int window_ydim = (projection_dim + 2) % 3;
   int xoffset = (GLUTwindow_width - nd->Scale(window_xdim) * projection_scale * nd->Resolution(window_xdim)) / 2;
   int yoffset = (GLUTwindow_height - nd->Scale(window_ydim) * projection_scale * nd->Resolution(window_ydim)) / 2;

   // convert from world to grid coordinates
   return R2Point((world_point.X() - xoffset) / (nd->Scale(window_xdim) * projection_scale), (world_point.Y() - yoffset) / (nd->Scale(window_ydim) * projection_scale));
}



static void UpdateSlices(void)
{
   // update selected affinities
   for (int dim = 0; dim <= 2; ++dim) {
      affinity_selected_slice[dim] = affinity_grid[dim]->Slice(projection_dim, selected_slice_index);
      affinity_range[dim] = affinity_selected_slice[dim]->Range();
   }

   // update selected images
   image_selected_slice = image_grid->Slice(projection_dim, selected_slice_index);
   image_range = image_selected_slice->Range();

   // update selected machine labels
   machine_labels_selected_slice = machine_labels_grid->Slice(projection_dim, selected_slice_index);
   machine_labels_range = machine_labels_selected_slice->Range();

   // update selected truths
   human_labels_selected_slice = human_labels_grid->Slice(projection_dim, selected_slice_index);
   human_labels_range = human_labels_selected_slice->Range();

   // update predictions
   predictions_selected_slice = predictions_grid->Slice(projection_dim, selected_slice_index);
   predictions_range = predictions_selected_slice->Range();

   // update pointers to selected slice
   if (show_projection_affinity) {
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
   }
   else if (show_projection_image) {
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
   }
   else if (show_projection_machine_labels) {
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
   }
   else if (show_projection_human_labels) {
      selected_slice = human_labels_selected_slice;
      selected_slice_range = human_labels_range;
   }
   else if (show_projection_predictions) {
      selected_slice = predictions_selected_slice;
      selected_slice_range = predictions_range;
   }
}



////////////////////////////////////////////////////////////////////////
// Input/output functions
////////////////////////////////////////////////////////////////////////

static int ReadData(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();

   // allocate new neuron data
   nd = new NeuronData();
   if (!nd) {
      fprintf(stderr, "Failed to allocate neuron data\n");
      return 0;
   }

   nd->ReadFile(input_name);

   // print statistics
   if (print_verbose) {
      printf("Read data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }

   // return success
   return 1;
}



static int ReadVoxelFiles(void) {
   // start statistics
   RNTime start_time;
   start_time.Read();

   // read affinity grids
   for (int dim = 0; dim <= 2; ++dim) {
      affinity_grid[dim] = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
      if (!affinity_grid[dim]) { fprintf(stderr, "Failed to allocate memory for R3Grid\n"); exit(-1); }
   }
   nd->ReadAffinities(affinity_grid);
   for (int dim = 0; dim <= 2; ++dim) {
      affinity_selected_slice[dim] = affinity_grid[dim]->Slice(projection_dim, selected_slice_index);
      affinity_range[dim] = affinity_selected_slice[dim]->Range();
   }

   // read image grid
   image_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!image_grid) { fprintf(stderr, "Failed to allocate memory for R3Grid\n"); exit(-1); }
   nd->ReadImages(image_grid);
   image_selected_slice = image_grid->Slice(projection_dim, selected_slice_index);
   image_range = image_selected_slice->Range();

   // read machine labels grid
   machine_labels_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!machine_labels_grid) { fprintf(stderr, "Failed to allocate memory for R3Grid\n"); exit(-1); }
   nd->ReadMachineLabels(machine_labels_grid);
   machine_labels_selected_slice = machine_labels_grid->Slice(projection_dim, selected_slice_index);
   machine_labels_range = machine_labels_selected_slice->Range();

   // read truth grid
   human_labels_grid = new R3Grid(nd->XResolution(), nd->YResolution(), nd->ZResolution());
   if (!human_labels_grid) { fprintf(stderr, "Failed to allocate memory for R3Grid\n"); exit(-1); }
   nd->ReadHumanLabels(human_labels_grid);
   human_labels_selected_slice = human_labels_grid->Slice(projection_dim, selected_slice_index);
   human_labels_range = human_labels_selected_slice->Range();

   // print statistics
   if (print_verbose) {
      printf("Read voxel grids...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
   }

   // get root filename
   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in predictions
   char predictions_filename[4096];
   sprintf(predictions_filename, "%s/%s/%s", output_directory, extension, root_filename);

   // read in meta/raw file
   predictions_grid = RNReadNeuronMetaRawFile(predictions_filename);
   predictions_selected_slice = predictions_grid->Slice(projection_dim, selected_slice_index);
   predictions_range = predictions_selected_slice->Range();

   // update pointers to selected slice
   if (show_projection_affinity) {
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
   }
   else if (show_projection_image) {
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
   }
   else if (show_projection_machine_labels) {
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
   }
   else if (show_projection_human_labels) {
      selected_slice = human_labels_selected_slice;
      selected_slice_range = human_labels_range;
   }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Drawing utility functions
////////////////////////////////////////////////////////////////////////

static RNRgb
Color(RNScalar value, RNInterval range)
{
   // Check for unknown value
   if (value == R2_GRID_UNKNOWN_VALUE) {
      if (color_type == 0) return RNRgb(1, 0.5, 0);
      else return RNblack_rgb;
   }

   // Normalize color
   RNScalar value_min = range.Min();
   RNScalar value_width = range.Diameter();
   RNScalar value_scale = (value_width > 0) ? 1.0 / value_width : 1.0;
   RNScalar normalized_value = value_scale * (value - value_min);

   // Compute color
   RNRgb c(0, 0, 0);
   if (color_type == 0) {
      c[0] = normalized_value;
      c[1] = normalized_value;
      c[2] = normalized_value;
   }
   else {
      if (normalized_value < 0.5) {
         c[0] = 1 - 2 * normalized_value;
         c[1] = 2 * normalized_value;
      }
      else {
         c[1] = 1 - 2 * (normalized_value - 0.5);
         c[2] = 2 * (normalized_value - 0.5);
      }
   }

   // Return color
   return c;
}



static void GLUTDrawText(const R2Point& position, const char *s)
{
   // draw text string s at position
   glRasterPos2d(position[0], position[1]);
   while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *(s++));
}



static void DrawSlice(void)
{
   R3Affine transformation = nd->Transformation();
   transformation.Push();

   RNLoadRgb(1.0, 1.0, 1.0);
   // draw the selected slice
   if (show_dimension_affinity) affinity_grid[affinity_dim]->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_human_labels) human_labels_grid->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_image) image_grid->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_machine_labels) machine_labels_grid->DrawSlice(projection_dim, selected_slice_index);
   if (show_dimension_predictions) predictions_grid->DrawSlice(projection_dim, selected_slice_index);

   transformation.Pop();
}



static void Draw2D(void)
{
   // prologue
   glDisable(GL_LIGHTING);

   // set projection matrix
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);

   // set model view matrix
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   // draw value at selected position
   if ((selected_slice_position.X() != RN_UNKNOWN) && (selected_slice_position.Y() != RN_UNKNOWN)) {
      int ix = (int)(selected_slice_position.X() + 0.5);
      int iy = (int)(selected_slice_position.Y() + 0.5);
      RNScalar value = selected_slice->GridValue(ix, iy);

      // create hover text
      char buffer[1024];
      if (value != R2_GRID_UNKNOWN_VALUE) sprintf(buffer, "%d %d: %g", ix, iy, value);
      else sprintf(buffer, "%d %d: %s", ix, iy, "Unknown");
      RNLoadRgb(RNblue_rgb);

      // move point to world location
      R2Point world_position = ConvertGridToWorld(selected_slice_position);
      GLUTDrawText(world_position + 2 * R2ones_vector, buffer);
   }

   // draw the actual slice
   for (int iy = 1; iy < selected_slice->YResolution(); ++iy) {
      glBegin(GL_TRIANGLE_STRIP);
      for (int ix = 1; ix < selected_slice->XResolution(); ++ix) {
         for (int k = -1; k <= 0; ++k) {
            int iv = nd->IndicesToIndex(ix, iy + k, selected_slice_index);

            // convert from world to grid coordinates
            R2Point grid_position = R2Point(ix, iy + k);
            R2Point world_position = ConvertGridToWorld(grid_position);

            // get the color for this area
            RNRgb color;

            if (show_prediction_results) {
               if (voxel_hash[iv] == 0) color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
               else if (voxel_hash[iv] == 1) color = RNgreen_rgb;
               else if (voxel_hash[iv] == 2) color = RNblue_rgb;
               else if (voxel_hash[iv] == 3) color = RNred_rgb;
               else if (voxel_hash[iv] == 4) color = RNmagenta_rgb;
            }
            else {
               if (show_projection_human_labels) {
                  if (selected_slice->GridValue(ix, iy + k) < 0.5) color = RNRgb(1, 0.5, 0);
                  else color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
               }
               else {
                  color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
               }
            }

            RNLoadRgb(color);

            // add the world position vertex
            glVertex2i(world_position.X(), world_position.Y());
         }
      }
      glEnd();
   }

   // reset projection matrix
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();

   // reset model view matrix
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();

   // epilogue
   glEnable(GL_LIGHTING);
}



static void Draw3D(void)
{
   // set viewing transformation
   viewer->Camera().Load();

   // set lights
   static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
   glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
   static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
   glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

   // prologue
   glDisable(GL_LIGHTING);

   // draw selected position as a sphere
   RNLoadRgb(RNblue_rgb);
   R3Sphere(selected_position, voxel_size / 2).Draw();

   // draw neuron data bounding box
   if (show_bbox) {
      RNLoadRgb(RNblack_rgb);
      world_box.Outline();
   }

   if (show_slice) DrawSlice();

   // draw selected position
   RNLoadRgb(RNblue_rgb);
   R3Sphere(selected_position, voxel_size).Draw();

   for (unsigned int iv = 1; iv < dijkstra_path.size(); ++iv) {
      RNLoadRgb(path_colors[iv]);
      int ix, iy, iz;
      nd->IndexToIndices(dijkstra_path[iv], ix, iy, iz);
      R3Point position = R3Point(ix, iy, iz);
      position.Transform(nd->Transformation());
      R3Sphere(position, voxel_size).Draw();
   }


   // epilogue
   glEnable(GL_LIGHTING);
}



////////////////////////////////////////////////////////////////////////
// GLUT interface functions
////////////////////////////////////////////////////////////////////////

void GLUTStop(void)
{
   // destroy window
   glutDestroyWindow(GLUTwindow);

   // delete the neuron data
   RNTime start_time;
   start_time.Read();

   delete nd;

   // print statistics
   if (print_verbose) {
      printf("Deleted data...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      fflush(stdout);
   }

   // exit
   exit(0);
}



void GLUTRedraw(void)
{
   // clear window
   glClearColor(background_color[0], background_color[1], background_color[2], 1.0);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // drawing varies on projection
   if (projection) Draw2D();
   else Draw3D();

   // swap buffers
   glutSwapBuffers();
}



void GLUTResize(int w, int h)
{
   // resize window
   glViewport(0, 0, w, h);

   // resize viewer viewport
   viewer->ResizeViewport(0, 0, w, h);

   // remember window size
   GLUTwindow_width = w;
   GLUTwindow_height = h;

   // redraw
   glutPostRedisplay();
}



void GLUTMotion2D(int x, int y)
{
   if (GLUTbutton[0]) {
      // query
      R2Point world_point = R2Point(x, y);
      selected_slice_position = ConvertWorldToGrid(world_point);

      // get window dimensions
      int window_xdim = (projection_dim + 1) % 3;
      int window_ydim = (projection_dim + 2) % 3;

      // get the grid x and y coordinates
      int gridx = (int)(selected_slice_position.X() + 0.5);
      int gridy = (int)(selected_slice_position.Y() + 0.5);

      // cannot select a point outside of the image
      if (gridx < 0 || gridx >= nd->Resolution(window_xdim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
      else if (gridy < 0 || gridy >= nd->Resolution(window_ydim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

      glutPostRedisplay();
   }
   else if (GLUTbutton[1]) {
      // zoom
   }
   else if (GLUTbutton[2]) {
      // move voxel
      R2Point world_point = R2Point(x, y);
      R2Point grid_point = ConvertWorldToGrid(world_point);

      int window_xdim = (projection_dim + 1) % 3;
      int window_ydim = (projection_dim + 2) % 3;

      if (grid_point.X() < 0) grid_point.SetCoord(RN_X, 0);
      if (grid_point.Y() < 0) grid_point.SetCoord(RN_Y, 0);
      if (grid_point.X() > nd->Resolution(window_xdim) - 1) grid_point.SetCoord(RN_X, nd->Resolution(window_xdim) - 1);
      if (grid_point.Y() > nd->Resolution(window_ydim) - 1) grid_point.SetCoord(RN_Y, nd->Resolution(window_ydim) - 1);

      R3Point intersection;
      if (projection_dim == RN_X) intersection = R3Point(selected_slice_index, grid_point.X(), grid_point.Y());
      else if (projection_dim == RN_Y) intersection = R3Point(grid_point.Y(), selected_slice_index, grid_point.X());
      else intersection = R3Point(grid_point.X(), grid_point.Y(), selected_slice_index);

      // update the global selected position
      selected_position = intersection;
      selected_position.Transform(nd->Transformation());

      // update the selected voxel
      selected_voxel[RN_X] = (int)(intersection.X() + 0.5);
      selected_voxel[RN_Y] = (int)(intersection.Y() + 0.5);
      selected_voxel[RN_Z] = (int)(intersection.Z() + 0.5);

      UpdatePath();

      glutPostRedisplay();
   }
}



void GLUTMotion3D(int x, int y)
{
   // compute mouse movement
   int dx = x - GLUTmouse[0];
   int dy = y - GLUTmouse[1];

   // world in hand navigation
   R3Point origin = world_box.Centroid();
   if (GLUTbutton[0]) {
      viewer->RotateWorld(1.0, origin, x, y, dx, dy);
      glutPostRedisplay();
   }
   else if (GLUTbutton[1]) {
      viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
      glutPostRedisplay();
   }
   else if (GLUTbutton[2]) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
         viewer->TranslateWorld(1.0, origin, x, y, dx, dy);
      }
      else {
         R3Vector towards = viewer->Camera().Towards();
         RNDimension dim = towards.MaxDimension();
         R3Vector normal = R3xyz_triad[dim];
         R3Plane plane(selected_position, normal);
         R3Ray viewer_ray = viewer->WorldRay(x, y);
         R3Point intersection;
         if (R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

         // confine point by x dimension
         if (selected_position.X() < world_box.XMin())
            selected_position.SetX(world_box.XMin());
         else if (selected_position.X() >= world_box.XMax())
            selected_position.SetX(world_box.XMax());
         // confine point by y dimension
         if (selected_position.Y() < world_box.YMin())
            selected_position.SetY(world_box.YMin());
         else if (selected_position.Y() >= world_box.YMax())
            selected_position.SetY(world_box.YMax());
         // confine point by z dimension
         if (selected_position.Z() < world_box.ZMin())
            selected_position.SetZ(world_box.ZMin());
         else if (selected_position.Z() >= world_box.ZMax())
            selected_position.SetZ(world_box.ZMax());

         // update selected voxel
         R3Point voxel = selected_position;
         voxel.InverseTransform(nd->Transformation());
         selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
         selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
         selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);

         UpdatePath();
      }
      glutPostRedisplay();
   }
}



void GLUTMotion(int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // different motion clicks for projection view
   if (projection) GLUTMotion2D(x, y);
   else GLUTMotion3D(x, y);

   // redisplay if a mouse was down
   if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;
}



void GLUTMouse2D(int button, int state, int x, int y)
{
   if (button == 2) {
      if (state == GLUT_DOWN) {
         // query 
         R2Point world_point = R2Point(x, y);
         selected_slice_position = ConvertWorldToGrid(world_point);

         // get window dimensions
         int window_xdim = (projection_dim + 1) % 3;
         int window_ydim = (projection_dim + 2) % 3;

         // get the grid x and y coordinates
         int gridx = (int)(selected_slice_position.X() + 0.5);
         int gridy = (int)(selected_slice_position.Y() + 0.5);

         // cannot select a point outside of the image
         if (gridx < 0 || gridx >= nd->Resolution(window_xdim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
         else if (gridy < 0 || gridy >= nd->Resolution(window_ydim)) selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);

         // move voxel
         R2Point grid_point = ConvertWorldToGrid(world_point);

         if (grid_point.X() < 0) grid_point.SetCoord(RN_X, 0);
         if (grid_point.Y() < 0) grid_point.SetCoord(RN_Y, 0);
         if (grid_point.X() > nd->Resolution(window_xdim) - 1) grid_point.SetCoord(RN_X, nd->Resolution(window_xdim) - 1);
         if (grid_point.Y() > nd->Resolution(window_ydim) - 1) grid_point.SetCoord(RN_Y, nd->Resolution(window_ydim) - 1);

         R3Point intersection;
         if (projection_dim == RN_X) intersection = R3Point(selected_slice_index, grid_point.X(), grid_point.Y());
         else if (projection_dim == RN_Y) intersection = R3Point(grid_point.Y(), selected_slice_index, grid_point.X());
         else intersection = R3Point(grid_point.X(), grid_point.Y(), selected_slice_index);

         // update the global selected position
         selected_position = intersection;
         selected_position.Transform(nd->Transformation());

         // update the selected voxel
         selected_voxel[RN_X] = (int)(intersection.X() + 0.5);
         selected_voxel[RN_Y] = (int)(intersection.Y() + 0.5);
         selected_voxel[RN_Z] = (int)(intersection.Z() + 0.5);

         UpdatePath();
      }
      else {
         selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
      }
   }
}



void GLUTMouse3D(int button, int state, int x, int y)
{
   if (button == 2) {
      if (state == GLUT_DOWN) {
         R3Vector towards = viewer->Camera().Towards();
         RNDimension dim = towards.MaxDimension();
         R3Vector normal = R3xyz_triad[dim];
         R3Plane plane(selected_position, normal);
         R3Ray viewer_ray = viewer->WorldRay(x, y);
         R3Point intersection;
         if (R3Intersects(viewer_ray, plane, &intersection)) selected_position = intersection;

         // confine point by x dimension
         if (selected_position.X() < world_box.XMin())
            selected_position.SetX(world_box.XMin());
         else if (selected_position.X() >= world_box.XMax())
            selected_position.SetX(world_box.XMax());
         // confine point by y dimension
         if (selected_position.Y() < world_box.YMin())
            selected_position.SetY(world_box.YMin());
         else if (selected_position.Y() >= world_box.YMax())
            selected_position.SetY(world_box.YMax());
         // confine point by z dimension
         if (selected_position.Z() < world_box.ZMin())
            selected_position.SetZ(world_box.ZMin());
         else if (selected_position.Z() >= world_box.ZMax())
            selected_position.SetZ(world_box.ZMax());

         // update selected voxel
         R3Point voxel = selected_position;
         voxel.InverseTransform(nd->Transformation());
         selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
         selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
         selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);

         UpdatePath();
      }
   }
}



void GLUTMouse(int button, int state, int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // process mouse button event
   if (projection) GLUTMouse2D(button, state, x, y);
   else GLUTMouse3D(button, state, x, y);

   // remember button state
   int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
   GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

   // remember modifiers
   GLUTmodifiers = glutGetModifiers();

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;

   // redraw
   glutPostRedisplay();
}



void GLUTSpecial(int key, int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;

   // remember modifiers
   GLUTmodifiers = glutGetModifiers();

   switch (key) {
   case GLUT_KEY_UP: {
      selected_slice_index++;
      if (selected_slice_index > nd->Resolution(projection_dim) - 1) selected_slice_index = nd->Resolution(projection_dim) - 1;
      char name[128];
      sprintf(name, "Current Slice: %d\n", selected_slice_index);
      glutSetWindowTitle(name);
      UpdateSlices();
      break; }

   case GLUT_KEY_DOWN: {
      selected_slice_index--;
      if (selected_slice_index < 0) selected_slice_index = 0;
      char name[128];
      sprintf(name, "Current Slice: %d\n", selected_slice_index);
      glutSetWindowTitle(name);
      UpdateSlices();
      break; }
   }

   // redraw
   glutPostRedisplay();
}



void GLUTKeyboard2D(unsigned char key, int x, int y)
{
   switch (key) {
   case 'A':
   case 'a': {
      show_projection_affinity = 1;
      show_projection_image = 0;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 0;
      show_projection_predictions = 0;
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
      break; }

   case 'I':
   case 'i': {
      show_projection_affinity = 0;
      show_projection_image = 1;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 0;
      show_projection_predictions = 0;
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'P':
   case 'p': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 0;
      show_projection_predictions = 1;
      selected_slice = predictions_selected_slice;
      selected_slice_range = predictions_range;
      break; }

   case 'S':
   case 's': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_machine_labels = 1;
      show_projection_human_labels = 0;
      show_projection_predictions = 0;
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
      break; }

   case 'T':
   case 't': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 1;
      show_projection_predictions = 0;
      selected_slice = human_labels_selected_slice;
      selected_slice_range = human_labels_range;
      break; }
   }
}



void GLUTKeyboard3D(unsigned char key, int x, int y)
{
   switch (key) {
   case 'A':
   case 'a': {
      show_dimension_affinity = 1;
      show_dimension_image = 0;
      show_dimension_machine_labels = 0;
      show_dimension_human_labels = 0;
      show_dimension_predictions = 0;
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
      break; }

   case 'B':
   case 'b': {
      show_bbox = 1 - show_bbox;
      break; }

   case 'I':
   case 'i': {
      show_dimension_affinity = 0;
      show_dimension_image = 1;
      show_dimension_machine_labels = 0;
      show_dimension_human_labels = 0;
      show_dimension_predictions = 0;
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'P':
   case 'p': {
      show_dimension_affinity = 0;
      show_dimension_image = 0;
      show_dimension_machine_labels = 0;
      show_dimension_human_labels = 0;
      show_dimension_predictions = 1;
      selected_slice = predictions_selected_slice;
      selected_slice_range = predictions_range;
      break; }

   case 'S':
   case 's': {
      show_dimension_affinity = 0;
      show_dimension_image = 0;
      show_dimension_machine_labels = 1;
      show_dimension_human_labels = 0;
      show_dimension_predictions = 0;
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
      break; }

   case 'T':
   case 't': {
      show_dimension_affinity = 0;
      show_dimension_image = 0;
      show_dimension_machine_labels = 0;
      show_dimension_human_labels = 1;
      show_dimension_predictions = 0;
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
      break; }

   case 'W':
   case 'w': {
      show_slice = 1 - show_slice;
      break; }
   }
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
   // invert y coordinate
   y = GLUTwindow_height - y;

   // remember mouse position
   GLUTmouse[0] = x;
   GLUTmouse[1] = y;

   // remember modifiers
   GLUTmodifiers = glutGetModifiers();

   // different keys based on projection or normal
   if (projection) {
      GLUTKeyboard2D(key, x, y);
   }
   else {
      GLUTKeyboard3D(key, x, y);
   }


   // keys regardless of projection status
   switch (key) {
   case '1': {
      affinity_dim = RN_X;
      break; }

   case '2': {
      affinity_dim = RN_Y;
      break; }

   case '3': {
      affinity_dim = RN_Z;
      break; }

   case 'R':
   case 'r': {
      show_prediction_results = 1 - show_prediction_results;
      break; }

   case 'X':
   case 'x': {
      projection_dim = RN_X;
      break; }

   case 'Y':
   case 'y': {
      projection_dim = RN_Y;
      break; }

   case 'Z':
   case 'z': {
      projection_dim = RN_Z;
      break; }

   case ENTER: {
      projection = 1 - projection;
      break; }

   case ESCAPE: {
      GLUTStop();
      break; }
   }

   // redraw
   glutPostRedisplay();
}



void GLUTInit(int *argc, char **argv)
{
   // open window
   glutInit(argc, argv);
   glutInitWindowPosition(10, 10);
   glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   GLUTwindow = glutCreateWindow("Neuron Visualizer");

   // initialize background color
   glClearColor(background_color[0], background_color[1], background_color[2], 1.0);

   // initialize lights
   static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
   glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
   static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
   glEnable(GL_LIGHT0);
   static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
   glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
   glEnable(GL_LIGHT1);
   glEnable(GL_NORMALIZE);
   glEnable(GL_LIGHTING);

   // initialize graphics mode
   glEnable(GL_DEPTH_TEST);

   // initialize GLUT callback functions
   glutDisplayFunc(GLUTRedraw);
   glutReshapeFunc(GLUTResize);
   glutKeyboardFunc(GLUTKeyboard);
   glutSpecialFunc(GLUTSpecial);
   glutMouseFunc(GLUTMouse);
   glutMotionFunc(GLUTMotion);

   // initialize font
#if (RN_OS == RN_WINDOWSNT)
   int font = glGenLists(256);
   wglUseFontBitmaps(wglGetCurrentDC(), 0, 256, font);
   glListBase(font);
#endif
}



void GLUTMainLoop(void)
{
   // run main loop -- never returns
   glutMainLoop();
}



////////////////////////////////////////////////////////////////////////
// Viewer functions
////////////////////////////////////////////////////////////////////////

static R3Viewer *CreateViewer(void)
{
   // start statistics
   RNTime start_time;
   start_time.Read();


   if (world_box.IsEmpty()) RNAbort("Error in CreateViewer - box is empty");
   RNLength r = world_box.DiagonalRadius();
   if (r < 0 || RNIsInfinite(r)) RNAbort("Error in CreateViewer - r must be positive finite");

   // set up camera view looking down the z axis
   static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.5);
   static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
   R3Point initial_camera_origin = world_box.Centroid() - initial_camera_towards * 2.5 * r;
   R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
   R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
   R3Viewer *viewer = new R3Viewer(camera, viewport);

   // print statistics
   if (print_verbose) {
      printf("Created viewer ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  Origin = %g %g %g\n", camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z());
      printf("  Towards = %g %g %g\n", camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z());
      printf("  Up = %g %g %g\n", camera.Up().X(), camera.Up().Y(), camera.Up().Z());
      printf("  Fov = %g %g\n", camera.XFOV(), camera.YFOV());
      printf("  Near = %g\n", camera.Near());
      printf("  Far = %g\n", camera.Far());
      fflush(stdout);
   }

   // return viewer
   return viewer;
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
         else if (!strcmp(*argv, "-debug")) print_debug = 1;
         else if (!strcmp(*argv, "-boundarymerge")) extension = "boundarymerge";
         else if (!strcmp(*argv, "-ordermerge")) extension = "ordermerge";
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
         else if (!extension) extension = *argv;
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      argv++; argc--;
   }

   // error if there is no input name
   if (!input_name) { fprintf(stderr, "Need to supply a neuron data file\n"); return 0; }

   // return success
   return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // parse program arguments
   if (!ParseArgs(argc, argv)) exit(-1);

   // read in the neuron data
   if (!ReadData()) exit(-1);

   char root_filename[4096];
   strncpy(root_filename, nd->Filename(), 4096);
   char *extp = strrchr(root_filename, '.');
   *extp = '\0';

   // read in distance file
   char distances_filename[4096];
   sprintf(distances_filename, "%s/%s_dijkstra_distances_boundary.distance", distances_directory, root_filename);

   FILE *distances_fp = fopen(distances_filename, "rb");
   if (!distances_fp) { fprintf(stderr, "Failed to read %s\n", distances_filename); return 0; }

   int distances_nvoxels;
   fread(&distances_nvoxels, sizeof(int), 1, distances_fp);
   rn_assertion(distances_nvoxels == nd->NVoxels());
   distances = new float[nd->NVoxels()];
   fread(distances, sizeof(float), nd->NVoxels(), distances_fp);

   // close file
   fclose(distances_fp);

   // read in the prev node file
   char prev_node_filename[4096];
   sprintf(prev_node_filename, "%s/%s_dijkstra_prev_node_boundary.distance", distances_directory, root_filename);

   FILE *prev_node_fp = fopen(prev_node_filename, "rb");
   if (!prev_node_fp) { fprintf(stderr, "Failed to read %s\n", prev_node_filename); return 0; }

   int prev_node_nvoxels;
   fread(&prev_node_nvoxels, sizeof(int), 1, prev_node_fp);
   rn_assertion(prev_node_nvoxels == nd->NVoxels());
   prev_node = new int[nd->NVoxels()];
   fread(prev_node, sizeof(int), nd->NVoxels(), prev_node_fp);

   // close file
   fclose(prev_node_fp);

   char meta_root_filename[4096];
   sprintf(meta_root_filename, "%s/%s/%s", output_directory, extension, root_filename);

   // read the prediction grid
   predictions_grid = RNReadNeuronMetaRawFile(meta_root_filename);
   if (!predictions_grid) { fprintf(stderr, "Failed to read %s\n", meta_root_filename); exit(-1); }

   // read in voxels
   if (print_verbose) { printf("Reading voxels..."); fflush(stdout); }
   RNTime voxel_time;
   voxel_time.Read();
   nd->ReadVoxels();
   if (print_verbose) { printf("done in %0.2f seconds.\n", voxel_time.Elapsed()); }

   // go through all predictions to find ones that don't exit
   int npredictions = (int)(predictions_grid->Maximum() + 0.5);
   RNBoolean *prediction_boundaries = new RNBoolean[npredictions];
   for (int ip = 0; ip < npredictions; ++ip)
      prediction_boundaries[ip] = FALSE;

   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (nd->Voxel(iv)->IsOnBoundary()) {
         int prediction_value = (int)(predictions_grid->GridValue(iv) + 0.5);
         prediction_boundaries[prediction_value] = TRUE;
      }
   }

   // create voxel hash
   if (print_verbose) { printf("Creating voxel hash...\n  "); fflush(stdout); }
   voxel_hash = new int[nd->NVoxels()];
   for (int iv = 0; iv < nd->NVoxels(); ++iv)
      voxel_hash[iv] = 0;

   // go through every voxel in this cellular
   for (int iv = 0; iv < nd->NVoxels(); ++iv) {
      if (print_verbose) RNProgressBar(iv, nd->NVoxels());

      NeuronVoxel *voxel = nd->Voxel(iv);
      if (voxel->IsExtracellular()) continue;

      NeuronSupervoxel *supervoxel = voxel->Supervoxel();
      int prediction_value = (int)(predictions_grid->GridValue(iv) + 0.5);

      // go through all the neighbors
      for (int in = 0; in < voxel->NNeighbors(); ++in) {
         // ignore the z direction
         if (in == 4 || in == 5) continue;
         NeuronVoxel *neighbor = voxel->Neighbor(in);
         if (!neighbor) continue;
         if (neighbor->IsExtracellular()) continue;

         NeuronSupervoxel *neighbor_supervoxel = neighbor->Supervoxel();
         if (supervoxel == neighbor_supervoxel) continue;

         int neighbor_prediction_value = (int)(predictions_grid->GridValue(neighbor->DataIndex()) + 0.5);

         int label = (supervoxel->MajorityHumanLabel() == neighbor_supervoxel->MajorityHumanLabel());
         int prediction = (neighbor_prediction_value == prediction_value);

         if (label && prediction) voxel_hash[iv] = 1;
         else if (!label && !prediction) voxel_hash[iv] = 2;
         else if (label && !prediction) voxel_hash[iv] = 3;
         else if (!label && prediction) voxel_hash[iv] = 4;
      }
   }

   nd->ReleaseVoxels();
   if (print_verbose) printf("\ndone!\n");

   delete[] prediction_boundaries;

   // set world box
   world_box = nd->WorldBox();
   selected_position = world_box.Centroid();

   // initialize selected voxel
   R3Point voxel = selected_position;
   voxel.InverseTransform(nd->Transformation());
   selected_voxel[RN_X] = (int)(voxel.X() + 0.5);
   selected_voxel[RN_Y] = (int)(voxel.Y() + 0.5);
   selected_voxel[RN_Z] = (int)(voxel.Z() + 0.5);
   selected_slice_index = selected_voxel[projection_dim];

   // set projection scale variables
   for (int dim = 0; dim <= 2; ++dim) {
      RNScalar ratio = GLUTwindow_height / (nd->Scale(dim) * nd->Resolution(dim));
      if (ratio < projection_scale) projection_scale = ratio;
   }

   /////////////////////////////////
   //// Read in the voxel files ////
   /////////////////////////////////

   if (!ReadVoxelFiles()) exit(-1);


   /////////////////////////////
   //// Read in the R3Grids ////
   /////////////////////////////

   // create viewer
   viewer = CreateViewer();
   if (!viewer) exit(-1);

   // initialize GLUT
   GLUTInit(&argc, argv);

   // run GLUT interface
   GLUTMainLoop();

   // return success
   return 1;
}