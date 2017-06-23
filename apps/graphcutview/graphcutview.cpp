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




// program variables

static NeuronData *nd = NULL;
static R3Viewer *viewer = NULL;
static R3Point selected_position;
static R3Box world_box;
static int selected_voxel[3];
static char root_filename[4096];



// voxel grids

static R3Grid *affinity_grid[3] = { NULL, NULL, NULL };
static R3Grid *human_labels_grid = NULL;
static R3Grid *image_grid = NULL;
static R3Grid *machine_labels_grid = NULL;



// projection variables

static R2Grid *affinity_selected_slice[3] = { NULL, NULL, NULL };
static R2Grid *human_labels_selected_slice = NULL;
static R2Grid *image_selected_slice = NULL;
static R2Grid *machine_labels_selected_slice = NULL;
static RNInterval affinity_range[3];
static RNInterval human_labels_range;
static RNInterval image_range;
static RNInterval machine_labels_range;
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
static RNScalar human_labels_color[] = { 1, 0, 0 };
static RNScalar machine_labels_color[] = { 0, 0, 0.5 };



// projection color variables

static int color_type = 0; // 0=gray, 1=red-green-blue



// display projection variables

static int show_projection_affinity = 0;
static int show_projection_image = 1;
static int show_projection_human_labels = 0;
static int show_projection_machine_labels = 0;
static int projection_dim = RN_Z;
static int affinity_dim = RN_X;
static int selected_slice_index = 0;



// display dimenson variables

static int show_bbox = 1;
static int show_slice = 0;
static int show_human_labels = 0;
static int show_machine_labels = 0;
static int projection = 0;
static int show_dimension_affinity = 0;
static int show_dimension_image = 1;
static int show_dimension_human_labels = 0;
static int show_dimension_machine_labels = 0;



// index variables

static int human_labels_index = 0;
static int supervoxel_index = 0;



// size variables

static RNScalar voxel_size = 5.0;



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



static void GLUTDrawMenuText(int offset, const char *s) {
   // OpenGL text prologue
   glDisable(GL_LIGHTING);
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D(0, GLUTwindow_width, 0, GLUTwindow_height);
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   RNLoadRgb(RNblack_rgb);
   // draw text string s by 13 * offset from the top
   glRasterPos2i(GLUTwindow_width - 196, GLUTwindow_height - 19 * offset);
   while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));

   // OpenGL text epilogue
   glPopMatrix();
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
   glMatrixMode(GL_MODELVIEW);
}



static void DrawMachineLabel(void)
{
   RNLoadRgb(machine_labels_color);
   // draw the supervoxel
   NeuronSupervoxel *supervoxel = nd->Supervoxel(supervoxel_index);
   supervoxel->Draw();
}



static void DrawHumanLabel(void)
{
   RNLoadRgb(human_labels_color);
   // draw the truth neuron
   NeuronHumanLabel *human_label = nd->HumanLabel(human_labels_index);
   human_label->Draw();
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
            // convert from world to grid coordinates
            R2Point grid_position = R2Point(ix, iy + k);
            R2Point world_position = ConvertGridToWorld(grid_position);

            // get the color for this area
            RNRgb color;
            if (show_projection_human_labels) {
               if (selected_slice->GridValue(ix, iy + k) < 0.5) color = RNred_rgb;
               else color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
            }
            else {
               color = Color(selected_slice->GridValue(ix, iy + k), selected_slice_range);
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

   if (show_human_labels) DrawHumanLabel();

   if (show_slice) DrawSlice();

   if (show_machine_labels) DrawMachineLabel();

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

   // draw menu
   int offset = 0;
   char menu_text[4096];
   sprintf(menu_text, "Human Label Index: %d", human_labels_index);
   GLUTDrawMenuText(++offset, menu_text);
   sprintf(menu_text, "Supervoxel Index: %d", supervoxel_index);
   GLUTDrawMenuText(++offset, menu_text);

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
   }
}



void GLUTMotion3D(int x, int y)
{
   // compute mouse movement
   int dx = x - GLUTmouse[0];
   int dy = y - GLUTmouse[1];

   // world in hand navigation
   R3Point origin = world_box.Centroid();
   if (GLUTbutton[0]) viewer->RotateWorld(1.0, origin, x, y, dx, dy);
   else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, origin, x, y, dx, dy);
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
      }
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
   if (button == 0) {
      if (state == GLUT_DOWN) {
         // query 
         R2Point world_point = R2Point(x, y);
         selected_slice_position = ConvertGridToWorld(world_point);

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
      else {
         selected_slice_position.Reset(RN_UNKNOWN, RN_UNKNOWN);
         glutPostRedisplay();
      }
   }
}



void GLUTMouse3D(int button, int state, int x, int y)
{
   if (button == 2) {
      if (glutGetModifiers() != GLUT_ACTIVE_SHIFT) {
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
         }
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
      UpdateSlices();
      break; }

   case GLUT_KEY_DOWN: {
      selected_slice_index--;
      if (selected_slice_index < 0) selected_slice_index = 0;
      UpdateSlices();
      break; }

   case GLUT_KEY_LEFT: {
      if (GLUTmodifiers == GLUT_ACTIVE_SHIFT) {
         if (supervoxel_index == 0) break;
         NeuronSupervoxel *supervoxel = nd->Supervoxel(supervoxel_index);
         supervoxel->ReleaseVoxels();
         supervoxel_index--;
         supervoxel = nd->Supervoxel(supervoxel_index);
         supervoxel->ReadVoxels();
      }
      else {
         if (human_labels_index == 0) break;
         NeuronHumanLabel *human_label = nd->HumanLabel(human_labels_index);
         human_label->ReleaseVoxels();
         human_labels_index--;
         human_label = nd->HumanLabel(human_labels_index);
         human_label->ReadVoxels();
      }
      break; }

   case GLUT_KEY_RIGHT: {
      if (GLUTmodifiers == GLUT_ACTIVE_SHIFT) {
         if (supervoxel_index == nd->NSupervoxels() - 1) break;
         NeuronSupervoxel *supervoxel = nd->Supervoxel(supervoxel_index);
         supervoxel->ReleaseVoxels();
         supervoxel_index++;
         supervoxel = nd->Supervoxel(supervoxel_index);
         supervoxel->ReadVoxels();
      }
      else {
         if (human_labels_index == nd->NHumanLabels() - 1) break;
         NeuronHumanLabel *human_label = nd->HumanLabel(human_labels_index);
         human_label->ReleaseVoxels();
         human_labels_index++;
         human_label = nd->HumanLabel(human_labels_index);
         human_label->ReadVoxels();
      }
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
      selected_slice = affinity_selected_slice[projection_dim];
      selected_slice_range = affinity_range[projection_dim];
      break; }

   case 'I':
   case 'i': {
      show_projection_affinity = 0;
      show_projection_image = 1;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 0;
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'S':
   case 's': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_machine_labels = 1;
      show_projection_human_labels = 0;
      selected_slice = machine_labels_selected_slice;
      selected_slice_range = machine_labels_range;
      break; }

   case 'T':
   case 't': {
      show_projection_affinity = 0;
      show_projection_image = 0;
      show_projection_machine_labels = 0;
      show_projection_human_labels = 1;
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
      selected_slice = image_selected_slice;
      selected_slice_range = image_range;
      break; }

   case 'S':
   case 's': {
      if (GLUTmodifiers == GLUT_ACTIVE_SHIFT) {
         show_machine_labels = 1 - show_machine_labels;
      }
      else {
         show_dimension_affinity = 0;
         show_dimension_image = 0;
         show_dimension_machine_labels = 1;
         show_dimension_human_labels = 0;
         selected_slice = machine_labels_selected_slice;
         selected_slice_range = machine_labels_range;
      }
      break; }

   case 'T':
   case 't': {
      if (GLUTmodifiers == GLUT_ACTIVE_SHIFT) {
         show_human_labels = 1 - show_human_labels;
      }
      else {
         show_dimension_affinity = 0;
         show_dimension_image = 0;
         show_dimension_machine_labels = 0;
         show_dimension_human_labels = 1;
         selected_slice = machine_labels_selected_slice;
         selected_slice_range = machine_labels_range;
      }
      break; }

   case 'W':
   case 'w': {
      show_slice = 1 - show_slice;
      break; }

   case SPACEBAR: {
      R3Ray ray = viewer->WorldRay(x, y);
      RNScalar closest = FLT_MAX;
      for (int is = 0; is < nd->NSupervoxels(); ++is) {
         NeuronSupervoxel *supervoxel = nd->Supervoxel(is);
         RNScalar dist = R3Distance(ray, supervoxel->Centroid());
         if (dist < closest) {
            closest = dist;
            supervoxel_index = is;
         }
      }
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
         else { fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
      }
      else {
         if (!input_name) input_name = *argv;
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

   // create output filename prefix
   strncpy(root_filename, input_name, 4096);
   char *endp = strrchr(root_filename, '.');
   if (endp) *endp = '\0';

   // create viewer
   viewer = CreateViewer();
   if (!viewer) exit(-1);

   // read in the necessary voxels
   NeuronHumanLabel *human_label = nd->HumanLabel(human_labels_index);
   human_label->ReadVoxels();
   NeuronSupervoxel *supervoxel = nd->Supervoxel(supervoxel_index);
   supervoxel->ReadVoxels();

   // initialize GLUT
   GLUTInit(&argc, argv);

   // run GLUT interface
   GLUTMainLoop();

   // return success
   return 1;
}