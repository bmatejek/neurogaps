// Source file for RNPyPlot class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RNDataStructures/RNDataStructures.h"



////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

RNPyPlot::
RNPyPlot(void) :
title(),
xlabel(),
ylabel(),
legend(),
extrema_type(NO_EXTREMA)
{
}



RNPyPlot::
~RNPyPlot(void)
{
}



RNPyHistogram::
RNPyHistogram(void) :
points(),
plot_type(HISTOGRAM),
nbins(100)
{
}



RNPyHistogram::
~RNPyHistogram(void)
{
}



RNPyPointPlot::
RNPyPointPlot(void) :
points(),
plot_type(LINE)
{
   // default is to ignore bounds
   xlim[0] = 1;
   xlim[1] = 0;
   ylim[0] = 1;
   ylim[1] = 0;
}



RNPyPointPlot::
~RNPyPointPlot(void)
{
}



RNPyImage::
RNPyImage(void) :
output_directory(),
plots()
{
   strncpy(extension, "png", 8);
}



RNPyImage::
~RNPyImage(void)
{
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void RNPyPlot::
SetTitle(const char *title)
{
   // set the title
   strncpy(this->title, title, 128);
}



void RNPyPlot::
SetXLabel(const char *xlabel)
{
   // set x label
   strncpy(this->xlabel, xlabel, 128);
}



void RNPyPlot::
SetYLabel(const char *ylabel)
{
   // set y label
   strncpy(this->ylabel, ylabel, 128);
}



void RNPyPlot::
SetLegend(const char *legend)
{
   // set the legend
   strncpy(this->legend, legend, 128);
}



void RNPyPlot::
SetXAxisMin(RNScalar xaxis_min)
{
   // set the xaxis min
   xlim[0] = xaxis_min;
}



void RNPyPlot::
SetXAxisMax(RNScalar xaxis_max)
{
   // set the xaxis max
   xlim[1] = xaxis_max;
}



void RNPyPlot::
SetYAxisMin(RNScalar yaxis_min)
{
   // set the yaxis min
   ylim[0] = yaxis_min;
}



void RNPyPlot::
SetYAxisMax(RNScalar yaxis_max)
{
   // set the yaxis max
   ylim[1] = yaxis_max;
}



void RNPyPlot::
SetExtremaType(enum DRAW_EXTREMA_TYPE extrema_type)
{
   // set extrema draw type
   this->extrema_type = extrema_type;
}



void RNPyHistogram::
InsertPoint(RNScalar point)
{
   // add this point to vector
   points.push_back(point);
}



void RNPyHistogram::
SetNBins(int nbins)
{
   // set the number of bins for the image
   this->nbins = nbins;
}



void RNPyPointPlot::
InsertPoint(R2Point point)
{
   // add this point to vector
   points.push_back(point);
}



void RNPyImage::
SetOutputDirectory(const char *output_directory)
{
   // set the output directory
   strncpy(this->output_directory, output_directory, 128);
}



void RNPyImage::
SetExtension(const char *extension)
{
   // set the extension 
   strncpy(this->extension, extension, 8);
}



void RNPyImage::
InsertPyPlot(RNPyPlot *plot)
{
   // insert this py plot into vector
   plots.push_back(plot);
}



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

int RNPyPointPlot::
WritePlot(FILE *fp) const
{
   rn_assertion(fp != NULL);

   // write the title for this plot
   fwrite(title, sizeof(char), 128, fp);

   // write the x axis label
   fwrite(xlabel, sizeof(char), 128, fp);

   // write the y axis label
   fwrite(ylabel, sizeof(char), 128, fp);

   // write the legend
   fwrite(legend, sizeof(char), 128, fp);

   // write the type of plot
   fwrite(&plot_type, sizeof(int), 1, fp);

   // write the extrema draw type
   fwrite(&extrema_type, sizeof(int), 1, fp);

   // write the axis limis
   fwrite(&(xlim[0]), sizeof(RNScalar), 1, fp);
   fwrite(&(xlim[1]), sizeof(RNScalar), 1, fp);
   fwrite(&(ylim[0]), sizeof(RNScalar), 1, fp);
   fwrite(&(ylim[1]), sizeof(RNScalar), 1, fp);

   // write the number of points
   int npoints = points.size();
   fwrite(&npoints, sizeof(int), 1, fp);

   for (int ip = 0; ip < npoints; ++ip) {
      RNScalar xcoord = points[ip].X();
      RNScalar ycoord = points[ip].Y();
      fwrite(&xcoord, sizeof(RNScalar), 1, fp);
      fwrite(&ycoord, sizeof(RNScalar), 1, fp);
   }

   // return success 
   return 1;
}



int RNPyHistogram::
WritePlot(FILE *fp) const
{
   rn_assertion(fp != NULL);

   // write the title for this plot
   fwrite(title, sizeof(char), 128, fp);

   // write the x axis label
   fwrite(xlabel, sizeof(char), 128, fp);

   // write the y axis label
   fwrite(ylabel, sizeof(char), 128, fp);

   // write the legend
   fwrite(legend, sizeof(char), 128, fp);

   // write the type of plot
   fwrite(&plot_type, sizeof(int), 1, fp);

   // write the extrema draw type
   fwrite(&extrema_type, sizeof(int), 1, fp);

   // write the axis limis
   fwrite(&(xlim[0]), sizeof(RNScalar), 1, fp);
   fwrite(&(xlim[1]), sizeof(RNScalar), 1, fp);
   fwrite(&(ylim[0]), sizeof(RNScalar), 1, fp);
   fwrite(&(ylim[1]), sizeof(RNScalar), 1, fp);

   // write the number of bins
   fwrite(&nbins, sizeof(int), 1, fp);

   // write the number of points
   int npoints = points.size();
   fwrite(&npoints, sizeof(int), 1, fp);

   for (int ip = 0; ip < npoints; ++ip)
      fwrite(&points[ip], sizeof(RNScalar), 1, fp);

   // return success
   return 1;
}



int RNPyImage::
WriteImageFile(const char output_filename[4096]) const
{
   // open file
   FILE *fp = fopen(output_filename, "wb");
   if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return 0; }

   // write the magic keyword
   char magic[16];
   sprintf(magic, "PYIMAGE");
   fwrite(magic, sizeof(char), 16, fp);

   // write the desired extension
   fwrite(extension, sizeof(char), 8, fp);

   // write the output directory
   fwrite(output_directory, sizeof(char), 128, fp);

   // write the number of plots to be expected
   int nplots = plots.size();
   fwrite(&nplots, sizeof(int), 1, fp);

   // write all of the plots
   for (int ip = 0; ip < nplots; ++ip) {
      RNPyPlot *plot = plots[ip];
      plot->WritePlot(fp);
   }

   // close file
   fclose(fp);

   // call the python function to create display
   char function_call[4096];
   sprintf(function_call, "createplot.py %s", output_filename);
   system(function_call);

   // return success
   return 1;
}


