// Include file for the RNPyPlot class

#ifndef __RN_PY_PLOT_H__
#define __RN_PY_PLOY_H__



/* useful constants */

enum PLOT_TYPE {
   LINE,
   HISTOGRAM,
   NTYPES
};



enum DRAW_EXTREMA_TYPE {
   NO_EXTREMA,
   MAX_EXTREMA,
   MIN_EXTREMA,
   NEXTREMA_OPTIONS
};



/////////////////////////////////////////////////////////////////////
// PY PLOT CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNPyPlot {
public:
   // constructors/destructors
   RNPyPlot(void);
   virtual ~RNPyPlot(void);

   // access functions
   const char *Title(void) const;
   const char *XLabel(void) const;
   const char *YLabel(void) const;
   const char *Legend(void) const;
   RNScalar XAxisMin(void) const;
   RNScalar XAxisMax(void) const;
   RNScalar YAxisMin(void) const;
   RNScalar YAxisMax(void) const;
   enum DRAW_EXTREMA_TYPE DrawExtremaType(void) const;
   virtual enum PLOT_TYPE PlotType(void) const;
   virtual int NPoints(void) const;

   // manipulation functions
   void SetTitle(const char *title);
   void SetXLabel(const char *xlabel);
   void SetYLabel(const char *ylabel);
   void SetLegend(const char *legend);
   void SetXAxisMin(RNScalar xaxis_min);
   void SetXAxisMax(RNScalar xaxis_max);
   void SetYAxisMin(RNScalar yaxis_min);
   void SetYAxisMax(RNScalar yaxis_max);
   void SetExtremaType(enum DRAW_EXTREMA_TYPE extrema_type);

protected:
   // I/O functions
   virtual int WritePlot(FILE *fp) const;

protected:
   friend class RNPyImage;
   friend class RNPyHistogram;
   friend class RNPyPointPlot;
   // instance variables
   char title[128];
   char xlabel[128];
   char ylabel[128];
   char legend[128];
   RNScalar xlim[2];
   RNScalar ylim[2];
   enum DRAW_EXTREMA_TYPE extrema_type;
};



/* inline functions */

inline const char *RNPyPlot::
Title(void) const
{
   // return the title
   return title;
}



inline const char *RNPyPlot::
XLabel(void) const
{
   // return the xlabel
   return xlabel;
}



inline const char *RNPyPlot::
YLabel(void) const
{
   // return the ylabel
   return ylabel;
}



inline const char *RNPyPlot::
Legend(void) const
{
   // return the legend title
   return legend;
}



inline RNScalar RNPyPlot::
XAxisMin(void) const
{
   // return min of x axis
   return xlim[0];
}



inline RNScalar RNPyPlot::
XAxisMax(void) const
{
   // return max of x axis
   return xlim[1];
}



inline RNScalar RNPyPlot::
YAxisMin(void) const
{
   // return min of y axis
   return ylim[0];
}



inline RNScalar RNPyPlot::
YAxisMax(void) const
{
   // reutrn max of y axis
   return ylim[1];
}



inline enum DRAW_EXTREMA_TYPE RNPyPlot::
DrawExtremaType(void) const
{
   // return the extrema type
   return extrema_type;
}



inline enum PLOT_TYPE RNPyPlot::
PlotType(void) const
{
   /* overridden */
   return NTYPES;
}



inline int RNPyPlot::
NPoints(void) const
{
   /* overridden */
   return -1;
}



inline int RNPyPlot::
WritePlot(FILE *fp) const
{
   /* overridden */
   return 0;
}



/////////////////////////////////////////////////////////////////////
// PY HISTOGRAM CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNPyHistogram : public RNPyPlot {
public:
   // constructors/destructors
   RNPyHistogram(void);
   virtual ~RNPyHistogram(void);

   // access functions
   virtual enum PLOT_TYPE PlotType(void) const;
   virtual int NPoints(void) const;
   RNScalar Point(int point_index) const;
   int NBins(void) const;

   // manipulation functions
   void InsertPoint(RNScalar point);
   void SetNBins(int nbins);

protected:
   // I/O functions
   virtual int WritePlot(FILE *fp) const;


private:
   // instance variables
   std::vector<RNScalar> points;
   enum PLOT_TYPE plot_type;
   int nbins;
};



/* inline functions */

inline enum PLOT_TYPE RNPyHistogram::
PlotType(void) const
{
   // return the plot type
   return plot_type;
}



inline int RNPyHistogram::
NPoints(void) const
{
   // return the number of points
   return points.size();
}



inline RNScalar RNPyHistogram::
Point(int point_index) const
{
   rn_assertion((0 <= point_index) && (point_index < (int)points.size()));
   // return this point index
   return points[point_index];
}



inline int RNPyHistogram::
NBins(void) const
{
   // return the number of bins
   return nbins;
}



/////////////////////////////////////////////////////////////////////
// PY POINT PLOT CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNPyPointPlot : public RNPyPlot {
public:
   // constructors/destructors
   RNPyPointPlot(void);
   virtual ~RNPyPointPlot(void);

   // access functions
   virtual enum PLOT_TYPE PlotType(void) const;
   virtual int NPoints(void) const;
   const R2Point& Point(int point_index) const;

   // manipulation functions
   void InsertPoint(R2Point point);

protected:
   // I/O functions
   virtual int WritePlot(FILE *fp) const;


private:
   // instance variables
   std::vector<R2Point> points;
   enum PLOT_TYPE plot_type;
};



/* inline functions */

inline enum PLOT_TYPE RNPyPointPlot::
PlotType(void) const
{
   // return the plot type
   return plot_type;
}



inline int RNPyPointPlot::
NPoints(void) const
{
   // return the number of points
   return points.size();
}



inline const R2Point& RNPyPointPlot::
Point(int point_index) const
{
   rn_assertion((0 <= point_index) && (point_index < (int)points.size()));
   // return the point
   return points[point_index];
}



/////////////////////////////////////////////////////////////////////
// PY IMAGE CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

class RNPyImage {
public:
   // constructors/destructors
   RNPyImage(void);
   ~RNPyImage(void);

   // access functions
   const char *OutputDirectory(void) const;
   const char *Extension(void) const;

   int NPlots(void) const;
   const RNPyPlot& Plot(int plot_index) const;

   // manipulation functions
   void SetOutputDirectory(const char *output_directory);
   void SetExtension(const char *extension);
   void InsertPyPlot(RNPyPlot *plot);

   // I/O functions
   int WriteImageFile(const char output_filename[4096]) const;


private:
   // instance variables
   char output_directory[128];
   char extension[8];
   std::vector<RNPyPlot *> plots;
};



/* inline functions */

inline const char *RNPyImage::
OutputDirectory(void) const
{
   // return the output directory
   return output_directory;
}



inline const char *RNPyImage::
Extension(void) const
{
   // return the extension
   return extension;
}



inline int RNPyImage::
NPlots(void) const
{
   // return the number of plots
   return plots.size();
}



inline const RNPyPlot &RNPyImage::
Plot(int plot_index) const
{
   rn_assertion((0 <= plot_index) && (plot_index < (int)plots.size()));
   // return this plot
   return *plots[plot_index];
}



#endif