// Source file for histogram



#include "RNHistogram.h"


RNHistogram::
RNHistogram(void) :
vector(),
sorted(FALSE)
{
}



RNHistogram::
~RNHistogram(void)
{
   vector.Empty();
}



void RNHistogram::
AddScalarQuantity(RNScalar value)
{
   sorted = FALSE;
   // add this value to the histogram
   vector.Insert(value);
}



int RNScalarSort(const RNScalar data1, const RNScalar data2)
{
   return (data1 < data2);
}



void RNHistogram::
SortVector(void)
{
   // sort either ascending or descending
   vector.QuickSort(RNScalarSort);
   sorted = TRUE;
}



int RNHistogram::
Histogram(int nbins, int *bins, RNScalar *boundaries, const char *output_filename)
{
   // create a histogram with nbins
   if (nbins < 0) return 0;
   if (!bins) return 0;
   if (!boundaries) return 0;

   // sort vector so that head and tail produce min and max
   if (!sorted) SortVector();

   // get minimum and maximum values
   RNScalar minimum_value = vector.Head();
   RNScalar maximum_value = vector.Tail();

   // get range and bin size
   RNScalar range = maximum_value - minimum_value;
   /* add a small constant to avoid double precision errors */
   RNScalar bin_size = range / nbins + 10e-6;

   // make sure bins and boundaries have defaults
   for (int ib = 0; ib < nbins; ++ib) {
      bins[ib] = 0;
      boundaries[ib] = minimum_value + bin_size * ib;
   }

   // create histogram
   for (int iv = 0; iv < vector.Size(); ++iv) {
      RNScalar value = vector[iv];
      int bin = (int)((value - minimum_value) / bin_size);
      rn_assertion((0 <= bin) && (bin < nbins));
      bins[bin]++;
   }

   // write results
   if (output_filename) {
      // open file
      FILE *fp = fopen(output_filename, "w");
      if (!fp) { fprintf(stderr, "Failed to open %s\n", output_filename); return 0; }

      for (int ib = 0; ib < nbins; ++ib) {
         fprintf(fp, "%f,%d\n", boundaries[ib], bins[ib]);
      }

      // close file 
      fclose(fp);
   }

   // return success
   return 1;
}

