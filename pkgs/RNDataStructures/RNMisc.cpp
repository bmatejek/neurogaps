#include "RNDataStructures.h"



void
RNProgressBar(int index, int nindices)
{
   rn_assertion((0 <= index) && (index < nindices));

   if (index == 0) { printf("0%%"); fflush(stdout); return; }
   else if (index == nindices - 1) { printf("100%%"); fflush(stdout); return; }

   RNScalar percentage = (100 * (RNScalar)index) / nindices;
   RNScalar prev_percentage = (100 * (RNScalar)(index - 1)) / nindices;

   int percent = (int)(percentage);
   int prev_percent = (int)(prev_percentage);

   if (percent == prev_percent) return;

   if (percent % 10 == 0) {
      printf("%d%%", percent);
      fflush(stdout);
   }
   else if (percent % 2 == 0) {
      printf(".");
      fflush(stdout);
   }
}



void
RNDeflateIntegerArray(int *entries, int nentries)
{
   rn_assertion(entries != NULL);

   // find the maximum entry in entries
   int max_entry = -1;
   for (int ie = 0; ie < nentries; ++ie) {
      rn_assertion(entries[ie] >= 0);
      if (max_entry < entries[ie]) {
         max_entry = entries[ie];
      }
   }

   // find the entries that exist
   RNBoolean *entry_exists = new RNBoolean[max_entry + 1];
   for (int ie = 0; ie < max_entry + 1; ++ie) {
      entry_exists[ie] = FALSE;
   }

   // see if this entry exists
   for (int ie = 0; ie < nentries; ++ie) {
      int entry = entries[ie];
      entry_exists[entry] = TRUE;
   }

   // create mapping
   int *mapping = new int[max_entry + 1];
   int num_skipped = 0;
   for (int ie = 0; ie < max_entry + 1; ++ie) {
      if (entry_exists[ie]) mapping[ie] = ie - num_skipped;
      else { num_skipped++; mapping[ie] = -1; }
   }

   // reset entries' values
   for (int ie = 0; ie < nentries; ++ie) {
      entries[ie] = mapping[entries[ie]];
   }

   // free memory
   delete[] entry_exists;
   delete[] mapping;
}



void 
RNBestFitLine(RNScalar *x, RNScalar *y, int n, RNScalar& alpha, RNScalar& beta, RNScalar& rsquared)
{
   // just checking
   rn_assertion(x != NULL);
   rn_assertion(y != NULL);
   rn_assertion(n > 0);

   // calcluate the average of x and y
   RNScalar avgx = 0.0;
   RNScalar avgy = 0.0;
   for (int i = 0; i < n; ++i) {
      avgx += x[i];
      avgy += y[i];
   }
   avgx /= n;
   avgy /= n;
   
   // calculate beta
   RNScalar beta_numerator = 0.0;
   RNScalar beta_denominator = 0.0;
   for (int i = 0; i < n; ++i) {
      beta_numerator += (x[i] - avgx) * (y[i] - avgy);
      beta_denominator += (x[i] - avgx) * (x[i] - avgx);
   }

   beta = beta_numerator / beta_denominator;
   
   // calculate alpha from beta
   alpha = avgy - beta * avgx;
   
   // calculate the coefficient of determination
   RNScalar avgxy = 0.0;
   RNScalar avgxsquared = 0.0;
   RNScalar avgysquared = 0.0;
   for (int i = 0; i < n; ++i) {
      avgxy += x[i] * y[i];
      avgxsquared += x[i] * x[i];
      avgysquared += y[i] * y[i];
   }
   avgxy /= n;
   avgxsquared /= n;
   avgysquared /= n;

   // rxy
   rsquared = (avgxy - avgx * avgy) / sqrt((avgxsquared - avgx * avgx) * (avgysquared - avgy * avgy));
   
   // coefficient of determination
   rsquared = rsquared * rsquared;
}
