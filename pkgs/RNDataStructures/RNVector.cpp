/* Source file for the vector class */

#ifndef __RN__VECTOR__C__
#define __RN__VECTOR__C__


/* Include files */

#include "RNDataStructures/RNDataStructures.h"



/* private constants */

static const int RN_VECTOR_MIN_ALLOCATED = 128;



/* public functions */

template <typename Type>
RNVector<Type>::RNVector(void) :
nallocated(RN_VECTOR_MIN_ALLOCATED),
nentries(0),
sorted(FALSE)
{
   entries = new Type *[RN_VECTOR_MIN_ALLOCATED];
}



template <typename Type>
RNVector<Type>::~RNVector(void)
{
   if (entries) {
      for (int ie = 0; ie < nentries; ++ie)
         delete entries[ie];
      delete[] entries;
   }
}



template <typename Type>
const RNBoolean RNVector<Type>::
IsEmpty(void) const
{
   // return if vector is empty
   return (nentries == 0);
}



template <typename Type>
const int RNVector<Type>::
NAllocated(void) const
{
   // return the number of allocated spots
   return nallocated;
}



template <typename Type>
const int RNVector<Type>::
NEntries(void) const
{
   // return the number of entries in the vector
   return nentries;
}



template <typename Type>
const int RNVector<Type>::
Size(void) const
{
   // return the size of the vector
   return nentries;
}



template <typename Type>
const Type& RNVector<Type>::
Head(void) const
{
   rn_assertion(nentries);
   // return head
   return *entries[0];
}



template <typename Type>
const Type& RNVector<Type>::
Kth(int k) const
{
   rn_assertion(k < nentries);
   // return kth entry
   return *entries[k];
}



template <typename Type>
const Type& RNVector<Type>::
operator[](int k) const
{
   rn_assertion(k < nentries);
   // return kth entry
   return *entries[k];
}



template <typename Type>
const Type& RNVector<Type>::
Tail(void) const
{
   rn_assertion(nentries);
   // return tail
   return *entries[nentries - 1];
}



template <typename Type>
void RNVector<Type>::
Insert(Type data)
{
   sorted = FALSE;
   // insert data at end
   InsertTail(data);
}



template <typename Type>
void RNVector<Type>::
InsertHead(Type data)
{
   sorted = FALSE;
   // insert data into vector at head
   InternalInsert(data, 0);
}



template <typename Type>
void RNVector<Type>::
InsertKth(Type data, int k)
{
   sorted = FALSE;
   // insert data into vector in kth position
   InternalInsert(data, k);
}



template <typename Type>
void RNVector<Type>::
InsertTail(Type data)
{
   sorted = FALSE;
   // insert data into vector at tail
   InternalInsert(data, nentries);
}


template <typename Type>
void RNVector<Type>::
Remove(void)
{
   sorted = FALSE;
   // remove tail entry
   RemoveTail();
}



template <typename Type>
void RNVector<Type>::
RemoveHead(void)
{
   sorted = FALSE;
   // remove head entry from array
   InternalRemove(0);
}



template <typename Type>
void RNVector<Type>::
RemoveKth(int k)
{
   sorted = FALSE;
   // remove kth entry
   InternalRemove(k);
}



template <typename Type>
void RNVector<Type>::
RemoveTail(void)
{
   sorted = FALSE;
   // remove tail entry from array
   InternalRemove(nentries - 1);
}



template <typename Type>
void RNVector<Type>::
Empty(RNBoolean deallocate)
{
   // remove all entries from vector
   Truncate(0);

   // deallocate memory
   if (deallocate) {
      if (entries) {
         for (int ie = 0; ie < nentries; ++ie)
            delete entries[ie];
         delete[] entries;
      }
      entries = NULL;
      nallocated = 0;
   }
}



template <typename Type>
void RNVector<Type>::
Truncate(int length)
{
   // remove tail entries from vector
   if (length < nentries) {
      for (int ie = length; ie < nentries; ++ie) {
         delete entries[ie];
      }
      nentries = length;
   }

   Resize(length);
}



template <typename Type>
void RNVector<Type>::
Shift(int delta)
{
   // shift all entries by delta
   Shift(0, 0, delta);
}



template <typename Type>
void RNVector<Type>::
Shift(int start, int length, int delta)
{
   // compute number of entries to shift
   if ((delta < 0) && (start < -delta)) start = -delta;
   int nshift = nentries - start;
   if (delta > 0) nshift -= delta;
   if (nshift <= 0) return;
   if ((length > 0) && (length < nshift)) nshift = length;

   // shift array entries
   if (delta < 0) {
      for (int ie = start; ie < (start + nshift); ++ie) {
         entries[ie + delta] = entries[ie];
      }
   }
   else if (delta > 0) {
      for (int ie = (start + nshift - 1); ie >= start; --ie) {
         entries[ie + delta] = entries[ie];
      }
   }
}



template <typename Type>
void RNVector<Type>::
Reverse(void)
{
   // reverse order of all entries
   Reverse(0, 0);
}



template <typename Type>
void RNVector<Type>::
Reverse(int start, int length)
{
   // compute number of entries to reserve
   int nreverse = nentries - start;
   if (nreverse <= 0) return;
   if ((length > 0) && (length < nreverse)) nreverse = length;
   if (nreverse <= 0) return;

   // reverse length at start
   int i, j;
   for (i = start, j = start + nreverse - 1; i < j; ++i, --j) {
      Swap(i, j);
   }
}



template <typename Type>
void RNVector<Type>::
Append(const RNVector& vector)
{
   // resize first
   Resize(NEntries() + vector.NEntries());

   // insert entries of vector
   for (int ie = 0; ie < vector.NEntries(); ++ie)
      Insert(vector[ie]);
}



template <typename Type>
void RNVector<Type>::
Shuffle(void)
{
   for (int i = 0; i < nentries; ++i) {
      int r = i + (int)(RNRandomScalar() * (nentries - i));
      // just checking
      rn_assertion((0 <= r) && (r < nentries));
      Swap(i, r);
   }
}


template <typename Type>
void RNVector<Type>::
QuickSort(int(*compare)(const Type data1, const Type data2))
{
   Shuffle();
   Sort(compare, 0, nentries - 1);
   sorted = TRUE;
}



template <typename Type>
RNBoolean RNVector<Type>::
IsSorted(void)
{
   return sorted;
}



template <typename Type>
void RNVector<Type>::
Sort(int(*compare)(const Type data1, const Type data2), int lo, int hi)
{
   if (hi <= lo) return;
   int j = Partition(compare, lo, hi);
   Sort(compare, lo, j - 1);
   Sort(compare, j + 1, hi);
}



template <typename Type>
int RNVector<Type>::
Partition(int(*compare)(const Type data1, const Type data2), int lo, int hi)
{
   int i = lo;
   int j = hi + 1;
   const Type v = *entries[lo];
   while (true) {
      // find item on lo to swap
      while (compare(*entries[++i], v))
         if (i == hi) break;

      // find item on hi to swap
      while (compare(v, *entries[--j]))
         if (j == lo) break;

      // check if pointers cross
      if (i >= j) break;

      Swap(i, j);
   }

   // put partioning item v at entry j
   Swap(lo, j);

   // now a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
   return j;
}



template <typename Type>
void RNVector<Type>::
Swap(int i, int j)
{
   Type *tmp = entries[i];
   entries[i] = entries[j];
   entries[j] = tmp;
}



template <typename Type>
void RNVector<Type>::
Resize(int length)
{
   // check if length is valid
   rn_assertion(nentries <= nallocated);

   // return if length less than nallocated
   if (length < nallocated) return;

   // adjust length to be next greater power of 2
   int tmplength = RN_VECTOR_MIN_ALLOCATED;
   while (tmplength < length)
      tmplength *= 2;
   length = tmplength;

   // allocate new entries
   Type **newentries = NULL;
   newentries = new Type *[length];
   rn_assertion(newentries);

   // copy old entries into new entries
   if (nentries > 0) {
      rn_assertion(entries);
      rn_assertion(newentries);
      for (int ie = 0; ie < nentries; ++ie) {
         newentries[ie] = entries[ie];
      }
   }

   // replace entries
   if (entries) delete[] entries;
   entries = newentries;

   // update nallocated
   nallocated = length;
}



template <typename Type>
void RNVector<Type>::
InternalInsert(Type data, int k)
{
   // make sure there is enough storage for new entry
   Resize(nentries + 1);

   // increment number of entries
   nentries++;

   // shift entries up one notch
   if (k < (nentries - 1)) Shift(k, 0, 1);

   // copy data into kth entry
   entries[k] = new Type(data);
}



template <typename Type>
void RNVector<Type>::
InternalRemove(int k)
{
   // delete the kth entry
   delete entries[k];

   // shift entries down one notch
   if (k < (nentries - 1)) Shift(k + 1, 0, -1);

   // decrement number of entries
   --nentries;
}

#endif
