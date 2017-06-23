// Include file for vector

#ifndef __RN__VECTOR__H__
#define __RN__VECTOR__H__



/////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
/////////////////////////////////////////////////////////////////////

template <typename Type>
class RNVector {
public:
   //////////////////////////////////
   //// CONSTRUCTORS/DESTRUCTORS ////
   //////////////////////////////////

   RNVector(void);
   ~RNVector(void);

   // vector property functions/operators
   const RNBoolean IsEmpty(void) const;
   const int NAllocated(void) const;
   const int NEntries(void) const;
   const int Size(void) const;

   // data access functions/operators
   const Type& Head(void) const;
   const Type& Kth(int k) const;
   const Type& operator[](int k) const;
   const Type& Tail(void) const;

   // insertion function/operators
   void Insert(Type data);
   void InsertHead(Type data);
   void InsertKth(Type data, int k);
   void InsertTail(Type data);

   // removal functions/operators
   void Remove(void);
   void RemoveHead(void);
   void RemoveKth(int k);
   void RemoveTail(void);

   // manipulation functions
   void Empty(RNBoolean deallocate = FALSE);
   void Truncate(int length);
   void Shift(int delta);
   void Shift(int start, int length, int delta);

   void Reverse(void);
   void Reverse(int start, int length);
   void Append(const RNVector& vector);
   void Shuffle(void);
   void QuickSort(int(*compare)(const Type data1, const Type data2));
   RNBoolean IsSorted(void);

   void Swap(int i, int j);
   void Resize(int length);

protected:
   // internal functions
   void InternalInsert(Type data, int k);
   void InternalRemove(int k);

   // sorting functions
   void Sort(int(*compare)(const Type data1, const Type data2), int lo, int hi);
   int Partition(int(*compare)(const Type data1, const Type data2), int lo, int hi);

private:
   Type **entries;
   int nallocated;
   int nentries;
   RNBoolean sorted;
};



// include files
#include "RNVector.cpp"



#endif
