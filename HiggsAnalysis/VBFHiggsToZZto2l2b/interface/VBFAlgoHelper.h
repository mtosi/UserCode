#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_VBFAlgoHelper_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_VBFAlgoHelper_h

// Various simple tools

#include <limits>
#include <iostream>
#include <cmath>

#include "PhysicsTools/Utilities/interface/PtComparator.h"
#include "PhysicsTools/Utilities/interface/EtComparator.h"

namespace {

  struct SortObject {
    double value;
    unsigned index;
  };
  
  inline bool sort_lt (const SortObject& a, const SortObject& b) {return a.value < b.value;}
  inline bool sort_gt (const SortObject& a, const SortObject& b) {return a.value > b.value;}
  
  template <class T> struct GetPt {inline double getValue (const T& a) {return a.pt();}};
  template <class T> struct GetEt {inline double getValue (const T& a) {return a.et();}};
  template <class T> struct GetPtRef {inline double getValue (const T& a) {return a->pt();}};
  template <class T> struct GetEtRef {inline double getValue (const T& a) {return a->et();}};
  
  template <class T, class GetValue>
    inline void sortGreater (std::vector <T>* container) {
    std::vector <SortObject> sortable (container->size());
    bool sorted = true;
    GetValue getter;
    for (unsigned i = 0; i < container->size(); i++) {
      sortable[i].value = getter.getValue ((*container)[i]);
      sortable[i].index = i;
      if (sorted && i && sort_gt (sortable[i-1], sortable[i])) sorted = false;
    }
    if (!sorted) { // needs sorting
      std::sort (sortable.begin(), sortable.end(), sort_gt);
      std::vector <T> result;
      result.reserve(container->size());
      for (unsigned i = 0; i < container->size(); i++) {
        result.push_back ((*container)[sortable[i].index]);
      }
      container->swap (result);
    }
  }
  
  template <class T>
    inline void sortByPt (std::vector <T>* container) {
    sortGreater <T, GetPt<T> > (container);
  }
  
  template <class T>
    inline void sortByEt (std::vector <T>* container) {
    sortGreater <T, GetEt<T> > (container);
  }
  
  
  template <class T>
    inline void sortByPtRef (std::vector <T>* container) {
    sortGreater <T, GetPtRef<T> > (container);
  }
  
  template <class T>
    inline void sortByEtRef (std::vector <T>* container) {
    sortGreater <T, GetEtRef<T> > (container);
  }
  
}

template <class T>
class GreaterByPtRef {
 public:
  int operator()(const T& a1, const T& a2) {
    if (!a1) return 0;
    if (!a2) return 1;
    NumericSafeGreaterByPt <typename T::value_type> comp;
    return comp.operator () (*a1, *a2);
  }
};
template <class T>
class GreaterByPtPtr {
 public:
  int operator()(const T* a1, const T* a2) {
    if (!a1) return 0;
    if (!a2) return 1;
    NumericSafeGreaterByPt <T> comp;
    return comp.operator () (*a1, *a2);
  }
};

template <class T>
class GreaterByEtRef {
 public:
  int operator()(const T& a1, const T& a2) {
    if (!a1) return 0;
    if (!a2) return 1;
    NumericSafeGreaterByEt <typename T::value_type> comp;
    return comp.operator () (*a1, *a2);
  }
};

#endif
