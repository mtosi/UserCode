#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_PhysVarProfile_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_PhysVarProfile_h 1

// &&& Design comments:
//     Here's a list of types accepted by ROOT:
//             - C : a character string terminated by the 0 character
//             - B : an 8 bit signed integer (Char_t)
//             - b : an 8 bit unsigned integer (UChar_t)
//             - S : a 16 bit signed integer (Short_t)
//             - s : a 16 bit unsigned integer (UShort_t)
//             - I : a 32 bit signed integer (Int_t)
//             - i : a 32 bit unsigned integer (UInt_t)
//             - F : a 32 bit floating point (Float_t)
//             - D : a 64 bit floating point (Double_t)
//             - L : a 64 bit signed integer (Long64_t)
//             - l : a 64 bit unsigned integer (ULong64_t)


// STL include files
#include <string>
#include <vector>

// ROOT include files
#include <TProfile.h>

// &&& Alert!  Is this a dependence on full framework?
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

namespace vbfhzz2l2b {

  class PhysVarProfile
  {
  public:
    PhysVarProfile( std::string name,
	       std::string title,
	       int         xnbins,
	       double      xlow,
	       double      xhigh,
	       double      ylow,
	       double      yhigh,
	       TFileDirectory * currDir = 0,
	       std::string units = "",
	       std::string type  = "D",
	       bool        saveHist = true,
	       bool        saveNtup = false,
	       bool        isMC     = false );

    virtual ~PhysVarProfile() { };  //!  Note we don't delete histograms!

    //--- Make one TH1 (may need more than one); all decorations should be done in this call.
    virtual void makeProfile(int nameId = -1);

    //--- Fill one of the histograms in histos_ vector.
    virtual void fill( double x,
		       double y,
		       unsigned int imulti = 1,
		       double weight = 1.0 );

    //--- Inline accessors.
    inline std::string name() { return name_ ; } //!< ROOT/cfg handle
    inline std::string type() { return type_ ; } //!< type of value_
    inline double xvalue() { return xvalue_ ; }    //!< current value
    inline double yvalue() { return yvalue_ ; }    //!< current value

    template <class T>
      void xvec(std::vector<T> & retVec)
    {
      retVec.resize( xvalueColl_.size() );
      for ( unsigned int i = 0; i < xvalueColl_.size(); i++ )
	retVec[i] = static_cast<T>( xvalueColl_[i] );
    } //!< vector of current values in a list
    template <class T>
      void yvec(std::vector<T> & retVec)
    {
      retVec.resize( yvalueColl_.size() );
      for ( unsigned int i = 0; i < yvalueColl_.size(); i++ )
	retVec[i] = static_cast<T>( yvalueColl_[i] );
    } //!< vector of current values in a list

    inline bool   saveHist() { return saveHist_ ; }
    inline bool   saveNtup() { return saveNtup_ ; }

    inline void   setSaveHist(bool flag) { saveHist_ = flag; } //!< save into a histogram
    inline void   setSaveNtup(bool flag) { saveNtup_ = flag; } //!< save into a ntuple

    inline void   setTFileDirectory(TFileDirectory * dir) { currDir_ = dir; };

    inline void   clearVec() { xvalueColl_.clear();  yvalueColl_.clear(); }

  private:
    //--- Stuff needed to book one histogram
    TFileDirectory * currDir_ ;  //!< ROOT thingy that makes/manages histograms
    std::string name_ ;   //!< ROOT handle, but used in cfg files as well...
    std::string type_ ;   //!< stuff for TBranch() constructor (e.g. "F4")
    std::string title_ ;  //!< nice, descriptive title
    int         xnbins_ ;  //!< num of bins
    double      xlow_ ;    //!< min value
    double      xhigh_ ;   //!< max value
    double      ylow_ ;    //!< min value
    double      yhigh_ ;   //!< max value
    std::string units_ ;  //!< "GeV/c^{2}" etc. for axis labels
    bool        isMC_;    //!< Is this for Gen Particles?

    //--- Cache to make histograms
    std::vector<TProfile *> histos_ ;   // maybe use a base class TH1* ?
    // &&& Should we template all this?

    //--- Internal cache
    double      xvalue_ ;        // our own cache
    double      yvalue_ ;        // our own cache
    void *      xvalue_ext_ ;    // cache is in a struct elsewhere, for TBranch
    void *      yvalue_ext_ ;    // cache is in a struct elsewhere, for TBranch
    std::vector<double> xvalueColl_; // our own cache of a list of values
    std::vector<double> yvalueColl_; // our own cache of a list of values
    std::vector<int>    indices_; // Used for GenParticles, this is the index map

    //--- Flags to control behavior
    // bool        active_ ;    // no clear use case to have this flag...
    bool        saveHist_ ;     // save info into a histogram
    bool        saveNtup_ ;     // save info into a ntuple

    int         verboseLevel_ ; // how much verbosity
  };

}

#endif
