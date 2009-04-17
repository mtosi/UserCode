
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PhysVarProfile.h"

#include <string>
#include <sstream>
#include <iostream>

using vbfhzz2l2b::PhysVarProfile;

// Constructor:
PhysVarProfile::PhysVarProfile( std::string name,
				std::string title,
				int         xnbins,
				double      xlow,
				double      xhigh,
				double      ylow,
				double      yhigh,
				TFileDirectory * currDir,
				std::string units,
				std::string type,
				bool        saveHist,
				bool        saveNtup,
				bool        isMC ) :
  currDir_ (currDir),
  name_   (name),
  type_   (type),
  title_  (title),
  xnbins_ (xnbins),
  xlow_   (xlow),
  xhigh_  (xhigh),
  ylow_   (ylow),
  yhigh_  (yhigh),
  units_  (units),
  isMC_   (isMC),
  //
  histos_    (),             // vector of 0 elements
  xvalue_     (-999999),
  yvalue_     (-999999),
  xvalue_ext_ (0),
  yvalue_ext_ (0),
  indices_   (0),
  //
  saveHist_ (saveHist),
  saveNtup_ (saveNtup),
  //
  verboseLevel_ (0)      // default = very verbose (&&& hardcoded)
{
  if (verboseLevel_ > 5)
    std::cout << "PhysVarProfile(" << name_ << "):: in constructor." << std::endl;

  // &&& Not sure what (if anything) will be needed to deal with the ntuple.
}


void PhysVarProfile::makeProfile(int nameId) {
  if (verboseLevel_ > 5)
    std::cout << "PhysVarProfile(" << name_ << "):: in makeProfile()." << std::endl;


  if (saveHist_ && currDir_) {
    if (verboseLevel_ > 5)
      std::cout << "PhysVarProfile(" << name_ << ")::makeProfile:"
		<< " making new histo." << std::endl;

    std::string name, title;
    std::stringstream partNum;
    if ( nameId < 0 )
      partNum << histos_.size()+1;
    else
      partNum << nameId;
    
    name  = name_  + "_" + partNum.str();
    title = title_ + " #" + partNum.str()+ ";" + units_ ;


    TProfile * h = currDir_->make<TProfile>( name.c_str(), title.c_str(),
					     xnbins_, xlow_, xhigh_,
					     ylow_, yhigh_ );
    // &&& Now decorate: set thickness, bkg color, axis labels, etc.
    // &&& lift the axis labels from RooFit which does a nice job on them

    histos_.push_back( h );
    if ( isMC_ ) 
      indices_.push_back( nameId );
  }
}


// Note: imulti follows Fortran-style indexing,
// so imulti==1 is the 1st element.  Note that the histos_ vector
// wakes up with no elements, so histos_.size() gives 0 in the first
// call to fill().
void
PhysVarProfile::fill( double x, double y, unsigned int imulti, double weight ) {
  // Once this PhysVarProfile was made, then there's no turning back
  // when it comes to fill() : we cache the value no matter what.

  xvalue_ = x;  // save in case we need to capture it for the ntuple.
  yvalue_ = y;  // save in case we need to capture it for the ntuple.

  xvalueColl_.push_back( x );
  yvalueColl_.push_back( y );

  if ( xvalue_ext_ ) {
    // Defined, use this one as well
    *( (double*)(xvalue_ext_) ) = xvalue_ ;
  }
  if ( yvalue_ext_ ) {
    // Defined, use this one as well
    *( (double*)(yvalue_ext_) ) = yvalue_ ;
  }

  if (verboseLevel_ > 5) {  // very verbose
    std::cout << "PhysVarProfile(" << name_ << ")::fill:"
	      << " xvalue = " << xvalue_ 
	      << " yvalue = " << yvalue_ 
	      << ". External storage are ";
    if (xvalue_ext_ && yvalue_ext_)  std::cout << "defined." << std::endl;
    else                             std::cout << "not defined." << std::endl;
  }



  if ( ! saveHist_ ) {
    std::cout << "PhysVarProfile(" << name_ << ")::fill: saveHist is not defined" << std::endl;
    return;
  }

  //--- If the requested histo location in multi-histogram
  //--- mode is larger than our current vector of histograms,
  //--- then grow it up to the location requested.
  //

  if ( !isMC_ ) {
    while ( imulti > histos_.size() ) {
      if (verboseLevel_ > 4) {
	std::cout << "PhysVarProfile(" << name_ << ")::fill: grow histo list by one."
		  << std::endl;
      }
      //--- Make another TProfile at index == current value of size()
      makeProfile();    // histos_ vector grows by one element
    }
    if ( verboseLevel_ > 4 )
      std::cout << "PhysVarProfile(" << name_ << ")::fill: About to fill " << imulti << std::endl;
    histos_[ imulti-1 ]->Fill( xvalue_, yvalue_, weight );
    if ( verboseLevel_ > 4 )
      std::cout << "PhysVarProfile(" << name_ << ")::fill: Done filling " << imulti << std::endl;
  } 
  else {
    
    std::vector<int>::const_iterator ifound = find( indices_.begin(), indices_.end(), imulti );
    if ( ifound == indices_.end() ) {
      makeProfile(imulti);
    } else {
      unsigned int index = ifound - indices_.begin();
      histos_[index]->Fill( xvalue_, yvalue_, weight );
    }
  }
}

