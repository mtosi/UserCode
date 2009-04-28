#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include <iostream>

#include "TLorentzVector.h"

namespace vbfhzz2l2b
{

  int bTaggerCode ( const std::string& bTagger ) {
    
    int code = -1;
    if ( bTagger == "HIGHEFF" )            code = HIGHEFF;
    else if ( bTagger == "HIGHPUR"       ) code = HIGHPUR;
    else if ( bTagger == "COMBSECVTX"    ) code = COMBSECVTX;
    else if ( bTagger == "COMBSECVTXMVA" ) code = COMBSECVTXMVA;
    else if ( bTagger == "SOFTMUON"      ) code = SOFTMUON;
    else if ( bTagger == "SOFTELECTRON"  ) code = SOFTELECTRON;
    else if ( bTagger == "JETPROB"       ) code = JETPROB;
    else
      std::cout << "[VBFHZZllbbUtils::bTaggerCode] --> WARNING: bTagger " << bTagger << " NOT IMPLEMENTED!" << std::endl;
    
    return code;
 }

  double resolution ( double & recValue, double & refValue ) {
    return (recValue-refValue)/refValue;
  }
  

void setMomentum (TLorentzVector & myvector, 
		  const reco::Candidate & gen) {
  myvector.SetPx (gen.px ()) ;
  myvector.SetPy (gen.py ()) ;
  myvector.SetPz (gen.pz ()) ;
  myvector.SetE (gen.energy ()) ;
}

double 
deltaPhi (double phi1, double phi2)
{

  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}


}

