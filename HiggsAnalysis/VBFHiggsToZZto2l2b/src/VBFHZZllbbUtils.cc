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
  
  //  double deltaPhi (double phi1, double phi2) {    
  //    double deltaphi=fabs(phi1-phi2);
  //    if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  //    if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  //    return deltaphi;
  //  }
  
  bool Bhadrontable(int pdgcode) {
    bool isBhadron = false;
    int Bmeson[53] = {511,521,10511,10521,513,523,10513,10523,20513,20523,515,
		      525,531,10531,533,10533,20533,535,541,10541,543,10543,
		      20543,545,551,10551,100551,110551,200551,210551,553,10553,
		      20553,30553,100553,110553,120553,130553,200553,210553,
		      220553,300553,9000553,9010553,555,10555,20555,100555,
		      110555,120555,200555,557,10557};
    
    int Bbaryon[35] = {5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,
		       5314,5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,
		       5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,5542,
		       5544,5554};
    
    for(int i=0;i<53;i++){
      if(abs(pdgcode) == Bmeson[i]){
	isBhadron = true;
	return isBhadron;
      }
    }
    for(int i=0;i<35;i++){
      if(abs(pdgcode) == Bbaryon[i]){
	isBhadron = true;
	return isBhadron;
      }
    }
    return isBhadron;
  }
  
}

