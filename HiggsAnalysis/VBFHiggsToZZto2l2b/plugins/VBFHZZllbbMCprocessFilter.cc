#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCprocessFilter.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include <iostream>

enum { FASTSIM = 0,
       FULLSIM = 1
};

//! constructor
VBFHZZllbbMCprocessFilter::VBFHZZllbbMCprocessFilter(const edm::ParameterSet& iConfig) :
  whichSim_ ( iConfig.getParameter<int>( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  signal_   ( iConfig.getParameter<int>( "signal"    ) )  // 0:Signal,  1:Background
{
}


// ------------------------------------------------------------------------------------


//! destructor
VBFHZZllbbMCprocessFilter::~VBFHZZllbbMCprocessFilter()
{}


// ------------------------------------------------------------------------------------


//! filtering method
bool 
VBFHZZllbbMCprocessFilter::filter (edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace vbfhzz2l2b;

  int processID = 0;
  
  if ( whichSim_ == FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    processID = mcEv->signal_process_id();
  }
  else if ( whichSim_ == FASTSIM ) {
    edm::Handle<int> genProcessID;
    try {
      iEvent.getByLabel( "genEventProcID", genProcessID ); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FastSim, but genEventProcID not found" << std::endl; }
    
    processID = *genProcessID;
  }
  else {
    std::cout << "--> WARNING: simulation not specificied!!" << std::endl;
  }

  std::cout << "processID: " << processID << std::endl;
  
  if (processID == HZZFusion_ || processID == HWWFusion_) return true ;
  // if (processID == HggFusion_) return true ;
  return false ;

}
	

