#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCprocessFilter.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

// utilities
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"

#include <iostream>

using namespace vbfhzz2l2b;


//! constructor
VBFHZZllbbMCprocessFilter::VBFHZZllbbMCprocessFilter(const edm::ParameterSet& iConfig) :
  whichSim_ ( iConfig.getParameter<int>( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  signal_   ( iConfig.getParameter<int>( "signal"    ) )  // 1:Signal,  0:Background
{
  vbfFlag_ = true;
  if (signal_)
    vbfFlag_ = iConfig.getParameter<bool>( "vbfFlag" );
  //  std::cout << "vbfFlag_: " << vbfFlag_ << std::endl;
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

  if (!signal_) return true;

  int processID = 0;
  
  if ( whichSim_ == vbfhzz2l2b::FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    processID = mcEv->signal_process_id();
  }
  else if ( whichSim_ == vbfhzz2l2b::FASTSIM ) {
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

  //  std::cout << "processID: " << processID << std::endl;
  
  if(vbfFlag_) {
    if (processID == HZZFusion_ || processID == HWWFusion_) 
      return true ; 
  } else
    if (processID == HggFusion_)
      return true ;
  
  return false ;
}
	

