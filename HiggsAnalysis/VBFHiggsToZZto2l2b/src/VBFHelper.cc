#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHelper.h"


using namespace vbfhzz2l2b;
using namespace reco;
using namespace std;

VBFHelper::VBFHelper( const edm::ParameterSet & iConfig )  :
  parameters_( iConfig ),
  physHistos_(0),
  verboseLevel_(100)
{
}


VBFHelper::~VBFHelper()
{
  if ( physHistos_ ) delete physHistos_;
}

PhysicsHistograms::KinAxisLimits VBFHelper::getAxisLimits(std::string name)
{
  PhysicsHistograms::KinAxisLimits axisLimits;
  if ( parameters_.exists(name) ) {
    edm::ParameterSet axisLimitsSet = parameters_.getParameter<edm::ParameterSet>(name);
    axisLimits.pt1 = axisLimitsSet.getParameter<double>("pt1");
    axisLimits.pt2 = axisLimitsSet.getParameter<double>("pt2");
    axisLimits.m1  = axisLimitsSet.getParameter<double>("m1");
    axisLimits.m2  = axisLimitsSet.getParameter<double>("m2");
  }
  return axisLimits;
}


void VBFHelper::bookHistos(edm::EDProducer * producer)
{
    // Initialize TFileService
  edm::Service<TFileService> fs;
  TFileDirectory summary = TFileDirectory( fs->mkdir("summary") );


  physHistos_ = new PhysicsHistograms ( getAxisLimits("muonAxis"),
					getAxisLimits("electronAxis"),
					getAxisLimits("jetAxis"),
					getAxisLimits("METAxis"),
					getAxisLimits("genParticleAxis")
					);

  // Get list of histograms to enable and disable
  string histos_to_disable =
    parameters_.getParameter<string>    ("disable");
  string histos_to_enable  =
    parameters_.getParameter<string>    ("enable");
  physHistos_->configure( histos_to_disable, histos_to_enable );

  // PhysicsHistograms takes ownership of the memory for these histograms here
  physHistos_->addHisto( h_runNumber_ = 
			new PhysVarHisto( "runNumber", "Run Number",
					  10000, 0, 10000, 
					  &summary ,"", "I") 
			);
  physHistos_->addHisto( h_eventNumber_ = 
			new PhysVarHisto( "eventNumber", "Event Number",
					  10000, 0, 10000, 
					  &summary ,"", "I") 
			);

  // Be sure to make the histograms!
  h_runNumber_->makeTH1();
  h_eventNumber_->makeTH1();

  // &&& Ntuple booking begin

  // Now that we know which variables are in the game, we could also
  // decide which ones to ntuplize
  string list_of_ntuple_vars =
    parameters_.getParameter<std::string>    ("ntuplize");

  if (list_of_ntuple_vars != "") {
    //
    //--- Collect all PhysVarHistos which need to store ntuple
    //--- variables and put them in here.
    physHistos_->select( list_of_ntuple_vars, ntVars_ );

    //--- Iterate over the list and "book" them via EDM
    std::vector< PhysVarHisto* >::iterator
      p    = ntVars_.begin(),
      pEnd = ntVars_.end();

    for ( ; p != pEnd; ++p ) {
      cout << "Adding ntuple variable " << (*p)->name() << endl;
      addNtupleVar( producer, (*p)->name(), (*p)->type() );
    }
    //
  } // end if

}


void VBFHelper::getHandles( edm::Event & event,
			       edm::Handle<reco::MuonCollection >           & muonHandle,	
			       edm::Handle<reco::GsfElectronCollection >    & electronHandle,
			       edm::Handle<reco::CaloJetCollection >        & caloJetHandle, 
			       edm::Handle<reco::CaloMETCollection >        & caloMETHandle,
			       edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle
			       )
{

  bool doMuon         = parameters_.getParameter<bool>("doMuon");
  bool doElectron     = parameters_.getParameter<bool>("doElectron");
  bool doJet          = parameters_.getParameter<bool>("doJet");
  bool doMET          = parameters_.getParameter<bool>("doMET");
  bool doGenParticles = parameters_.getParameter<bool>("doGenParticles");
  
  
  edm::InputTag muonName         = parameters_.getParameter<edm::InputTag>("muonSrc"       );
  edm::InputTag electronName     = parameters_.getParameter<edm::InputTag>("electronSrc"   );
  edm::InputTag jetName          = parameters_.getParameter<edm::InputTag>("jetSrc"        );
  edm::InputTag METName          = parameters_.getParameter<edm::InputTag>("METSrc"        );
  edm::InputTag genParticlesName = parameters_.getParameter<edm::InputTag>("genParticleSrc");

  

  if ( doMuon         ) event.getByLabel(muonName        , muonHandle);
  if ( doElectron     ) event.getByLabel(electronName    , electronHandle);
  if ( doJet          ) event.getByLabel(jetName         , caloJetHandle);
  if ( doMET          ) event.getByLabel(METName         , caloMETHandle);
  if ( doGenParticles ) event.getByLabel(genParticlesName, genParticlesHandle );

}


void VBFHelper::fillHistograms(edm::Event & event,
				  edm::Handle<reco::MuonCollection >           & muonHandle,	 
				  edm::Handle<reco::GsfElectronCollection >    & electronHandle, 
				  edm::Handle<reco::CaloJetCollection >        & caloJetHandle,  
				  edm::Handle<reco::CaloMETCollection >        & caloMETHandle,  
				  edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle    
				  )
{
  physHistos_->clearVec();
  if ( muonHandle.isValid()         ) physHistos_->fillCollection(*muonHandle);
  if ( electronHandle.isValid()     ) physHistos_->fillCollection(*electronHandle);
  if ( caloJetHandle.isValid()      ) physHistos_->fillCollection(*caloJetHandle);
  if ( caloMETHandle.isValid()      ) physHistos_->fillCollection(*caloMETHandle);
  if ( genParticlesHandle.isValid() ) physHistos_->fillCollection(*genParticlesHandle);



  // save the list of ntuple varibles to the event record
  saveNtuple( event, ntVars_ );
}



// &&& Design task: add all data types supported by PhysVarHisto
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
void
VBFHelper::addNtupleVar( edm::EDProducer * producer, std::string name, std::string type )
{
  if      (type == "D") {
    producer->produces<double>( name ).setBranchAlias( name );
  }
  else if (type == "F") {
    producer->produces<float>( name ).setBranchAlias( name );
  }
  else if (type == "I") {
    producer->produces<int>( name ).setBranchAlias( name );
  }
  else if (type == "i") {
    producer->produces<unsigned int>( name ).setBranchAlias( name );
  }
  else if (type == "S") {
    producer->produces<short>( name ).setBranchAlias( name );
  }
  else if (type == "s") {
    producer->produces<unsigned short>( name ).setBranchAlias( name );
  }
  else if (type == "L") {
    producer->produces<long>( name ).setBranchAlias( name );
  }
  else if (type == "l") {
    producer->produces<unsigned long>( name ).setBranchAlias( name );
  }
  else if (type == "vD") {
    producer->produces<vector<double> >( name ).setBranchAlias( name );
  }
  else if (type == "vF") {
    producer->produces<vector<float> >( name ).setBranchAlias( name );
  }
  else if (type == "vI") {
    producer->produces<vector<int> >( name ).setBranchAlias( name );
  }
  else if (type == "vi") {
    producer->produces<vector<unsigned int> >( name ).setBranchAlias( name );
  }
  else if (type == "vS") {
    producer->produces<vector<short> >( name ).setBranchAlias( name );
  }
  else if (type == "vs") {
    producer->produces<vector<unsigned short> >( name ).setBranchAlias( name );
  }
  else if (type == "vL") {
    producer->produces<vector<long> >( name ).setBranchAlias( name );
  }
  else if (type == "vl") {
    producer->produces<vector<unsigned long> >( name ).setBranchAlias( name );
  }
  else {
    std::cout << "VBFHelper::addNtupleVar (ERROR): "
	      << "unknown type " << type << std::endl;

    // &&& Throw an exception in order to abort the job!
  }
}







void
VBFHelper::saveNtuple( edm::Event & event,
			  const std::vector<vbfhzz2l2b::PhysVarHisto*> & ntvars )
{
  //  Ntuplization

  if ( verboseLevel_ > 0 )
    cout << "About to save ntuple" << endl;
  if ( ntvars.size() ) {

    //--- Iterate over the list and "fill" them via EDM
    std::vector< PhysVarHisto* >::const_iterator
      p    = ntvars.begin(),
      pEnd = ntvars.end();

    for ( ; p != pEnd; ++p ) {

      if      ((*p)->type() == "D") {
	saveNtupleVar<double>( event, (*p)->name(), (*p)->value() );
      }
      else if ((*p)->type() == "F") {
	saveNtupleVar<float>( event, (*p)->name(), (*p)->value() );
      }
      else if ((*p)->type() == "I") {
	saveNtupleVar<int>( event, (*p)->name(), static_cast<int>((*p)->value()) );
      }
      else if ((*p)->type() == "i") {
	saveNtupleVar<unsigned int>( event, (*p)->name(), static_cast<unsigned int>((*p)->value()) );
      }
      else if ((*p)->type() == "S") {
	saveNtupleVar<short>( event, (*p)->name(), static_cast<short>((*p)->value()) );
      }
      else if ((*p)->type() == "s") {
	saveNtupleVar<unsigned short>( event, (*p)->name(), static_cast<unsigned short>((*p)->value()) );
      }
      else if ((*p)->type() == "L") {
	saveNtupleVar<long>( event, (*p)->name(), static_cast<long>((*p)->value()) );
      }
      else if ((*p)->type() == "l") {
	saveNtupleVar<unsigned long>( event, (*p)->name(), static_cast<unsigned long>((*p)->value()) );
      }
      else if ((*p)->type() == "vD") {
	vector<double> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<double>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vF") {
	vector<float> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<float>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vI") {
	vector<int> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<int>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vi") {
	vector<unsigned int> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<unsigned int>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vS") {
	vector<short> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<short>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vs") {
	vector<unsigned short> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<unsigned short>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vL") {
	vector<long> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<long>( event, (*p)->name(), retvec );
      }
      else if ((*p)->type() == "vl") {
	vector<unsigned long> retvec;
	(*p)->vec( retvec );
	saveNtupleVec<unsigned long>( event, (*p)->name(), retvec );
      }
      else {
	std::cout << "VBFHelper::addNtupleVar (ERROR): "
		  << "unknown type " << std::endl;

	// &&& Throw an exception in order to abort the job!
      }


      (*p)->clearVec(); // Clear ntuple cache here as well

    }
    //
  } // end if
}


template <class T>
void
VBFHelper::saveNtupleVar( edm::Event & event,
			     std::string name, T value )
{
  std::auto_ptr<T> aptr( new T (value ) );
  if ( verboseLevel_ > 0 )
    cout << "Putting variable " << name << " with value " << value << " in event" << endl;
  event.put( aptr, name );
}


template <class T>
void
VBFHelper::saveNtupleVec( edm::Event & event,
			       std::string name, const vector<T> & value )
{
  std::auto_ptr<vector<T> > aptr( new vector<T> ( value ) );
  if ( verboseLevel_ > 0 ) {
    cout << "Putting variable " << name << " with values : " << endl;
    typename vector<T>::const_iterator i = aptr->begin(), iend = aptr->end();
    for ( ; i != iend; ++i ) {
      cout << *i << " ";
    }
    cout << endl;
  }
  event.put( aptr, name );
}

