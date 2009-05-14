#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/loadSimpleNtple.h"


void loadSimpleNtuple(int index,TChain* chain,
		      TClonesArray* evtClass,
		      TClonesArray* jetClass,
		      TClonesArray* muonClass,
		      TClonesArray* electronClass,
		      TClonesArray* zhadClass) {

  using namespace vbfhzz2l2b;
  using namespace vbfhzz2l2b::SimpleNtpleObj;
  
  chain->GetEntry(index);
  EVT* event = (EVT*) evtClass->At(0);

  
  EVT_Run            = event->Run;            // run number
  EVT_Event          = event->Event;          // event number
  EVT_Ilum           = event->Ilum;           // instantaneous luminosity (e30)
  EVT_eventID        = event->eventID;        // event ID
                                              // VBF:123 or 124
                                              // ggF: 102
  EVT_nPV            = event->nPV;            // number of primary vertex
  EVT_trigpath       = event->trigpath;       // Z_BB Trigger Path: 1*main + 10*test1 + 100*test2 (if exists). Ex 101 means main + 2nd test trigger where fired.
  
  EVT_P4bquark1      = event->P4bquark1;      // first  b-quark
  EVT_P4bquark2      = event->P4bquark2;      // second b-quark
  EVT_indjetb1       = event->indjetb1;       // first b-quark associated jet index
  EVT_indjetb2       = event->indjetb2;       // second b-quark associated jet index
  
  EVT_pthat          = event->pthat;          // pthat
  
  EVT_njet           = event->njet;           // number of jets in the event
  EVT_ngoodjet       = event->ngoodjet;       // number of "good" jets in the event
  EVT_nbtag          = event->nbtag;          // number of tag in the event
  EVT_Zvertex        = event->Zvertex;        // Z of the reconstructed primary vertex
  EVT_P2met          = event->P2met;          // Missing Et vector

  EVT_nmuon          = event->nmuon;          // number of muons in the event	 
  EVT_nelectron      = event->nelectron;      // number of electrons in the event
  EVT_nZhad          = event->nZhad;          // number of hadronic Z            
  
  
  // jet block
  // ---------
  for ( int jetIndex = 0; jetIndex < EVT_njet; jetIndex++ ) {
    JET* jet = (JET*) jetClass->At(jetIndex);
    JET_P4jetl0.push_back     ( jet->P4jetl0     );         // 4-momentum of jet for Jet correction L0
    JET_P4jetl23.push_back    ( jet->P4jetl23    );         // 4-momentum of jet for Jet correction L2+L3
    JET_etjetl0.push_back     ( jet->etjetl0     );         // Et of jet for Jet Correction L0
    JET_etjetl23.push_back    ( jet->etjetl23    );         // Et of jet for Jet Correction L2+L3
    JET_etajet.push_back      ( jet->etajet      );         // Eta of jet
    JET_detetajet.push_back   ( jet->detetajet   );         // Detector eta of jet
    JET_rapidityjet.push_back ( jet->rapidityjet );         // Rapidity of jet
    JET_emfjet.push_back      ( jet->emfjet      );         // EMfraction of jet
    JET_chfjet.push_back      ( jet->chfjet      );         // Charge fraction of jet (p/E)
    JET_massjet.push_back     ( jet->massjet     );         // Mass of jet
    JET_phijet.push_back      ( jet->phijet      );         // Phi of jet
    JET_etaetajet.push_back   ( jet->etaetajet   );         // Eta-Eta moment
    JET_etaphijet.push_back   ( jet->etaphijet   );         // Eta-Phi moment
    JET_phiphijet.push_back   ( jet->phiphijet   );         // Phi-Phi moment
    JET_metprjjet.push_back   ( jet->metprjjet   );         // MET projected in the direction to jet
						    
    JET_ntrk_tot.push_back    ( jet->ntrk_tot    );      // The number of total track
    JET_P4trk_tot.push_back   ( jet->P4trk_tot   );      // Sum of 4-momentum of total track 
    JET_pt_tot.push_back      ( jet->pt_tot      );      // Pt of total track
    JET_trkmass_tot.push_back ( jet->trkmass_tot );      // Mass calculated from total track             
    JET_sumPt_tot.push_back   ( jet->sumPt_tot   );       
    JET_ch_tot.push_back      ( jet->ch_tot      );      // Sum charge of total track
						    
    JET_bDiscr.push_back      ( jet->bDiscr      );
    JET_ntrk_tag.push_back    ( jet->ntrk_tag    );      // The number of SVX track
    JET_P4trk_tag.push_back   ( jet->P4trk_tag   );      // Sum of 4-momentum of SVX track 
    JET_pt_tag.push_back      ( jet->pt_tag      );      // Pt of SVX track
    JET_trkmass_tag.push_back ( jet->trkmass_tag );      // Mass calculated from SVX track
    JET_ch_tag.push_back      ( jet->ch_tag      );      // Sum charge of SVX track
    JET_lxy_tag.push_back     ( jet->lxy_tag     );      // Lxy 
    JET_L3d_tag.push_back     ( jet->L3d_tag     );      // 3 dimensional B hadron flight direction
    JET_tag.push_back         ( jet->tag         );      // tag -1:negatively tag
                                                         //      0:taggable
                                                         //      1:positively tag

    JET_SecVtxVertex.push_back ( jet->SecVtxVertex ); // secondary vertex 

    JET_nele.push_back     ( jet->nele     );      // The number of electron into jet
    JET_P4ele1.push_back   ( jet->P4ele1   );      // 4-momentum of electron 1st into jet
    JET_P4ele2.push_back   ( jet->P4ele2   );      // 4-momentum of electron 2nd into jet
    JET_pt_ele1.push_back  ( jet->pt_ele1  );      // Pt of first electron into jet
    JET_pt_ele2.push_back  ( jet->pt_ele2  );      // Pt of second electron into jet
    JET_ptr_ele1.push_back ( jet->ptr_ele1 );      // Ptrel of first electron into jet
    JET_ptr_ele2.push_back ( jet->ptr_ele2 );      // Ptrel of second electron into jet
    JET_nmuon.push_back    ( jet->nmuon    );      // The number of muon into jet
    JET_P4muo1.push_back   ( jet->P4muo1   );      // 4-momentum of muon 1st into jet 
    JET_P4muo2.push_back   ( jet->P4muo2   );      // 4-momentum of muon 2nd into jet
    JET_pt_mu1.push_back   ( jet->pt_mu1   );      // Pt of first muon into jet
    JET_pt_mu2.push_back   ( jet->pt_mu2   );      // Pt of second muon into jet
    JET_ptr_mu1.push_back  ( jet->ptr_mu1  );      // Ptrel of first muon into jet
    JET_ptr_mu2.push_back  ( jet->ptr_mu2  );      // Ptrel of second muon into jet
		           			     
    JET_nbqrk.push_back   ( jet->nbqrk   );           // 
    JET_dr_bqrk.push_back ( jet->dr_bqrk );          // Delta R between jet and nearest b quark
    JET_P4bqrk.push_back  ( jet->P4bqrk  );          // 4-momentum of nearest b quark (pttrue)
    JET_pt_bqrk.push_back ( jet->pt_bqrk );          //
    JET_assjet.push_back  ( jet->assjet  );          // The relation between quark and jet 0:not match, 
                                                     //                                    1:first b-quark matched,
                                                     //                                    2:second b-quark matched    
    JET_pt_bhad.push_back ( jet->pt_bhad );          // Estimated b-hadron pt
    JET_pt_neu.push_back  ( jet->pt_neu  );          // neutral Pt of B hadron

  }

  // hadronic Z block
  // ----------------
  for ( int zhadIndex = 0; zhadIndex < EVT_nZhad; zhadIndex++ ) {
    ZHAD* zhad = (ZHAD*) zhadClass->At(zhadIndex);

    ZHAD_P4Zhadl0.push_back             ( zhad->P4Zhadl0             ); // 4-momentum of hadronic Z for Jet correction L0
    ZHAD_P4Zhadl23.push_back            ( zhad->P4Zhadl23            );  // 4-momentum of hadronic Z for jet correction L2+L3
    ZHAD_etZhadl0.push_back             ( zhad->etZhadl0             ); // Et of hadronic Z for jet Correction L0
    ZHAD_etZhadl23.push_back            ( zhad->etZhadl23            ); // Et of hadronic Z for jet Correction L2+L3
    ZHAD_etaZhad.push_back              ( zhad->etaZhad              ); // Eta of hadronic Z
    ZHAD_detetaZhad.push_back           ( zhad->detetaZhad           ); // Detector eta of hadronic Z
    ZHAD_rapidityZhad.push_back         ( zhad->rapidityZhad         ); // Rapidity of hadronic Z
    ZHAD_emfZhad.push_back              ( zhad->emfZhad              ); // EMfraction of hadronic Z
    ZHAD_chfZhad.push_back              ( zhad->chfZhad              ); // Charge fraction of hadronic Z (p/E)
    ZHAD_phiZhad.push_back              ( zhad->phiZhad              ); // Phi of hadronic Z
    ZHAD_etaetaZhad.push_back           ( zhad->etaetaZhad           ); // Eta-Eta moment
    ZHAD_etaphiZhad.push_back           ( zhad->etaphiZhad           ); // Eta-Phi moment
    ZHAD_phiphiZhad.push_back           ( zhad->phiphiZhad           ); // Phi-Phi moment
    ZHAD_metprjZhad.push_back           ( zhad->metprjZhad           ); // MET projected in the direction to hadronic Z
    ZHAD_massZhad_l0.push_back          ( zhad->massZhad_l0          ); // invariant mass corrected L0
    ZHAD_massZhad_l23.push_back         ( zhad->massZhad_l23         ); // invariant mass corrected L2+L3
    ZHAD_collinearityZhad_l0.push_back  ( zhad->collinearityZhad_l0  ); //
    ZHAD_collinearityZhad_l23.push_back ( zhad->collinearityZhad_l23 ); //

  }

  // muon block
  // ---------
  for ( int muonIndex = 0; muonIndex < EVT_nmuon; muonIndex++ ) {
    MUON* muon = (MUON*) muonClass->At(muonIndex);

    MUON_P4muon.push_back            ( muon->P4muon            );  // 4-momentum of muon
    MUON_isolMuonSumPt.push_back     ( muon->isolMuonSumPt     );  // sum of tracks Pt associated to isolated muon
    MUON_isolMuonTrkNumber.push_back ( muon->isolMuonTrkNumber );  // number of tracks associated to isolated muon

  }

  // electron block
  // ---------
  for ( int electronIndex = 0; electronIndex < EVT_nelectron; electronIndex++ ) {
    ELECTRON* electron = (ELECTRON*) electronClass->At(electronIndex);

    ELECTRON_P4electron.push_back       ( electron->P4electron       );   // 4-momentum of electron
    ELECTRON_isolEleSumPt.push_back     ( electron->isolEleSumPt     );   // sum of tracks Pt associated to isolated electron
    ELECTRON_isolEleTrkNumber.push_back ( electron->isolEleTrkNumber );   // number of tracks associated to isolated electron
    ELECTRON_eleId.push_back            ( electron->eleId            );   // electron ID

  }



}

