#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <vector>

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleNtpleObj.h"

using namespace std;


  int    EVT_Run;                   // run number
  int    EVT_Event;                 // event number
  double EVT_Ilum;                  // instantaneous luminosity (e30)
  int    EVT_eventID;               // event ID
                                    // VBF:123 or 124
                                    // ggF: 102
  int    EVT_nPV;                   // number of primary vertices
  int    EVT_trigpath;              // Z_BB Trigger Path: 1*main + 10*test1 + 100*test2 (if exists). 
                                    // ex. 101 means main + 2nd test trigger where fired.
  double  EVT_pthat;                // pthat

  TLorentzVector EVT_P4bquark1;     // first  b-quark
  TLorentzVector EVT_P4bquark2;     // second b-quark
  int            EVT_indjetb1;      // first b-quark associated jet index
  int            EVT_indjetb2;      // second b-quark associated jet index

  int      EVT_njet;                // number of jets in the event
  int      EVT_ngoodjet;            // number of "good" jets in the event
  int      EVT_nbtag;               // number of tag in the event
  double   EVT_Zvertex;             // Z of the reconstructed primary vertex
  TVector2 EVT_P2met;               // Missing Et vector

  int EVT_nmuon;     // number of muons in the event	 
  int EVT_nelectron; // number of electrons in the event
  int EVT_nZhad;     // number of hadronic Z            


  // JET
  //-----
  std::vector<TLorentzVector> JET_P4jetl0;      // 4-momentum of jet for Jet correction L0
  std::vector<TLorentzVector> JET_P4jetl23;     // 4-momentum of jet for Jet correction L2+L3
  std::vector<double>         JET_etjetl0;      // Et of jet for Jet Correction L0
  std::vector<double>         JET_etjetl23;     // Et of jet for Jet Correction L2+L3
  std::vector<double>         JET_etajet;       // Eta of jet
  std::vector<double>         JET_detetajet;    // Detector eta of jet
  std::vector<double>         JET_rapidityjet;  // Rapidity of jet
  std::vector<double>         JET_emfjet;       // EMfraction of jet
  std::vector<double>         JET_chfjet;       // Charge fraction of jet (p/E)
  std::vector<double>         JET_massjet;      // Mass of jet
  std::vector<double>         JET_phijet;       // Phi of jet
  std::vector<double>         JET_etaetajet;    // Eta-Eta moment
  std::vector<double>         JET_etaphijet;    // Eta-Phi moment
  std::vector<double>         JET_phiphijet;    // Phi-Phi moment
  std::vector<double>         JET_metprjjet;    // MET projected in the direction to jet
					        
  std::vector<int>            JET_ntrk_tot;     // The number of total track
  std::vector<TLorentzVector> JET_P4trk_tot;    // Sum of 4-momentum of total track 
  std::vector<double>         JET_pt_tot;       // Pt of total track
  std::vector<double>         JET_trkmass_tot;  // Mass calculated from total track
  std::vector<double>         JET_sumPt_tot;    // 
  std::vector<double>         JET_ch_tot;       // Sum charge of total track
  					        
  std::vector<double>         JET_bDiscr;       
  std::vector<int>            JET_ntrk_tag;     // The number of SVX track
  std::vector<TLorentzVector> JET_P4trk_tag;    // Sum of 4-momentum of SVX track 
  std::vector<double>         JET_pt_tag;       // Pt of SVX track
  std::vector<double>         JET_trkmass_tag;  // Mass calculated from SVX track
  std::vector<double>         JET_ch_tag;       // Sum charge of SVX track
  std::vector<double>         JET_lxy_tag;      // Lxy 
  std::vector<TVector3>       JET_L3d_tag;      // 3 dimensional B hadron flight direction
  std::vector<int>            JET_tag;          // tag -1:negatively tag
                                                //      0:taggable
                                                //      1:positively ta
  std::vector<TVector3>       JET_SecVtxVertex;  // secondary vertex

  std::vector<int>            JET_nele;         // The number of electron into jet
  std::vector<int>            JET_nmuon;        // The number of muon into jet
  std::vector<TLorentzVector> JET_P4ele1;       // 4-momentum of electron 1st into jet
  std::vector<TLorentzVector> JET_P4ele2;       // 4-momentum of electron 2nd into jet
  std::vector<double>         JET_pt_ele1;      // Pt of first electron into jet
  std::vector<double>         JET_pt_ele2;      // Pt of second electron into jet
  std::vector<double>         JET_ptr_ele1;     // Ptrel of first electron into jet
  std::vector<double>         JET_ptr_ele2;     // Ptrel of second electron into jet
  std::vector<TLorentzVector> JET_P4muo1;       // 4-momentum of muon 1st into jet
  std::vector<TLorentzVector> JET_P4muo2;       // 4-momentum of muon 2nd into jet
  std::vector<double>         JET_pt_mu1;       // Pt of first muon into jet
  std::vector<double>         JET_pt_mu2;       // Pt of second muon into jet
  std::vector<double>         JET_ptr_mu1;      // Ptrel of first muon into jet
  std::vector<double>         JET_ptr_mu2;      // Ptrel of second muon into jet
  					        
  std::vector<int>            JET_nbqrk;        // 
  std::vector<double>         JET_dr_bqrk;      // Delta R between jet and nearest b quark
  std::vector<TLorentzVector> JET_P4bqrk;       // 4-momentum of nearest b quark (pttrue)
  std::vector<double>         JET_pt_bqrk;      // Pt of nearest b quark
  std::vector<int>            JET_assjet;       // The relation between quark and jet 0:not match, 
                                                //                                    1:first b-quark matched,
                                                //                                    2:second b-quark matched    
  std::vector<double>         JET_pt_bhad;      // Estimated b-hadron pt
  std::vector<double>         JET_pt_neu;       // neutral Pt of B hadron


  // ZHAD
  // ----
  std::vector<TLorentzVector> ZHAD_P4Zhadl0;             // 4-momentum of hadronic Z for Jet correction L0
  std::vector<TLorentzVector> ZHAD_P4Zhadl23;            // 4-momentum of hadronic Z for jet correction L2+L3
  std::vector<double>         ZHAD_etZhadl0;             // Et of hadronic Z for jet Correction L0
  std::vector<double>         ZHAD_etZhadl23;            // Et of hadronic Z for jet Correction L2+L3
  std::vector<double>         ZHAD_etaZhad;              // Eta of hadronic Z
  std::vector<double>         ZHAD_detetaZhad;           // Detector eta of hadronic Z
  std::vector<double>         ZHAD_rapidityZhad;         // Rapidity of hadronic Z
  std::vector<double>         ZHAD_emfZhad;              // EMfraction of hadronic Z
  std::vector<double>         ZHAD_chfZhad;              // Charge fraction of hadronic Z (p/E)
  std::vector<double>         ZHAD_phiZhad;              // Phi of hadronic Z
  std::vector<double>         ZHAD_etaetaZhad;           // Eta-Eta moment
  std::vector<double>         ZHAD_etaphiZhad;           // Eta-Phi moment
  std::vector<double>         ZHAD_phiphiZhad;           // Phi-Phi moment
  std::vector<double>         ZHAD_metprjZhad;           // MET projected in the direction to hadronic Z
  std::vector<double>         ZHAD_massZhad_l0;          // invariant mass corrected L0
  std::vector<double>         ZHAD_massZhad_l23;         // invariant mass corrected L2+L3
  std::vector<double>         ZHAD_collinearityZhad_l0;  //
  std::vector<double>         ZHAD_collinearityZhad_l23; //

  
  // MUON
  // ----
  std::vector<TLorentzVector> MUON_P4muon;              // 4-momentum of muon
  std::vector<double>         MUON_isolMuonSumPt;       // sum of tracks Pt associated to isolated muon
  std::vector<int>            MUON_isolMuonTrkNumber;   // number of tracks associated to isolated muon


  // ELECTRON
  // --------
  std::vector<TLorentzVector> ELECTRON_P4electron;         // 4-momentum of electron
  std::vector<double>         ELECTRON_isolEleSumPt;       // sum of tracks Pt associated to isolated electron
  std::vector<int>            ELECTRON_isolEleTrkNumber;   // number of tracks associated to isolated electron
  std::vector<int>            ELECTRON_eleId;              // electron ID


void loadSimpleNtuple(int index,TChain* chain,
		      TClonesArray* evtClass,
		      TClonesArray* jetClass,
		      TClonesArray* muonClass,
		      TClonesArray* electronClass,
		      TClonesArray* zhadClass);

