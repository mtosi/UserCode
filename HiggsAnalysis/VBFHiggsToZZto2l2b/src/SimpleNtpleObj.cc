#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleNtpleObj.h"

// Event constructor
// -----------------
EVT::EVT() :
  Run(0),                 // run number
  Event(0),               // event number
  Ilum(0.),               // instantaneous luminosity (e30)
  eventID(0),             // event ID
                          // VBF:123 or 124
                          // ggF: 102
  nPV(0),                 // Number of vertices 
  trigpath(-1),           // Z_BB Trigger Path: 1*main + 10*test1 + 100*test2 (if exists).
                          // ex. 101 means main + 2nd test trigger where fired.
  pthat(0),               // pthat
  P4bquark1(0.,0.,0.,0.), // first  b-quark
  P4bquark2(0.,0.,0.,0.), // second b-quark
  indjetb1(-1),           // first b-quark associated jet index
  indjetb2(-1),           // second b-quark associated jet index
  njet(0),                // number of jets in the event
  ngoodjet(0),            // number of "good" jets in the event
  nbtag(0),               // number of tag in the event
  Zvertex(0.),            // Z of the reconstructed primary vertex
  P2met(0.,0.),           // Missing Et vector
  nmuon(0),
  nelectron(0),
  nZhad(0)

{
}

// Jet
// ---
JET::JET() :
  P4jetl0(0.,0.,0.,0.),   // 4-momentum of jet for Jet correction L0
  P4jetl23(0.,0.,0.,0.),  // 4-momentum of jet for Jet correction L2+L3
  etjetl0(0.),            // Et of jet for Jet Correction L0       
  etjetl23(0.),           // Et of jet for Jet Correction L2+L3
  etajet(0.),             // Eta of jet
  detetajet(0.),          // Detector eta of jet
  rapidityjet(0.),        // Rapidity of jet
  emfjet(0.),             // EMfraction of jet
  chfjet(0.),             // Charge fraction of jet (p/E)
  massjet(0.),            // Mass of jet
  phijet(0.),             // Phi of jet
  etaetajet(0.),          // Eta-Eta moment
  etaphijet(0.),          // Eta-Phi moment
  phiphijet(0.),          // Phi-Phi moment
  metprjjet(0.),          // MET projected in the direction to jet
  ntrk_tot(0),            // The number of total track
  P4trk_tot(0.,0.,0.,0.), // Sum of 4-momentum of total track 
  pt_tot(0.),             // Pt of total track
  trkmass_tot(0.),        // Mass calculated from total track
  sumPt_tot(0.),          // Pt of total track
  ch_tot(0.),             // Sum charge of total track    
  bDiscr(0.),
  ntrk_tag(0),            // The number of SVX track
  P4trk_tag(0.,0.,0.,0.), // Sum of 4-momentum of SVX track 
  pt_tag(0.),             // Pt of SVX track
  trkmass_tag(0.),        // Mass calculated from SVX track
  ch_tag(0.),             // Sum charge of SVX track
  lxy_tag(0.),            // Lxy 
  L3d_tag(0.,0.,0.),      // 3 dimensional B hadron flight direction
  tag(0),                 // tag -1:negatively tag
                          //      0:taggable
                          //      1:positively tag
  SecVtxVertex(0.,0.,0.), // secondary vertex
  nele(0),                // The number of electron into jet
  nmuon(0),               // The number of muon into jet
  P4ele1(0.,0.,0.,0.),    // 4-momentum of electron 1st into jet
  P4ele2(0.,0.,0.,0.),    // 4-momentum of electron 2nd into jet
  pt_ele1(0.),            // Pt of first muon into jet
  pt_ele2(0.),            // Pt of second muon into jet 
  ptr_ele1(0.),           // Ptrel of first muon into jet
  ptr_ele2(0.),           // Ptrel of second muon into jet 
  P4muo1(0.,0.,0.,0.),    // 4-momentum of muon 1st into jet
  P4muo2(0.,0.,0.,0.),    // 4-momentum of muon 2nd into jet
  pt_mu1(0.),             // Pt of first muon into jet
  pt_mu2(0.),             // Pt of second muon into jet 
  ptr_mu1(0.),            // Ptrel of first muon into jet
  ptr_mu2(0.),            // Ptrel of second muon into jet 
  nbqrk(0),               // 
  dr_bqrk(0.),            // Delta R between jet and nearest b quark
  P4bqrk(0.,0.,0.,0.),    // 4-momentum of nearest b quark (pttrue)
  pt_bqrk(0.),            // Pt of of nearest b quark
  assjet(0),              // The relation between quark and jet 0:not match, 
                          //                                    1:first b-quark matched,
                          //                                    2:second b-quark matched    
  pt_bhad(0.),            // Estimated b-hadron pt
  pt_neu(0.)              // neutral Pt of B hadron
{
}

// ZHAD
ZHAD::ZHAD() :
  P4Zhadl0(0.,0.,0.,0.),    // 4-momentum of hadronic Z for Jet correction L0
  P4Zhadl23(0.,0.,0.,0.),   // 4-momentum of hadronic Z for jet correction L2+L3
  etZhadl0(0.),             // Et of hadronic Z for jet Correction L0
  etZhadl23(0.),            // Et of hadronic Z for jet Correction L2+L3
  etaZhad(0.),              // Eta of hadronic Z
  detetaZhad(0.),           // Detector eta of hadronic Z
  rapidityZhad(0.),         // Rapidity of hadronic Z
  emfZhad(0.),              // EMfraction of hadronic Z
  chfZhad(0.),              // Charge fraction of hadronic Z (p/E)
  phiZhad(0.),              // Phi of hadronic Z
  etaetaZhad(0.),           // Eta-Eta moment
  etaphiZhad(0.),           // Eta-Phi moment
  phiphiZhad(0.),           // Phi-Phi moment
  metprjZhad(0.),           // MET projected in the direction to hadronic Z
  massZhad_l0(0.),          // invariant mass corrected L0
  massZhad_l23(0.),         // invariant mass corrected L2+L3
  collinearityZhad_l0 (0.), //
  collinearityZhad_l23(0.)  //
{
}

// MUON
MUON::MUON() :
  P4muon(0.,0.,0.,0.), // 4-momentum of muon
  isolMuonSumPt(0.),   // sum of tracks Pt associated to isolated muon
  isolMuonTrkNumber(0) // number of tracks associated to isolated muon
{
}

// ELECTRON
ELECTRON::ELECTRON() :
  P4electron(0.,0.,0.,0.), // 4-momentum of electron
  isolEleSumPt(0.),        // sum of tracks Pt associated to isolated electron
  isolEleTrkNumber(0),     // number of tracks associated to isolated electron
  eleId(0)                 // electron ID
{
}

// MET
MET::MET() :
  P4met(0.,0.,0.,0.) // 4-momentum of met
{
}
