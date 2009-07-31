#ifndef SIMPLENTPLEOBJ_H
#define SIMPLENTPLEOBJ_H


#include "TObject.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

namespace vbfhzz2l2b {

  namespace SimpleNtpleObj {

    class EVT;
    class JET;
    class MUON;
    class ELECTRON;
    class ZHAD;

    // event variables
    class EVT : public TObject {
    public:
      EVT();
      ~EVT() {};
      //+++
      int    Run;                   // run number
      int    Event;                 // event number
      double Ilum;                  // instantaneous luminosity (e30)
      int    eventID;               // event ID
                                    // VBF:123 or 124
                                    // ggF: 102
      int    nPV;                   // number of primary vertices
      int    trigpath;              // Z_BB Trigger Path: 1*main + 10*test1 + 100*test2 (if exists). 
      // ex. 101 means main + 2nd test trigger where fired.
      double  pthat;                // pthat
      
      TLorentzVector P4bquark1;     // first  b-quark
      TLorentzVector P4bquark2;     // second b-quark
      int            indjetb1;      // first b-quark associated jet index
      int            indjetb2;      // second b-quark associated jet index
      
      int      njet;                // number of jets in the event
      int      ngoodjet;            // number of "good" jets in the event
      int      nbtag;               // number of tag in the event
      double   Zvertex;             // Z of the reconstructed primary vertex
      TVector2 P2met;               // Missing Et vector
      
      int nmuon;
      int nelectron;
      int nZhad;

      ClassDef(EVT,1)
    };

    // JET
    //-----
    class JET : public TObject {
    public:
      JET();
      ~JET() {};
      TLorentzVector P4jetl0;     // 4-momentum of jet for Jet correction L0
      TLorentzVector P4jetl23;    // 4-momentum of jet for Jet correction L2+L3
      double         etjetl0;     // Et of jet for Jet Correction L0
      double         etjetl23;    // Et of jet for Jet Correction L2+L3
      double         etajet;      // Eta of jet
      double         detetajet;   // Detector eta of jet
      double         rapidityjet; // Rapidity of jet
      double         emfjet;      // EMfraction of jet
      double         chfjet;      // Charge fraction of jet (p/E)
      double         massjet;     // Mass of jet
      double         phijet;      // Phi of jet
      double         etaetajet;   // Eta-Eta moment
      double         etaphijet;   // Eta-Phi momentum
      double         phiphijet;   // Phi-Phi momentum
      double         metprjjet;   // MET projected in the direction to jet

      int            ntrk_tot;    // The number of total track
      TLorentzVector P4trk_tot;   // Sum of 4-momentum of total track 
      double         pt_tot;      // Pt of total track
      double         trkmass_tot; // Mass calculated from total track
      double         sumPt_tot;   // 
      double         ch_tot;      // Sum charge of total track
      
      double         bDiscr;
      int            ntrk_tag;    // The number of SVX track
      TLorentzVector P4trk_tag;   // Sum of 4-momentum of SVX track 
      double         pt_tag;      // Pt of SVX track
      double         trkmass_tag; // Mass calculated from SVX track
      double         ch_tag;      // Sum charge of SVX track
      double         lxy_tag;     // Lxy 
      TVector3       L3d_tag;     // 3 dimensional B hadron flight direction
      int            tag;         // tag -1:negatively tag
                                  //      0:taggable
                                  //      1:positively ta
      TVector3   SecVtxVertex;    // secondary vertex
      
      int            nele;        // The number of electron into jet
      int            nmuon;       // The number of muon into jet
      TLorentzVector P4ele1;      // 4-momentum of electron 1st into jet
      TLorentzVector P4ele2;      // 4-momentum of electron 2nd into jet
      double         pt_ele1;     // Pt of first electron into jet
      double         pt_ele2;     // Pt of second electron into jet
      double         ptr_ele1;    // Ptrel of first electron into jet
      double         ptr_ele2;    // Ptrel of second electron into jet
      TLorentzVector P4muo1;      // 4-momentum of muon 1st into jet
      TLorentzVector P4muo2;      // 4-momentum of muon 2nd into jet
      double         pt_mu1;      // Pt of first muon into jet
      double         pt_mu2;      // Pt of second muon into jet
      double         ptr_mu1;     // Ptrel of first muon into jet
      double         ptr_mu2;     // Ptrel of second muon into jet
      
      int            nbqrk;       // 
      double         dr_bqrk;     // Delta R between jet and nearest b quark
      TLorentzVector P4bqrk;      // 4-momentum of nearest b quark (pttrue)
      double         pt_bqrk;     // Pt of nearest b quark
      int            assjet;      // The relation between quark and jet 0:not match, 
                                  //                                    1:first b-quark matched,
                                  //                                    2:second b-quark matched    
      double         pt_bhad;     // Estimated b-hadron pt
      double         pt_neu;      // neutral Pt of B hadron
      
      //  ClassDef(JET,1);
    
    };

// ZHAD
class ZHAD : public TObject {
public:
  ZHAD();
  ~ZHAD() {};
  TLorentzVector P4Zhadl0;             // 4-momentum of hadronic Z for Jet correction L0
  TLorentzVector P4Zhadl23;            // 4-momentum of hadronic Z for jet correction L2+L3
  double         etZhadl0;             // Et of hadronic Z for jet Correction L0
  double         etZhadl23;            // Et of hadronic Z for jet Correction L2+L3
  double         etaZhad;              // Eta of hadronic Z
  double         detetaZhad;           // Detector eta of hadronic Z
  double         rapidityZhad;         // Rapidity of hadronic Z
  double         emfZhad;              // EMfraction of hadronic Z
  double         chfZhad;              // Charge fraction of hadronic Z (p/E)
  double         phiZhad;              // Phi of hadronic Z
  double         etaetaZhad;           // Eta-Eta moment
  double         etaphiZhad;           // Eta-Phi moment
  double         phiphiZhad;           // Phi-Phi moment
  double         metprjZhad;           // MET projected in the direction to hadronic Z
  double         massZhad_l0;          // invariant mass corrected L0
  double         massZhad_l23;         // invariant mass corrected L2+L3
  double         collinearityZhad_l0;  //
  double         collinearityZhad_l23; //

  //  ClassDef(ZHAD,1);
};

// MUON
class MUON : public TObject {
public:
  MUON();
  ~MUON() {};
  TLorentzVector P4muon;            // 4-momentum of muon
  double         isolMuonSumPt;     // sum of tracks Pt associated to isolated muon
  int            isolMuonTrkNumber; // number of tracks associated to isolated muon

  //  ClassDef(MUON,1);
};

// ELECTRON
class ELECTRON : public TObject {
public:
  ELECTRON();
  ~ELECTRON() {};
  TLorentzVector P4electron;         // 4-momentum of electron
  double         isolEleSumPt;       // sum of tracks Pt associated to isolated electron
  int            isolEleTrkNumber;   // number of tracks associated to isolated electron
  int            eleId;              // electron ID

  //  ClassDef(ELECTRON,1);
};
  }
}
#endif // SIMPLENTPLEOBJ_H
