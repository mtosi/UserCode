/**

  This macro will draw histograms multiplied by their cross-section
  from a list of root files and write them to a target root file.
  The target file is newly created and must not be
  identical to one of the source files.

 */


#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "Riostream.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/findObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/mergeObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/THStackLegend.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ListFashionAttributedHisto.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"

using namespace vbfhzz2l2b;

void plotMergedHisto(TString outputFileName = "plotMergedOutput.root") {

  // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

  // ************************************************************
  std::vector<TString> sgnFileSuffix;
  sgnFileSuffix.push_back("VBFH160ZZ");
  sgnFileSuffix.push_back("VBFH200ZZ");
  sgnFileSuffix.push_back("VBFH400ZZ");
  sgnFileSuffix.push_back("VBFH800ZZ");
  std::vector<TString> bkgZZFileSuffix;
  bkgZZFileSuffix.push_back("ZZ0Jets");
  bkgZZFileSuffix.push_back("ZZ1Jets");
  bkgZZFileSuffix.push_back("ZZ2Jets");
  std::vector<TString> bkgZbbFileSuffix;
  bkgZbbFileSuffix.push_back("Zbb0Jets");
  bkgZbbFileSuffix.push_back("Zbb1Jets");
  bkgZbbFileSuffix.push_back("Zbb2Jets");
  // ************************************************************
  // List of Files
  TList * sgnFileList    = new TList();
  std::vector<TString>::const_iterator sgnFileSuffix_itr = sgnFileSuffix.begin();
  for ( ; sgnFileSuffix_itr != sgnFileSuffix.end(); ++sgnFileSuffix_itr )
    sgnFileList->Add( TFile::Open(*sgnFileSuffix_itr+"output.root"  ) );
  TList * bkgZZFileList  = new TList();
  std::vector<TString>::const_iterator bkgZZFileSuffix_itr = bkgZZFileSuffix.begin();
  for ( ; bkgZZFileSuffix_itr != bkgZZFileSuffix.end(); ++bkgZZFileSuffix_itr )
    bkgZZFileList->Add( TFile::Open(*bkgZZFileSuffix_itr+"output.root"  ) );
  TList * bkgZbbFileList = new TList();
  std::vector<TString>::const_iterator bkgZbbFileSuffix_itr = bkgZbbFileSuffix.begin();
  for ( ; bkgZbbFileSuffix_itr != bkgZbbFileSuffix.end(); ++bkgZbbFileSuffix_itr )
    bkgZbbFileList->Add( TFile::Open(*bkgZbbFileSuffix_itr+"output.root"  ) );
    // ************************************************************
  // List of xSec (in pb)
  // stored in 
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"
  std::vector<double> sgnCrossSections;
  std::vector<double> bkgZZCrossSections;
  std::vector<double> bkgZbbCrossSections;
  sgnCrossSections.push_back(1000*xSec_VBFH160_*BR_H160ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH200_*BR_H200ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH400_*BR_H400ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH800_*BR_H800ZZ_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ0Jets_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ1Jets_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ2Jets_);
  bkgZbbCrossSections.push_back(1000*xSec_Zbb0Jets_);
  bkgZbbCrossSections.push_back(1000*xSec_Zbb1Jets_);
  bkgZbbCrossSections.push_back(1000*xSec_Zbb2Jets_);
  // *************************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("recZllMass_leptonicZcutsAFTERZdileptonMass"      );
  varNameVector.push_back("recZllMass_leptonicZcutsAFTERZleptonDeltaR"      );
  varNameVector.push_back("recZllDeltaR"                                    );
  varNameVector.push_back("recJetsInvMass"                                  );
  varNameVector.push_back("recZjjMass_hadronicZcutsAFTERZdijetMass"         );
  varNameVector.push_back("recZjjMass_hadronicZcutsAFTERZjetDeltaR"         );
  varNameVector.push_back("recZjjDeltaR"                                    );
  varNameVector.push_back("forwardDiJetMass"                                );
  varNameVector.push_back("recFWDjj_hadronicZcutsAFTERjetDeltaR"            );
  varNameVector.push_back("recZllShiftedEta"                                );
  varNameVector.push_back("recZleptonsShitedEta_hadronicZcutsAFTERjetDeltaR");
  varNameVector.push_back("pTbalance"                                       );
  varNameVector.push_back("recFWDjjPtBalance_hadronicZcutsAFTERjetDeltaR"   );
  varNameVector.push_back("recHZZdeltaR"                                    );
  varNameVector.push_back("recHZZdeltaR_FWDcutsAFTERpTbalance"              );
  varNameVector.push_back("recHlljj"                                        );
  varNameVector.push_back("recHlljj_final"                                  );
  varNameVector.push_back("recHlljj_forwardDiJetMass"                       );
  // *************************************************************************

  int indexSample = 0;
  std::vector<TH1*> sgnEventsNumberTH1vector;
  if (sgnFileSuffix.size() != 0) {
    findObj( sgnEventsNumberTH1vector, 
	     sgnFileList, 
	     "eventsNumber", 
	     &sgnFileSuffix
	     );
    indexSample = 0;
    std::vector<TH1*>::const_iterator sgnEventsNumberTH1vector_itr = sgnEventsNumberTH1vector.begin();  
    for ( ; sgnEventsNumberTH1vector_itr != sgnEventsNumberTH1vector.end(); ++sgnEventsNumberTH1vector_itr,
	    indexSample++) { 
      double sampleEvents = (*sgnEventsNumberTH1vector_itr)->GetMaximum(); 
      if(sampleEvents!=0.) sgnCrossSections[indexSample] = sgnCrossSections[indexSample]/sampleEvents;
    }
  }
  std::vector<TH1*> bkgZZEventsNumberTH1vector;
  findObj( bkgZZEventsNumberTH1vector, 
	   bkgZZFileList, 
	   "eventsNumber", 
	   &bkgZZFileSuffix
	   );
  std::vector<TH1*>::const_iterator bkgZZEventsNumberTH1vector_itr = bkgZZEventsNumberTH1vector.begin();  
  indexSample = 0;
  for ( ; bkgZZEventsNumberTH1vector_itr != bkgZZEventsNumberTH1vector.end(); ++bkgZZEventsNumberTH1vector_itr,
	                                                                      indexSample++) { 
    double sampleEvents = (*bkgZZEventsNumberTH1vector_itr)->GetMaximum();
    bkgZZCrossSections[indexSample] = bkgZZCrossSections[indexSample]/sampleEvents;
  }
  std::vector<TH1*> bkgZbbEventsNumberTH1vector;
  findObj( bkgZbbEventsNumberTH1vector, 
	   bkgZbbFileList, 
	   "eventsNumber", 
	   &bkgZbbFileSuffix
	   );
  std::vector<TH1*>::const_iterator bkgZbbEventsNumberTH1vector_itr = bkgZbbEventsNumberTH1vector.begin();  
  indexSample = 0;
  for ( ; bkgZbbEventsNumberTH1vector_itr != bkgZbbEventsNumberTH1vector.end(); ++bkgZbbEventsNumberTH1vector_itr,
	                                                                        indexSample++) { 
    double sampleEvents = (*bkgZbbEventsNumberTH1vector_itr)->GetMaximum();
    bkgZbbCrossSections[indexSample] = bkgZbbCrossSections[indexSample]/sampleEvents;
  }


  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );

  std::vector<TString>::const_iterator varNameVector_itr = varNameVector.begin();
  for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
    std::vector<TH1*> TH1vector;
    std::vector<TH1*> bkgZZTH1vector;
    std::vector<TH1*> bkgZbbTH1vector;
    findObj( TH1vector, 
	     sgnFileList, 
	     *varNameVector_itr, 
	     &sgnFileSuffix,
	     "TH1",
	     &sgnCrossSections
	     );
    findObj( bkgZZTH1vector, 
	     bkgZZFileList, 
	     *varNameVector_itr, 
	     &bkgZZFileSuffix,
	     "TH1",
	     &bkgZZCrossSections
	     );
    findObj( bkgZbbTH1vector, 
	     bkgZbbFileList, 
	     *varNameVector_itr, 
	     &bkgZbbFileSuffix,
	     "TH1",
	     &bkgZbbCrossSections
	     );

    //    TH1vector.push_back(mergeObj((*varNameVector_itr)+"_ZZNjets", bkgZZTH1vector, &bkgZZCrossSections, "kTRUE","kTRUE"));
    //    TH1vector.push_back(mergeObj((*varNameVector_itr)+"_ZbbNjets",bkgZbbTH1vector,&bkgZbbCrossSections,"kTRUE","kTRUE"));
    TH1vector.push_back(mergeObj((*varNameVector_itr)+"_ZZNjets", bkgZZTH1vector));
    TH1vector.push_back(mergeObj((*varNameVector_itr)+"_ZbbNjets",bkgZbbTH1vector));

    std::cout << "TH1vector.size(): " << TH1vector.size() << std::endl;

    // THStackLegend histogram
    TString Stack( "Stack_" );
    THStackLegend<TH1> * StackLegend_ = new THStackLegend<TH1>( *varNameVector_itr, 0.6,0.68,0.98,0.97 );

    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameVector_itr, *varNameVector_itr );
    directory->cd();
    std::vector<TH1*>::const_iterator TH1vector_itr = TH1vector.begin();
    for ( ; TH1vector_itr != TH1vector.end(); ++TH1vector_itr ) {
      TH1 * histo = (TH1*)(*TH1vector_itr)->Clone();
      TString histoName = histo->GetName();
      TString sample = histoName.Remove(0,(*varNameVector_itr).Length()+1);
      FashionAttributedHisto<TH1D> * dressedHisto;
      bool dressed = false;
      if (      sample.Contains("VBFH160ZZ") ) { dressedHisto = new H160<TH1D>((TH1D*)histo);       dressed = true; }
      else if ( sample.Contains("VBFH200ZZ") ) { dressedHisto = new H200<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( sample.Contains("VBFH400ZZ") ) { dressedHisto = new H400<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( sample.Contains("VBFH800ZZ") ) { dressedHisto = new H800<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( sample.Contains("ZZ" ) && sample.Contains("jets") ) { dressedHisto = new ZZNjets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( sample.Contains("Zbb") && sample.Contains("jets") ) { dressedHisto = new ZbbNjets<TH1D>((TH1D*)histo); dressed = true; }
      else if ( sample.Contains("WZ" ) && sample.Contains("jets") ) { dressedHisto = new WZNjets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( sample.Contains("tt" ) && sample.Contains("jets") ) { dressedHisto = new ttNjets<TH1D>((TH1D*)histo);  dressed = true; }
      else dressedHisto = new FashionAttributedHisto<TH1D>((TH1D*)histo);
      char nev[50];
      sprintf(nev,sample+": %.2f",dressedHisto->Integral());
      StackLegend_->Add(dressedHisto,nev,false,"lf",false,dressed);
      //      StackLegend_->Add(dressedHisto,sample,false,"fl",false,dressed);

      dressedHisto->Write();
    }
    StackLegend_->SetLogY();
    StackLegend_->Print("nostack",*varNameVector_itr+"MergedMultiPlot.jpg");
    StackLegend_->Write("nostackhisto");
    StackLegend_->Draw("nostackhisto");
    //    StackLegend_->SavePrimitive(Stack+*varNameVector_itr+".C","nostackhisto");

    //    TH1vector.clear();
    //    bkgZZTH1vector.clear();
    //    bkgZbbTH1vector.clear();
  }

}


