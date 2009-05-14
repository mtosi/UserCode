#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "Riostream.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/findObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/THStackLegend.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ListFashionAttributedHisto.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"

using namespace vbfhzz2l2b;

void plotEfficiencyCutsVSmass(TString outputFileName = "plotEfficiencyVSmassOutput.root") {

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
  std::vector<TString> FileSuffix;
  FileSuffix.push_back("VBFH160ZZ");
  FileSuffix.push_back("VBFH200ZZ");
  FileSuffix.push_back("VBFH400ZZ");
  FileSuffix.push_back("VBFH800ZZ");
  FileSuffix.push_back("ZZ0Jets");
  FileSuffix.push_back("ZZ1Jets");
  FileSuffix.push_back("ZZ2Jets");
  FileSuffix.push_back("Zbb0Jets");
  FileSuffix.push_back("Zbb1Jets");
  FileSuffix.push_back("Zbb2Jets");
  // ************************************************************
  // List of Files
  TList * FileList = new TList();
  std::vector<TString>::const_iterator FileSuffix_itr = FileSuffix.begin();
  for ( ; FileSuffix_itr != FileSuffix.end(); ++FileSuffix_itr )
    FileList->Add( TFile::Open(*FileSuffix_itr+"output.root"  ) );
    // ************************************************************
  // List of xSec (in pb)
  // stored in 
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"
  std::vector<double> crossSections;
  crossSections.push_back(1000*xSec_VBFH160_*BR_H160ZZ_);
  crossSections.push_back(1000*xSec_VBFH200_*BR_H200ZZ_);
  crossSections.push_back(1000*xSec_VBFH400_*BR_H400ZZ_);
  crossSections.push_back(1000*xSec_VBFH800_*BR_H800ZZ_);
  crossSections.push_back(1000*xSec_ZZ0Jets_);
  crossSections.push_back(1000*xSec_ZZ1Jets_);
  crossSections.push_back(1000*xSec_ZZ2Jets_);
  crossSections.push_back(1000*xSec_Zbb0Jets_);
  crossSections.push_back(1000*xSec_Zbb1Jets_);
  crossSections.push_back(1000*xSec_Zbb2Jets_);
  // *************************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("recHlljj"                                );
  varNameVector.push_back("recHlljj_leptonicZcutsAFTERZdileptonMass");
//  varNameVector.push_back("recHlljj_leptonicZcutsAFTERZleptonDeltaR");
//  varNameVector.push_back("recHlljj_hadronicZcutsAFTERZdijetMass"   );
//  varNameVector.push_back("recHlljj_hadronicZcutsAFTERZjetDeltaR"   );
//  varNameVector.push_back("recHlljj_FWDcutsAFTERleptonsShiftedEta"  );
//  varNameVector.push_back("recHlljj_FWDcutsAFTERjetsShiftedEta"     );
//  varNameVector.push_back("recHlljj_FWDcutsAFTERforwardDiJetMass"   );
//  varNameVector.push_back("recHlljj_FWDcutsAFTERpTbalance"          );
//  varNameVector.push_back("recHlljj_final"                          );
  // ************************************************************
  std::vector<TString> varNameCutVector;
  varNameCutVector.push_back("ZdileptonMass"    );
//  varNameCutVector.push_back("ZleptonDeltaR"    );
//  varNameCutVector.push_back("ZdijetMass"       );
//  varNameCutVector.push_back("ZjetDeltaR"       );
//  varNameCutVector.push_back("leptonsShiftedEta");
//  varNameCutVector.push_back("jetsShiftedEta"   );
//  varNameCutVector.push_back("forwardDiJetMass" );
//  varNameCutVector.push_back("pTbalance"        );
//  varNameCutVector.push_back("final"            );
  // ************************************************************

  std::vector<TH1*> eventsNumberTH1vector;
  findObj( eventsNumberTH1vector, 
	   FileList, 
	   "eventsNumber", 
	   &FileSuffix
	   );
  int indexSample = 0;
  std::vector<TH1*>::const_iterator eventsNumberTH1vector_itr = eventsNumberTH1vector.begin();
  for ( ; eventsNumberTH1vector_itr != eventsNumberTH1vector.end(); ++eventsNumberTH1vector_itr,
	                                                            indexSample++                ) {
    double sampleEvents = (*eventsNumberTH1vector_itr)->GetMaximum();
    if(sampleEvents != 0.) crossSections[indexSample]=crossSections[indexSample]/sampleEvents;
  }

  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );
  int indexCut = 0;
  std::vector<TString>::const_iterator varNameVector_itr;
  std::vector<TString>::const_iterator varNameCutVector_itr = varNameCutVector.begin();
  std::cout << "varNameCutVector.size(): " << varNameCutVector.size() << std::endl;
  for ( ; varNameCutVector_itr != varNameCutVector.end(); ++varNameCutVector_itr,
	                                                  indexCut++ ) {
    std::cout << "varNameCutVector_itr: " << *varNameCutVector_itr << std::endl;
    std::vector<TH1*> afterTH1vector;
    std::vector<TH1*> beforeTH1vector;
    varNameVector_itr = varNameVector.begin();
    for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
      TString varName = *varNameVector_itr;
      if (varName.Contains(*varNameCutVector_itr)) {
	std::cout << "varName: " << varName << std::endl;
	findObj( afterTH1vector,
		 FileList, 
		 *varNameVector_itr, 
		 &FileSuffix
		 );
	findObj( beforeTH1vector,
		 FileList, 
		 *(varNameVector_itr-1), 
		 &FileSuffix
		 );
      }
    }
    std::cout << "beforeTH1vector.size: " << beforeTH1vector.size() << std::endl;
    std::cout << "afterTH1vector.size: " << afterTH1vector.size() << std::endl;
    // THStackLegend histogram
    TString Stack( "Stack_" );
    THStackLegend<TH1> * StackLegend_ = new THStackLegend<TH1>( *varNameCutVector_itr, 0.75,0.15,0.98,0.4 );
    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameCutVector_itr, *varNameCutVector_itr );
    directory->cd();

    int    nbins = beforeTH1vector[0]->GetNbinsX();
    double xmin  = beforeTH1vector[0]->GetXaxis()->GetXmin();
    double xmax  = beforeTH1vector[0]->GetXaxis()->GetXmax();
    FileSuffix_itr = FileSuffix.begin();
    indexSample = 0;
    for ( ; FileSuffix_itr != FileSuffix.end(); ++FileSuffix_itr,
	                                        indexSample++     ) {
      TH1D * histo = new TH1D(*varNameCutVector_itr+"efficiency",
			      *varNameCutVector_itr+" efficiency",
			      nbins,xmin,xmax);
      std::cout << "FileSuffix_itr: " << *FileSuffix_itr << std::endl;

      for ( int ibin = 1; ibin < nbins; ibin++ ) {
	double eventsBefore = beforeTH1vector[indexSample]->GetBinContent(ibin);
	double eventsAfter  = afterTH1vector[indexSample]->GetBinContent(ibin);
	double efficiency = 0.;
	//	if( eventsAfter >= 5. ) {
	if( eventsBefore != 0. ) {
	  efficiency = eventsAfter/eventsBefore;
	  histo->SetBinContent(ibin,efficiency);
	  histo->SetBinError(ibin,TMath::Sqrt(efficiency*(1-efficiency)/eventsBefore));
	}
      }
      std::cout << "varNameCutVector_itr: " << *varNameCutVector_itr << std::endl;
      FashionAttributedHisto<TH1D> * dressedHisto;
      bool dressed = false;
      if (      *FileSuffix_itr == "VBFH160ZZ" ) { dressedHisto = new H160<TH1D>((TH1D*)histo);     dressed = true; }
      else if ( *FileSuffix_itr == "VBFH200ZZ" ) { dressedHisto = new H200<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "VBFH400ZZ" ) { dressedHisto = new H400<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "VBFH800ZZ" ) { dressedHisto = new H800<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "ZZ0Jets"   ) { dressedHisto = new ZZ0jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "ZZ1Jets"   ) { dressedHisto = new ZZ1jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "ZZ2Jets"   ) { dressedHisto = new ZZ2jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "Zbb0Jets"  ) { dressedHisto = new Zbb0jets<TH1D>((TH1D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "Zbb1Jets"  ) { dressedHisto = new Zbb1jets<TH1D>((TH1D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "Zbb2Jets"  ) { dressedHisto = new Zbb2jets<TH1D>((TH1D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "WZ0Jets"   ) { dressedHisto = new WZ0jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "WZ1Jets"   ) { dressedHisto = new WZ1jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "WZ2Jets"   ) { dressedHisto = new WZ2jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt0Jets"   ) { dressedHisto = new tt0jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt1Jets"   ) { dressedHisto = new tt1jets<TH1D>((TH1D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt2Jets"   ) { dressedHisto = new tt2jets<TH1D>((TH1D*)histo);  dressed = true; }
      else dressedHisto = new FashionAttributedHisto<TH1D>((TH1D*)histo); 
      StackLegend_->Add(dressedHisto,*FileSuffix_itr,false,"p",false,dressed);
      
      dressedHisto->Write();
    }
    //    StackLegend_->SetLogY();
    StackLegend_->Print("nostackp",*varNameCutVector_itr+"EfficiencyVSmassPlot.jpg");
    StackLegend_->Write("nostack");
    StackLegend_->Draw("nostack");
    
  }
    
}
