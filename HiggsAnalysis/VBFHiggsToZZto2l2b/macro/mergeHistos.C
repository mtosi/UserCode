/**

  This macro will add histograms multiplied by their cross-section
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

#include "AnalysisExamples/AnalysisClasses/interface/findObj.h"
#include "AnalysisExamples/AnalysisClasses/interface/mergeObj.h"
#include "AnalysisExamples/AnalysisClasses/interface/THStackLegend.h"
#include "AnalysisExamples/AnalysisObjects/interface/xSecLO.h"
#include "AnalysisExamples/AnalysisObjects/interface/xSecNLO.h"
#include "AnalysisExamples/AnalysisObjects/interface/BR.h"

using namespace anaobj;

void merge(TString outputFileName = "mergeOutput.root") {

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
  //  FileSuffix.push_back("VBFH160ZZ");
  //  FileSuffix.push_back("VBFH200ZZ");
  //  FileSuffix.push_back("VBFH400ZZ");
  //  FileSuffix.push_back("VBFH800ZZ");
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
  for ( ; FileSuffix_itr != FileSuffix.end(); ++FileSuffix_itr ) {
    FileList->Add( TFile::Open(*FileSuffix_itr+"output.root"  ) );
  }
  // ************************************************************
  // List of xSec (in pb)
  // stored in 
  // "AnalysisExamples/AnalysisObjects/interface/xSecLO.h"
  // "AnalysisExamples/AnalysisObjects/interface/xSecNLO.h"
  // "AnalysisExamples/AnalysisObjects/interface/BR.h"
  std::vector<double> crossSections;
  //  crossSections.push_back(xSec_VBFH160_*BR_H160ZZ_);
  //  crossSections.push_back(xSec_VBFH200_*BR_H200ZZ_);
  //  crossSections.push_back(xSec_VBFH400_*BR_H400ZZ_);
  //  crossSections.push_back(xSec_VBFH800_*BR_H800ZZ_);
  crossSections.push_back(xSec_ZZ0Jets_);
  crossSections.push_back(xSec_ZZ1Jets_);
  crossSections.push_back(xSec_ZZ2Jets_);
  crossSections.push_back(xSec_Zbb0Jets_*BR_Zll_);
  crossSections.push_back(xSec_Zbb1Jets_*BR_Zll_);
  crossSections.push_back(xSec_Zbb2Jets_*BR_Zll_);
  // ************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("recHlljj");
  varNameVector.push_back("pTbalance");
  //varNameVector.push_back("recZllMass");
  //varNameVector.push_back("recZjjMass");
  //varNameVector.push_back("recZllMassResolution");
  //varNameVector.push_back("recZjjMassResolution");
  varNameVector.push_back("recHZZdeltaPhi");
  varNameVector.push_back("recHZZdeltaEta");
  varNameVector.push_back("recHZZdeltaR");
  varNameVector.push_back("DeltaPhiZbbbar");
  varNameVector.push_back("DeltaPhiZelectrons");
  varNameVector.push_back("DeltaPhiZmuons");
  // ************************************************************
  
  std::vector<TH1*> TH1vector;
  std::vector<TH1*> eventsNumberTH1vector;
  
  findObj( eventsNumberTH1vector, 
	   FileList, 
	   "eventsNumber", 
	   &FileSuffix
	   );
  
  std::vector<double>::const_iterator crossSections_itr = crossSections.begin();
  std::vector<TH1*>::const_iterator eventsNumberTH1vector_itr  = eventsNumberTH1vector.begin();  
  int index = 0;
  for ( ; eventsNumberTH1vector_itr != eventsNumberTH1vector.end(); ++eventsNumberTH1vector_itr,
	                                                            ++crossSections_itr,
	                                                            index++) { 
    TH1 * histo = (TH1*)(*eventsNumberTH1vector_itr)->Clone();
    crossSections[index] = *crossSections_itr/double(histo->GetMaximum());
    delete histo;
  }
  
  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );
  
  std::vector<TString>::const_iterator varNameVector_itr = varNameVector.begin();
  for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
    findObj( TH1vector, 
	     FileList, 
	     *varNameVector_itr, 
	     &FileSuffix,
	     "TH1",
	     &crossSections
	     );
    
    TH1 * mergedHisto = mergeObj("merged_"+*varNameVector_itr,
				 TH1vector,
				 &crossSections);
    // THStackLegend histogram
    TString Stack( "Stack_" );
    THStackLegend<TH1> * StackLegend_ = new THStackLegend<TH1>( *varNameVector_itr );
    
    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameVector_itr, *varNameVector_itr );
    directory->cd();
    StackLegend_->Add(mergedHisto,"merged",true);

    mergedHisto->Write();
    std::vector<TH1*>::const_iterator    TH1vector_itr  = TH1vector.begin();
    std::vector<TString>::const_iterator FileSuffix_itr = FileSuffix.begin();
    for ( ; TH1vector_itr != TH1vector.end(); ++TH1vector_itr,
	                                      ++FileSuffix_itr) {
      TH1 * histo = (TH1*)(*TH1vector_itr)->Clone();
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
      StackLegend_->Add(dressedHisto,*FileSuffix_itr,true,"lf",false,dressed);

      dressedHisto->Write();
    }
    StackLegend_->Write("nostack");
    StackLegend_->Draw("nostack");
    
    TH1vector.clear();
  }
  
}


