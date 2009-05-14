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
#include "TH2.h"
#include "Riostream.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/findObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/THStackLegend.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ListFashionAttributedHisto.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BR.h"

using namespace vbfhzz2l2b;

void plotHisto2d(TString outputFileName = "plot2dOutput.root") {

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
  for ( ; FileSuffix_itr != FileSuffix.end(); ++FileSuffix_itr ) {
    FileList->Add( TFile::Open(*FileSuffix_itr+"output.root"  ) );
  }
  // ************************************************************
  // List of xSec (in pb)
  // stored in 
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecNLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BR.h"
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
  // ************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("pTbalanceVSMEt");
  // ************************************************************

  std::vector<TH2*> TH2vector;
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
  }

  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );

  std::vector<TString>::const_iterator varNameVector_itr = varNameVector.begin();
  for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
    findObj( TH2vector, 
	     FileList, 
	     *varNameVector_itr, 
	     &FileSuffix,
	     "TH2",
	     &crossSections
	     );

    // THStackLegend histogram
    TString Stack( "Stack_" );
    THStackLegend<TH2> * StackLegend_ = new THStackLegend<TH2>( *varNameVector_itr, 0.6,0.68,0.98,0.97 );
    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameVector_itr, *varNameVector_itr );
    directory->cd();
    std::vector<TH2*>::const_iterator    TH2vector_itr  = TH2vector.begin();
    std::vector<TString>::const_iterator FileSuffix_itr = FileSuffix.begin();
    crossSections_itr = crossSections.begin();
    for ( ; TH2vector_itr != TH2vector.end(); ++TH2vector_itr, 
	                                      ++FileSuffix_itr,
	                                      ++crossSections_itr) {
      TH2 * histo = (TH2*)(*TH2vector_itr)->Clone();
      FashionAttributedHisto<TH2D> * dressedHisto;
      bool dressed = false;
      if (      *FileSuffix_itr == "VBFH160ZZ" ) { dressedHisto = new H160<TH2D>((TH2D*)histo);     dressed = true; }
      else if ( *FileSuffix_itr == "VBFH200ZZ" ) { dressedHisto = new H200<TH2D>((TH2D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "VBFH400ZZ" ) { dressedHisto = new H400<TH2D>((TH2D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "VBFH800ZZ" ) { dressedHisto = new H800<TH2D>((TH2D*)histo);	    dressed = true; }
      else if ( *FileSuffix_itr == "ZZ0Jets"   ) { dressedHisto = new ZZ0jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "ZZ1Jets"   ) { dressedHisto = new ZZ1jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "ZZ2Jets"   ) { dressedHisto = new ZZ2jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "Zbb0Jets"  ) { dressedHisto = new Zbb0jets<TH2D>((TH2D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "Zbb1Jets"  ) { dressedHisto = new Zbb1jets<TH2D>((TH2D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "Zbb2Jets"  ) { dressedHisto = new Zbb2jets<TH2D>((TH2D*)histo); dressed = true; }
      else if ( *FileSuffix_itr == "WZ0Jets"   ) { dressedHisto = new WZ0jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "WZ1Jets"   ) { dressedHisto = new WZ1jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "WZ2Jets"   ) { dressedHisto = new WZ2jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt0Jets"   ) { dressedHisto = new tt0jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt1Jets"   ) { dressedHisto = new tt1jets<TH2D>((TH2D*)histo);  dressed = true; }
      else if ( *FileSuffix_itr == "tt2Jets"   ) { dressedHisto = new tt2jets<TH2D>((TH2D*)histo);  dressed = true; }
      else dressedHisto = new FashionAttributedHisto<TH2D>((TH2D*)histo);
      dressedHisto->SetMarkerSize(0.4);
      StackLegend_->Add(dressedHisto,*FileSuffix_itr,false,"p",false,dressed);
      TCanvas * canvas = new TCanvas();
      canvas->cd();
      dressedHisto->Draw("lego2");
      dressedHisto->Write();
    }
    StackLegend_->Print("nostack",*varNameVector_itr+"MultiPlot.jpg");
    StackLegend_->Write("nostack");
    StackLegend_->Draw("nostack");
    //    StackLegend_->SavePrimitive(Stack+*varNameVector_itr+".C","nostack");

    TH2vector.clear();

  }

}


