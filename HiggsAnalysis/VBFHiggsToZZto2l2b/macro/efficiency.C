#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "Riostream.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/findObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/THStackLegend.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ListFashionAttributedHisto.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"

using namespace vbfhzz2l2b;

void plotEfficiencyCuts(TString outputFileName = "plotOutput.root") {

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
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"
  std::vector<double> crossSections;
  crossSections.push_back(1000*xSec_VBFH160_*BR_H160ZZ_);
  crossSections.push_back(1000*xSec_VBFH200_*BR_H200ZZ_);
  crossSections.push_back(1000*xSec_VBFH400_*BR_H400ZZ_);
  crossSections.push_back(1000*xSec_VBFH800_*BR_H800ZZ_);
  crossSections.push_back(1000*xSec_ZZ0Jets_*BR_Zll_);
  crossSections.push_back(1000*xSec_ZZ1Jets_*BR_Zll_);
  crossSections.push_back(1000*xSec_ZZ2Jets_*BR_Zll_);
  crossSections.push_back(1000*xSec_Zbb0Jets_*BR_Zll_);
  crossSections.push_back(1000*xSec_Zbb1Jets_*BR_Zll_);
  crossSections.push_back(1000*xSec_Zbb2Jets_*BR_Zll_);
  // ************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("ZbeforeDeltaRcutNumber");
  varNameVector.push_back("ZafterDeltaRcutNumber");
  varNameVector.push_back("ZbeforeLeptonShiftedEtaCutNumber");
  varNameVector.push_back("ZafterLeptonShiftedEtaCutNumber");
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

    TH1D * efficiencyHZZ = new TH1D("efficiencyHZZ","efficiency",100,0.,1000.);
    efficiencyHZZ->SetBinContent(16,0.40);
    efficiencyHZZ->SetBinContent(20,0.44);
    efficiencyHZZ->SetBinContent(40,0.60);
    efficiencyHZZ->SetBinContent(80,0.80);
    efficiencyHZZ->SetBinError(16,0.03);
    efficiencyHZZ->SetBinError(20,0.03);
    efficiencyHZZ->SetBinError(40,0.04);
    efficiencyHZZ->SetBinError(80,0.04);

    efficiencyHZZ->GetXaxis()->SetTitle("m_{H} (GeV/c^{2})");
    efficiencyHZZ->GetXaxis()->SetLabelSize(0.05);
    efficiencyHZZ->GetXaxis()->SetBinLabel(16,"160");
    efficiencyHZZ->GetXaxis()->SetBinLabel(20,"200");
    efficiencyHZZ->GetXaxis()->SetBinLabel(40,"400");
    efficiencyHZZ->GetXaxis()->SetBinLabel(80,"800");

    TLine * efficiencyZbb = new TLine(0.,0.7,1000.,0.7);
    TLine * efficiencyZZ  = new TLine(0.,0.65,1000.,0.65);

    TLegend * legend = new TLegend(0.65,0.3,0.85,0.6);
    legend->AddEntry(efficiencyHZZ,"HZZ","lpe");
    legend->AddEntry(efficiencyZbb,"Zbb+n jets","l");
    legend->AddEntry(efficiencyZZ, "ZZ+n jets","l");

    efficiencyHZZ->Draw();
    efficiencyZbb->Draw("same");
    efficiencyZZ->Draw("same");
    legend->Draw();
//    // HZZ efficiency histogram
//    THStackLegend<TH1> * StackLegend_ = new THStackLegend<TH1>( *varNameVector_itr );
//
//    outputFile->cd();
//    TDirectory * directory = outputFile->mkdir( *varNameVector_itr, *varNameVector_itr );
//    directory->cd();
//    std::vector<TH1*>::const_iterator    TH1vector_itr  = TH1vector.begin();
//    std::vector<TString>::const_iterator FileSuffix_itr = FileSuffix.begin();
//    for ( ; TH1vector_itr != TH1vector.end(); ++TH1vector_itr, 
//	                                      ++FileSuffix_itr) {
//      TH1 * histo = (TH1*)(*TH1vector_itr)->Clone();
//      FashionAttributedHisto<TH1D> * dressedHisto;
//      bool dressed = false;
//      if (      *FileSuffix_itr == "VBFH160ZZ" ) { dressedHisto = new H160<TH1D>((TH1D*)histo);     dressed = true; }
//      else if ( *FileSuffix_itr == "VBFH200ZZ" ) { dressedHisto = new H200<TH1D>((TH1D*)histo);	    dressed = true; }
//      else if ( *FileSuffix_itr == "VBFH400ZZ" ) { dressedHisto = new H400<TH1D>((TH1D*)histo);	    dressed = true; }
//      else if ( *FileSuffix_itr == "VBFH800ZZ" ) { dressedHisto = new H800<TH1D>((TH1D*)histo);	    dressed = true; }
//      else if ( *FileSuffix_itr == "ZZ0Jets"   ) { dressedHisto = new ZZ0jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "ZZ1Jets"   ) { dressedHisto = new ZZ1jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "ZZ2Jets"   ) { dressedHisto = new ZZ2jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "Zbb0Jets"  ) { dressedHisto = new Zbb0jets<TH1D>((TH1D*)histo); dressed = true; }
//      else if ( *FileSuffix_itr == "Zbb1Jets"  ) { dressedHisto = new Zbb1jets<TH1D>((TH1D*)histo); dressed = true; }
//      else if ( *FileSuffix_itr == "Zbb2Jets"  ) { dressedHisto = new Zbb2jets<TH1D>((TH1D*)histo); dressed = true; }
//      else if ( *FileSuffix_itr == "WZ0Jets"   ) { dressedHisto = new WZ0jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "WZ1Jets"   ) { dressedHisto = new WZ1jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "WZ2Jets"   ) { dressedHisto = new WZ2jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "tt0Jets"   ) { dressedHisto = new tt0jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "tt1Jets"   ) { dressedHisto = new tt1jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else if ( *FileSuffix_itr == "tt2Jets"   ) { dressedHisto = new tt2jets<TH1D>((TH1D*)histo);  dressed = true; }
//      else dressedHisto = new FashionAttributedHisto<TH1D>((TH1D*)histo);
//      StackLegend_->Add(dressedHisto,*FileSuffix_itr,true,"l",false,dressed);
//
//      dressedHisto->Write();
//    }
//    StackLegend_->Write("nostack");
//    StackLegend_->Draw("nostack");
//    //    StackLegend_->SavePrimitive(Stack+*varNameVector_itr+".C","nostack");
//
//    TH1vector.clear();
//
  }
}

  

