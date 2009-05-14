#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "Riostream.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/mergeObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/findObj.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/THStackLegend.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ListFashionAttributedHisto.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BR.h"

using namespace vbfhzz2l2b;

void plotEfficiencyCuts(TString outputFileName = "plotEfficiencyOutput.root") {

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
  // ************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("recHlljj"                                );
  varNameVector.push_back("recHlljj_leptonicZcutsAFTERZdileptonMass");
  varNameVector.push_back("recHlljj_leptonicZcutsAFTERZleptonDeltaR");
  varNameVector.push_back("recHlljj_hadronicZcutsAFTERZdijetMass"   );
  varNameVector.push_back("recHlljj_hadronicZcutsAFTERZjetDeltaR"   );
  varNameVector.push_back("recHlljj_FWDcutsAFTERleptonsShiftedEta"  );
  varNameVector.push_back("recHlljj_FWDcutsAFTERjetsShiftedEta"     );
  varNameVector.push_back("recHlljj_FWDcutsAFTERforwardDiJetMass"   );
  varNameVector.push_back("recHlljj_FWDcutsAFTERpTbalance"          );
  varNameVector.push_back("recHlljj_final"                          );
  // ************************************************************
  std::vector<TString> varNameCutVector;
  varNameCutVector.push_back("ZdileptonMass"    );
  varNameCutVector.push_back("ZleptonDeltaR"    );
  varNameCutVector.push_back("ZdijetMass"       );
  varNameCutVector.push_back("ZjetDeltaR"       );
  varNameCutVector.push_back("leptonsShiftedEta");
  varNameCutVector.push_back("jetsShiftedEta"   );
  varNameCutVector.push_back("forwardDiJetMass" );
  varNameCutVector.push_back("pTbalance"        );
  varNameCutVector.push_back("final"            );
  // ************************************************************

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


  // output latex file
  ofstream latexfile;
  latexfile.open("efficiency.tex", ofstream::out | ofstream::app);

  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );
  int indexCut = 0;
  std::vector<TString>::const_iterator varNameVector_itr;
  std::vector<TString>::const_iterator varNameCutVector_itr = varNameCutVector.begin();
  for ( ; varNameCutVector_itr != varNameCutVector.end(); ++varNameCutVector_itr,
	                                                  indexCut++ ) {
    std::vector<TH1*> afterSgnTH1vector;
    std::vector<TH1*> beforeSgnTH1vector;
    std::vector<TH1*> afterBkgTH1vector;
    std::vector<TH1*> beforeBkgTH1vector;
    std::vector<TH1*> afterbkgZZTH1vector;
    std::vector<TH1*> beforebkgZZTH1vector;
    std::vector<TH1*> afterbkgZbbTH1vector;
    std::vector<TH1*> beforebkgZbbTH1vector;
    varNameVector_itr = varNameVector.begin();
    for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
      TString varName = *varNameVector_itr;
      if (varName.Contains(*varNameCutVector_itr)) {
	std::cout << "varName: " << varName << std::endl;
	findObj( afterbkgZZTH1vector, 
		 bkgZZFileList, 
		 *varNameVector_itr, 
		 &bkgZZFileSuffix
		 );
	findObj( beforebkgZZTH1vector, 
		 bkgZZFileList, 
		 *(varNameVector_itr-1), 
		 &bkgZZFileSuffix
		 );
	findObj( afterbkgZbbTH1vector, 
		 bkgZbbFileList, 
		 *varNameVector_itr, 
		 &bkgZbbFileSuffix
		 );
	findObj( beforebkgZbbTH1vector, 
		 bkgZbbFileList, 
		 *(varNameVector_itr-1), 
		 &bkgZbbFileSuffix
		 );

	findObj( afterSgnTH1vector,
		 sgnFileList, 
		 *varNameVector_itr, 
		 &sgnFileSuffix
		 );
	findObj( beforeSgnTH1vector,
		 sgnFileList, 
		 *(varNameVector_itr-1), 
		 &sgnFileSuffix
		 );

	afterBkgTH1vector.push_back( mergeObj((*varNameVector_itr)    +"_ZZNjets", afterbkgZZTH1vector,  &bkgZZCrossSections, "kTRUE","kTRUE"));
	beforeBkgTH1vector.push_back(mergeObj((*(varNameVector_itr-1))+"_ZZNjets", beforebkgZZTH1vector, &bkgZZCrossSections, "kTRUE","kTRUE"));
	afterBkgTH1vector.push_back( mergeObj((*varNameVector_itr)    +"_ZbbNjets",afterbkgZbbTH1vector, &bkgZbbCrossSections,"kTRUE","kTRUE"));
	beforeBkgTH1vector.push_back(mergeObj((*(varNameVector_itr-1))+"_ZbbNjets",beforebkgZbbTH1vector,&bkgZbbCrossSections,"kTRUE","kTRUE"));
      }
    }

    int    nbins = 100;
    double xmin  = 0.;
    double xmax  = 1000.;
    // sgn efficiency histogram
    TH1D * efficiencySgnTH1 = new TH1D(*varNameCutVector_itr+"efficiencyHZZ",*varNameCutVector_itr+" efficiency",nbins,xmin,xmax);
    std::vector<TH1*>::const_iterator beforeSgnTH1vector_itr = beforeSgnTH1vector.begin();
    std::vector<TH1*>::const_iterator afterSgnTH1vector_itr = afterSgnTH1vector.begin();
    for ( ; beforeSgnTH1vector_itr != beforeSgnTH1vector.end(); ++beforeSgnTH1vector_itr,
	                                                        ++afterSgnTH1vector_itr ) {
      TString histoName = (*afterSgnTH1vector_itr)->GetName();
      TString sample = (histoName.Remove(0,histoName.Index((*varNameCutVector_itr))+(*varNameCutVector_itr).Length()+1));
      latexfile << " " << (*varNameCutVector_itr) << endl;
      latexfile << " " << sample << endl;
      latexfile << "\\hline" << endl;
      double eventsBefore = (*beforeSgnTH1vector_itr)->Integral();
      double eventsAfter  = (*afterSgnTH1vector_itr)->Integral();
      double efficiency       = 0.;
      double efficiencyError2 = 0.;
      if( eventsBefore != 0. ) {
	efficiency = eventsAfter/eventsBefore;
	efficiencyError2 = efficiency*(1-efficiency)/eventsBefore;
      }
      int ibin = 0;
      if (      sample.Contains("VBFH160") ) ibin = 16;
      else if ( sample.Contains("VBFH200") ) ibin = 20;
      else if ( sample.Contains("VBFH400") ) ibin = 40;
      else if ( sample.Contains("VBFH800") ) ibin = 80;
      efficiencySgnTH1->SetBinContent(ibin,efficiency);
      efficiencySgnTH1->SetBinError(ibin,TMath::Sqrt(efficiencyError2));
      latexfile << " " << efficiency << " " << TMath::Sqrt(efficiencyError2) << endl;
      std::cout << sample << ": " << efficiency << "+-" << TMath::Sqrt(efficiencyError2) << std::endl;
    }
    // bkg efficiency histo vector
    vector<TH1D*> efficiencyBkgTH1vector;
    std::vector<TH1*>::const_iterator beforeBkgTH1vector_itr = beforeBkgTH1vector.begin();
    std::vector<TH1*>::const_iterator afterBkgTH1vector_itr = afterBkgTH1vector.begin();
    for ( ; beforeBkgTH1vector_itr != beforeBkgTH1vector.end(); ++beforeBkgTH1vector_itr,
	                                                        ++afterBkgTH1vector_itr ) {
      TString histoName = (*afterBkgTH1vector_itr)->GetName();
      TString sample = (histoName.Remove(0,histoName.Index((*varNameCutVector_itr))+(*varNameCutVector_itr).Length()+1));
      TH1D * efficiencyBkgTH1 = new TH1D(*varNameCutVector_itr+"_efficiency_"+sample,*varNameCutVector_itr+" efficiency "+sample,1,xmin,xmax);
      double eventsBefore = (*beforeBkgTH1vector_itr)->Integral();
      double eventsAfter  = (*afterBkgTH1vector_itr)->Integral();
      double efficiency       = 0.;
      double efficiencyError2 = 0.;
      if( eventsBefore != 0. ) {
	efficiency = eventsAfter/eventsBefore;
	efficiencyError2 = efficiency*(1-efficiency)/eventsBefore;
      }
      efficiencyBkgTH1->SetBinContent(1,efficiency);
      efficiencyBkgTH1->SetBinError(1,TMath::Sqrt(efficiencyError2));
      efficiencyBkgTH1vector.push_back((TH1D*)efficiencyBkgTH1->Clone());
      delete efficiencyBkgTH1;
      latexfile << " " << efficiency << " " << TMath::Sqrt(efficiencyError2) << endl;
      std::cout << sample << ": " << efficiency << "+-" << TMath::Sqrt(efficiencyError2) << std::endl;
    }
    latexfile << "\\hline" << endl;

    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameCutVector_itr, *varNameCutVector_itr );
    directory->cd();
    TCanvas * canvas = new TCanvas(*varNameCutVector_itr,*varNameCutVector_itr,1000,800);
    canvas->cd();
    TLegend * legend = new TLegend(0.65,0.3,0.85,0.6);

    FashionAttributedHisto<TH1D> * dressedSgnHisto = new H400<TH1D>((TH1D*)(efficiencySgnTH1->Clone()));
    dressedSgnHisto->SetMaximum(1.);
    dressedSgnHisto->GetXaxis()->SetTitle("m_{H} (GeV/c^{2})");
    dressedSgnHisto->GetXaxis()->SetTitleSize(0.04);
    dressedSgnHisto->GetYaxis()->SetTitle("efficiency (%)");
    dressedSgnHisto->GetYaxis()->SetTitleSize(0.04);
    dressedSgnHisto->GetXaxis()->SetLabelSize(0.05);
    dressedSgnHisto->GetXaxis()->SetBinLabel(16,"160");
    dressedSgnHisto->GetXaxis()->SetBinLabel(20,"200");
    dressedSgnHisto->GetXaxis()->SetBinLabel(40,"400");
    dressedSgnHisto->GetXaxis()->SetBinLabel(80,"800");
    legend->AddEntry(dressedSgnHisto,"HZZ","lpe");
    dressedSgnHisto->Draw("pe");
    
    vector<TH1D*>::const_iterator efficiencyBkgTH1vector_itr = efficiencyBkgTH1vector.begin();
    for ( ; efficiencyBkgTH1vector_itr != efficiencyBkgTH1vector.end(); ++efficiencyBkgTH1vector_itr ) {
      TString histoName = (*efficiencyBkgTH1vector_itr)->GetName();
      TString sample = histoName.Remove(0,histoName.Index(TString("efficiency"))+(TString("efficiency")).Length()+1);
      FashionAttributedHisto<TH1D> * dressedBkgHisto;
      if      ( sample.Contains("ZZ" ) && sample.Contains("jets") ) dressedBkgHisto = new ZZNjets<TH1D>( (TH1D*)(*efficiencyBkgTH1vector_itr)->Clone());
      else if ( sample.Contains("Zbb") && sample.Contains("jets") ) dressedBkgHisto = new ZbbNjets<TH1D>((TH1D*)(*efficiencyBkgTH1vector_itr)->Clone());
      else if ( sample.Contains("WZ" ) && sample.Contains("jets") ) dressedBkgHisto = new WZNjets<TH1D>( (TH1D*)(*efficiencyBkgTH1vector_itr)->Clone());
      else if ( sample.Contains("tt" ) && sample.Contains("jets") ) dressedBkgHisto = new ttNjets<TH1D>( (TH1D*)(*efficiencyBkgTH1vector_itr)->Clone());
      else dressedBkgHisto = new FashionAttributedHisto<TH1D>((TH1D*)(*efficiencyBkgTH1vector_itr)->Clone()); 
      dressedBkgHisto->SetMaximum(1.);
      dressedBkgHisto->Draw("e2lsame");
      legend->AddEntry(dressedBkgHisto,sample,"lf");
      dressedBkgHisto->Write();
    }
    legend->Draw();

    efficiencySgnTH1->Write();

    canvas->Write();
    canvas->Print(*varNameCutVector_itr+"EfficiencyPlot.jpg");

  }
  latexfile.close();

}

