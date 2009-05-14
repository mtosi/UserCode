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
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecNLO.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BR.h"

using namespace vbfhzz2l2b;

void plotDummy(TString outputFileName = "plotDummy.root") {

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
  std::vector<double>  sgnSampleEvents;
  sgnFileSuffix.push_back("TMVA200");sgnSampleEvents.push_back(90000.);
  //  sgnFileSuffix.push_back("TMVA800");sgnSampleEvents.push_back(96000.);
  std::vector<TString> bkgZZFileSuffix;
  std::vector<double> bkgZZSampleEvents;
  bkgZZFileSuffix.push_back("TMVA_ZZ0JETS");bkgZZSampleEvents.push_back(10822.);
  bkgZZFileSuffix.push_back("TMVA_ZZ1JETS");bkgZZSampleEvents.push_back(5546.);
  bkgZZFileSuffix.push_back("TMVA_ZZ2JETS");bkgZZSampleEvents.push_back(9572.);
  // ************************************************************
  // List of Files
  TList * sgnFileList    = new TList();
  std::vector<TString>::const_iterator sgnFileSuffix_itr = sgnFileSuffix.begin();
  for ( ; sgnFileSuffix_itr != sgnFileSuffix.end(); ++sgnFileSuffix_itr )
    sgnFileList->Add( TFile::Open(*sgnFileSuffix_itr+".root"  ) );
  TList * bkgZZFileList  = new TList();
  std::vector<TString>::const_iterator bkgZZFileSuffix_itr = bkgZZFileSuffix.begin();
  for ( ; bkgZZFileSuffix_itr != bkgZZFileSuffix.end(); ++bkgZZFileSuffix_itr )
    bkgZZFileList->Add( TFile::Open(*bkgZZFileSuffix_itr+".root"  ) );
    // ************************************************************
  // List of xSec (in pb)
  // stored in 
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/xSecNLO.h"
  // "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BR.h"
  std::vector<double> sgnCrossSections;
  std::vector<double> bkgZZCrossSections;
  sgnCrossSections.push_back(1000*xSec_VBFH200_*BR_H200ZZ_);
  //sgnCrossSections.push_back(1000*xSec_VBFH800_*BR_H800ZZ_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ0Jets_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ1Jets_);
  bkgZZCrossSections.push_back(1000*xSec_ZZ2Jets_);

  // *************************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("recHlljj_final"                                  );
  // *************************************************************************

  std::vector<TH1*> sgnEventsNumberTH1vector;
  if (sgnFileSuffix.size() != 0) {
    //    findObj( sgnEventsNumberTH1vector, 
    //	     sgnFileList, 
    //	     "eventsNumber", 
    //	     &sgnFileSuffix
    //	     );
    for ( int indexSample = 0; indexSample < sgnFileSuffix.size();indexSample++ ) { 
      if(sgnSampleEvents[indexSample]!=0.) sgnCrossSections[indexSample] = sgnCrossSections[indexSample]/sgnSampleEvents[indexSample];
    }
    for ( int indexSample = 0; indexSample < bkgZZFileSuffix.size();indexSample++ ) { 
      if(bkgZZSampleEvents[indexSample]!=0.) bkgZZCrossSections[indexSample] = bkgZZCrossSections[indexSample]/bkgZZSampleEvents[indexSample];
    }
  }

  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );

  std::vector<TString>::const_iterator varNameVector_itr = varNameVector.begin();
  for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
    std::vector<TH1*> TH1vector;
    std::vector<TH1*> bkgZZTH1vector;
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
    TH1vector.push_back(mergeObj((*varNameVector_itr)+"_ZZNjets", bkgZZTH1vector));
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
      if ( sample.Contains("200") ) { dressedHisto = new H200<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( sample.Contains("800") ) { dressedHisto = new H800<TH1D>((TH1D*)histo);	    dressed = true; }
      else if ( sample.Contains("ZZ" ) && sample.Contains("jets") ) { dressedHisto = new ZZNjets<TH1D>((TH1D*)histo);  dressed = true; }
      else dressedHisto = new FashionAttributedHisto<TH1D>((TH1D*)histo);
      char nev[50];
      sprintf(nev,sample+": %.2f",dressedHisto->Integral());
      StackLegend_->Add(dressedHisto,nev,false,"lf",false,dressed);
      //      StackLegend_->Add(dressedHisto,sample,false,"fl",false,dressed);

      dressedHisto->Write();
    }
    StackLegend_->SetLogY();
    StackLegend_->Print("nostack",*varNameVector_itr+"Dummy200.jpg");
    StackLegend_->Write("nostackhisto");
    StackLegend_->Draw("nostackhisto");
    //    StackLegend_->SavePrimitive(Stack+*varNameVector_itr+".C","nostackhisto");
    
    //    TH1vector.clear();
    //    bkgZZTH1vector.clear();
    //    bkgZbbTH1vector.clear();
  }
  
  }
  

