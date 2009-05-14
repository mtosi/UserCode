/**

  This macro will add histograms multiplied by their cross-section
  from a list of root files and write them
  to a target root file. The target file is newly created and must not be
  identical to one of the source files.


  Author: Sven A. Schmidt, sven.schmidt@cern.ch
  Date:   13.2.2001

  Editing Author: Michael B. Anderson, mbanderson@hep.wisc.edu
  Date:  July 12, 2007

  This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
  which had a problem with directories more than one level deep.
  (see macro hadd_old.C for this previous implementation).

  The macro from Sven has been enhanced by
     Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
   to automatically add Trees (via a chain of trees).

  To use this macro, modify the file names in function hadd.

  NB: This macro is provided as a tutorial.
      Use $ROOTSYS/bin/hadd to merge many histogram files

 */


#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

template<class T1>
T1* mergeObj( TString               outputObjname,
	      std::vector<T1*>    & inputVector
	      ) {
  T1 * outputObj = (T1*)(inputVector[0])->Clone();
  outputObj->SetName(outputObjname);
  outputObj->Sumw2();
  typename std::vector<T1*>::const_iterator inputVector_itr = inputVector.begin()+1;
  for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr ) {
    T1 * histo = (T1*)(*inputVector_itr)->Clone();
    outputObj->Add(histo);
    delete histo;
  }
  return outputObj;
}

template<class T1>
T1* mergeObj( TString               outputObjname,
	      std::vector<T1*>    & inputVector,
	      std::vector<double> * xSecVector
	      ) {
  T1 * outputObj = (T1*)(inputVector[0])->Clone();
  outputObj->Sumw2();
  if (xSecVector) outputObj->Scale((*xSecVector)[0]);
  outputObj->SetName(outputObjname);
  typename std::vector<T1*>::const_iterator inputVector_itr = inputVector.begin()+1;
  int index = 1;
  for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr, index++ ) {
    T1 * histo = (T1*)(*inputVector_itr)->Clone();
    histo->Scale((*xSecVector)[index]);
    outputObj->Add(histo);
    delete histo;
  }
  return outputObj;
}

template<class T1>
T1* mergeObj( TString               outputObjname,
	      std::vector<T1*>    & inputVector,
	      std::vector<double> * xSecVector,
	      std::vector<double> * totNumberVector,
	      bool                  rebinning = kTRUE
	      ) {
  T1 * outputObj;
  //  if (rebinning) outputObj = ((T1*)(inputVector[0])->Clone())->Rebin();
  if (rebinning) outputObj = ((T1*)(inputVector[0])->Clone())->Rebin(5);
  else outputObj = (T1*)(inputVector[0])->Clone();
  outputObj->SetName(outputObjname);
  int nbins = outputObj->GetNbinsX();
  for ( int ibin = 1; ibin <= nbins; ibin++ ) {
    typename std::vector<T1*>::const_iterator inputVector_itr = inputVector.begin();
    int index = 0;
    double output_ibin_content = 0.;
    double output_ibin_error2  = 0.;
    for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr, index++ ) {
      //      double ibin_content = (T1*)((*inputVector_itr)->Rebin())->GetBinContent(ibin);
      //      double ibin_error   = (T1*)((*inputVector_itr)->Rebin())->GetBinError(ibin);
      double ibin_content = (T1*)((*inputVector_itr)->Rebin(5))->GetBinContent(ibin);
      double ibin_error   = (T1*)((*inputVector_itr)->Rebin(5))->GetBinError(ibin);
      output_ibin_content += ibin_content/(*totNumberVector)[index]*(*xSecVector)[index];
      output_ibin_error2  += pow(ibin_error,2)*ibin_content*
	((*totNumberVector)[index]-ibin_content)/pow((*totNumberVector)[index],2);
    }
    outputObj->SetBinContent(ibin,output_ibin_content);
    outputObj->SetBinError(ibin,TMath::Sqrt(output_ibin_error2));
  }
  return outputObj;
}

template<class T1>
T1* mergeObj( TString               outputObjname,
	      std::vector<T1*>    & inputVector,
	      std::vector<double> * xSecVector,
	      bool                  errorEnforce = kTRUE,
	      bool                  rebinning    = kTRUE
	      ) {
  T1 * outputObj;
  //  if (rebinning) outputObj = ((T1*)(inputVector[0])->Clone())->Rebin();
  if (rebinning) outputObj = ((T1*)(inputVector[0])->Clone())->Rebin(5);
  else outputObj = (T1*)(inputVector[0])->Clone();
  outputObj->SetName(outputObjname);
  if (!errorEnforce) {
    outputObj->Sumw2();
    if (xSecVector) outputObj->Scale((*xSecVector)[0]);
    typename std::vector<T1*>::const_iterator inputVector_itr = inputVector.begin()+1;
    int index = 1;
    for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr, index++ ) {
      T1 * histo = (T1*)(*inputVector_itr)->Clone();
      histo->Scale((*xSecVector)[index]);
      outputObj->Add(histo);
      delete histo;
    }
  } else {
    int nbins = outputObj->GetNbinsX();
    for ( int ibin = 1; ibin <= nbins; ibin++ ) {
      typename std::vector<T1*>::const_iterator inputVector_itr = inputVector.begin();
      int index = 0;
      double output_ibin_content = 0.;
      double output_ibin_error2  = 0.;
      for ( ; inputVector_itr != inputVector.end(); ++inputVector_itr, index++ ) {
	T1 * obj;
	//if (rebinning) obj = ((T1*)(*inputVector_itr)->Clone())->Rebin();
	if (rebinning) obj = ((T1*)(*inputVector_itr)->Clone())->Rebin(5);
	else obj = (T1*)(*inputVector_itr)->Clone();
	double ibin_content = obj->GetBinContent(ibin);
	double ibin_error   = obj->GetBinError(ibin);
	output_ibin_content += ibin_content*(*xSecVector)[index];
	output_ibin_error2  += pow(ibin_error*(*xSecVector)[index],2);
      }
      outputObj->SetBinContent(ibin,output_ibin_content);
      outputObj->SetBinError(ibin,TMath::Sqrt(output_ibin_error2));
    }
  }
  return outputObj;
}
