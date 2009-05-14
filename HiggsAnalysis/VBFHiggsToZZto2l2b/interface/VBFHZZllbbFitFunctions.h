#ifndef VBFHZZllbbFITFUNCTIONS_H
#define VBFHZZllbbFITFUNCTIONS_H

#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TMath.h>
#include <TChain.h>
#include <TProfile2D.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include <TLatex.h>
 
// system include files
#include <iostream>
#include <fstream>
#include <string> 
#include <memory>

//###############//
// Fit functions //
//###############//

namespace vbfhzz2l2b
{

  // Erf function
  double fit_erf(double *x, double *par){
    return par[0]+par[1]*TMath::Erf((x[0]-par[2])/par[3]);
  }


  // power 7 function
  double fit_pol7(double *x, double *par){
    const int Npow = 7;
    double fpol = 0;
    for ( int i=0;i<=Npow;i++ ){
      fpol += par[i] * pow(x[0],i);
    }
    return fpol;
  }

  // power 7 function
  double fit_pol7_rejectRange(double *x, double *par, double & xMin, double & xMax){
    if (x[0] > xMin && x[0] < xMax) {
      TF1::RejectPoint();
      return 0;
    }
    return fit_pol7(x,par);
  }

  // UA2 function
  double fit_UA2(double *x, double *par){
    return  par[0]*pow(x[0],-par[1])*exp(-par[2]*x[0])*exp(-par[3]*x[0]*x[0]);
  }


  double fit_UA2_rejectRange(double *x, double *par, double & xMin, double & xMax){
    if (x[0] > xMin && x[0] < xMax) {
      TF1::RejectPoint();
      return 0;
    }
    return  fit_UA2(x,par);
  }
  
  
double fit_UA2_tot(double *x, double *par){
  // Signal and BG function
  const float pi = 3.1415926;
  double argom =0;
  if (par[2] !=0) argom = (x[0] - par[5])/par[6];
  double fitval =  par[0]*pow(x[0],-par[1])*exp(-par[2]*x[0])*exp(-par[3]*x[0]*x[0]) + par[4]*exp(-0.5*argom*argom)/(sqrt(2*pi)*par[6]);
  return fitval;
}


double fit_gaus(double *x, double *par){  
  double arg =0;
  if ( par[2]!=0 ) arg = (x[0] - par[1])/par[2];      
  double fitval =par[0]*exp(-0.5*arg*arg)/(sqrt(2*TMath::Pi())*par[2]);
  return fitval;
}


double fit2gaus(double *x, double *par){ // fit w/ two gausian
  const float pi = 3.1415926;
  double arg1 = 0;
  double arg2 = 0;
  if ( par[2]!=0 ) arg1 = (x[0] - par[1])/par[2];      
  if ( par[4]!=0 ) arg2 = (x[0] - par[3])/par[4];      
  double fitval = par[0] * exp(-0.5*arg1*arg1)/(sqrt(2*pi)*par[2]) * exp(-0.5*arg2*arg2)/(sqrt(2*pi)*par[4]);
  return fitval;
}

double fit2gausSum(double *x, double *par){ // fit w/ sum of two gausian
  const float pi = 3.1415926;
  double arg1 = 0;
  double arg2 = 0;
  if ( par[2]!=0 ) arg1 = (x[0] - par[1])/par[2];
  if ( par[5]!=0 ) arg2 = (x[0] - par[4])/par[5];      
  double fitval1 = par[0] * exp(-0.5*arg1*arg1)/(sqrt(2*pi)*par[2]);
  double fitval2 = par[3] * exp(-0.5*arg2*arg2)/(sqrt(2*pi)*par[5]);
  return fitval1 + fitval2;
}

double fit3gausSum(double *x, double *par){ // fit w/ sum of three gausian
  const float pi = 3.1415926;
  double arg1 = 0;
  double arg2 = 0;
  double arg3 = 0;
  if ( par[2]!=0 ) arg1 = (x[0] - par[1])/par[2];
  if ( par[5]!=0 ) arg2 = (x[0] - par[4])/par[5];
  if ( par[8]!=0 ) arg3 = (x[0] - par[7])/par[8];
  double fitval1 = par[0] * exp(-0.5*arg1*arg1)/(sqrt(2*pi)*par[2]);
  double fitval2 = par[3] * exp(-0.5*arg2*arg2)/(sqrt(2*pi)*par[5]);
  double fitval3 = par[6] * exp(-0.5*arg3*arg3)/(sqrt(2*pi)*par[8]);
  return par[9] * (fitval1 + fitval2 + fitval3);
}

double fitlandau(double *x, double *par){
  return par[0]*TMath::Landau(x[0],par[1],par[2]);
}


double fitlandauUA2(double *x, double *par){
  return par[0]*TMath::Landau(x[0],par[1],par[2]) + par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]);
  //  if ( x[0]>50 && x[0]<75 ) return par[0]*TMath::Landau(x[0],par[1],par[2]);
  //  else if ( x[0]>=75 ) return  par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]);
  //  else return 0;
}


double fitlandauUA2_test(double *x, double *par){
  if ( x[0]>50 && x[0]<75 ) return par[0]*TMath::Landau(x[0],par[1],par[2]);
  else if ( x[0]>=75 ) return  par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]);
  else return 0;
}


double fitlandauUA2_test_tot(double *x, double *par){
  const float pi = 3.1415926;
  double argom =0;
  if (par[2] !=0) argom = (x[0] - par[8])/par[9];
  if ( x[0]>50 && x[0]<75 ) return par[0]*TMath::Landau(x[0],par[1],par[2]) + par[7]*exp(-0.5*argom*argom)/(sqrt(2*pi)*par[9]);
  else if ( x[0]>=75 ) return  par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]) + par[7]*exp(-0.5*argom*argom)/(sqrt(2*pi)*par[9]);
  else return 0;
}


double fitlandauUA2_rejectRange(double *x, double *par){
  // Reject points in Z dimass region
  if (x[0] > 75 && x[0] < 120) {
    TF1::RejectPoint();
    return 0;
  }
  //  return par[0]*TMath::Landau(x[0],par[1],par[2]) + par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]);
  return par[7] * (par[0]*TMath::Landau(x[0],par[1],par[2]) + par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]));
}


double fitlandauUA2_tot(double *x, double *par){
  // Signal and BG function
  const float pi = 3.1415926;
  double argom =0;
  //  if (par[2] !=0) argom = (x[0] - par[8])/par[9];
  //  return par[0]*TMath::Landau(x[0],par[1],par[2]) + par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0]) 
  //    + par[7]*exp(-0.5*argom*argom)/(sqrt(2*pi)*par[9]);

  if (par[2] !=0) argom = (x[0] - par[9])/par[10];
  return par[7] * (par[0]*TMath::Landau(x[0],par[1],par[2]) + par[3]*pow(x[0],-par[4])*exp(-par[5]*x[0])*exp(-par[6]*x[0]*x[0])) 
    + par[8]*exp(-0.5*argom*argom)/(sqrt(2*pi)*par[10]);
}


double fit3exp(double *x, double *par){
  return par[0]*exp(par[1]*x[0]) * par[2]*exp(par[3]*x[0]) * par[4]*exp(par[5]*x[0]);
}


Double_t langaufun(Double_t *x, Double_t *par) {
  // Landau convoluted to gaussian (from ROOT tutorial):
  //----------------------------------------------------
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}


Double_t poigaufun(Double_t *x, Double_t *par) {
  // Poisson convoluted to Gaussian
  //-------------------------------
  //Fit parameters:
  //par[0]: lambda poisson
  //par[1]: Total area (integral -inf to inf, normalization constant)
  //par[2]: Width (sigma) of convoluted Gaussian function

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t fpoi;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  // Range of convolution integral
  xlow = x[0] - sc * par[2];
  xupp = x[0] + sc * par[2];

  step = (xupp-xlow) / np;

  // Convolution integral of Poisson and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fpoi = TMath::Poisson(xx,par[0]);
    //	 fpoi = exp(xx*-par[0]);
    sum += fpoi * TMath::Gaus(x[0],xx,par[2]);

    xx = xupp - (i-.5) * step;
    fpoi = TMath::Poisson(xx,par[0]);
    //	 fpoi = exp(xx*-par[0]);
    sum += fpoi * TMath::Gaus(x[0],xx,par[2]);
  }

  return (par[1] * step * sum * invsq2pi / par[2]);
}


TF1 *poigaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Fit parameters:
   // par[0]: lambda poisson
   // par[1]: Total area (integral -inf to inf, normalization constant)
   // par[2]: Width (sigma) of convoluted Gaussian function

   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[3]  reasonable start values for the fit
   //   parlimitslo[3]  lower parameter limits
   //   parlimitshi[3]  upper parameter limits
   //   fitparams[3]    returns the final fit parameters
   //   fiterrors[3]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn2_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,poigaufun,fitrange[0],fitrange[1],3);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Lambda","Area","GSigma");
   
   for (i=0; i<3; i++) {
     //      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<3; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
}


double fitPIV(double *x, double *par){
  // Pearson IV
  //-----------
  return par[0]*pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2]));
}


double fitPIV_rp(double *x, double *par){
  if (x[0]>50 && x[0]<80) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0]*pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2]));
}


double fitUA2PIV(double *x, double *par){
  // Pearson IV + UA2 functions
  //---------------------------
  double y = par[0]*pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2])) + par[5]*pow(x[0],-par[6])*exp(-par[7]*x[0])*exp(-par[8]*x[0]*x[0]);
  if ( y>0 ) return y;
  else return 0;
}


double fitpol5PIV(double *x, double *par){
  // Pearson IV + pol5
  //------------------
  const int Npow = 5;
  double fpiv = par[0]*pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2])); // PIV
  double fpol = 0; 
  for ( int i=0;i<=Npow;i++ ){
    fpol += par[5+i] * pow(x[0],i); // pol 5
  }
  double y = fpiv + fpol;
  if ( y>0 ) return y;
  else return 0;
}


//double fiterfPIV(double *x, double *par){
//  // Pearson IV + erf
//  //-----------------
//  double fpiv = par[0]*pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2])); // PIV
//  double ferf = par[5]+par[6]*TMath::Erf((x[0]-par[7])/par[8]); // erf
//  double y = fpiv + ferf;
//  if ( y>0 ) return y;
//  else return 0;
//}

double fiterfPIV(double *x, double *par){
  // Pearson IV + erf (with normalization)
  //--------------------------------------
  double fpiv = pow(1+(pow((-x[0]-par[1])/par[2],2)),-par[3])*exp(-par[4]*atan((-x[0]-par[1])/par[2])); // PIV
  double ferf = par[5]+par[6]*TMath::Erf((x[0]-par[7])/par[8]); // erf
  double y = par[0] * (fpiv + ferf);
  if ( y>0 ) return y;
  else return 0;
}

double funcerfPIV3gausum(double *x, double *par){ 
  // Function sum of fit3gausum and fiterfPIV (for display purpose mostly)
  //----------------------------------------------------------------------
  // Sum of 3 gaussian
  //------------------
  const float pi = 3.1415926;
  double arg1 = 0;
  double arg2 = 0;
  double arg3 = 0;
  if ( par[2]!=0 ) arg1 = (x[0] - par[1])/par[2];
  if ( par[5]!=0 ) arg2 = (x[0] - par[4])/par[5];
  if ( par[8]!=0 ) arg3 = (x[0] - par[7])/par[8];
  double fitval1 = par[0] * exp(-0.5*arg1*arg1)/(sqrt(2*pi)*par[2]);
  double fitval2 = par[3] * exp(-0.5*arg2*arg2)/(sqrt(2*pi)*par[5]);
  double fitval3 = par[6] * exp(-0.5*arg3*arg3)/(sqrt(2*pi)*par[8]);

  // Pearson IV + erf (with normalization)
  //--------------------------------------
  double fpiv = pow(1+(pow((-x[0]-par[11])/par[12],2)),-par[13])*exp(-par[14]*atan((-x[0]-par[11])/par[12])); // PIV
  double ferf = par[15]+par[16]*TMath::Erf((x[0]-par[17])/par[18]); // erf
  double y = par[10] * (fpiv + ferf);
  if ( y>0 ) return y + par[9] * (fitval1 + fitval2 + fitval3);
  else return 0;
}

void fitXsig(TH1F *hist, float X){
  // Gaussian fit function
  //   - fits +- X-sigma around mass peak
  //-------------------------------------

  double mean, sigma;

  // Get Maximum
  //------------
  int maxbin = hist->GetMaximumBin();
  mean = hist->GetBinCenter(maxbin);

  // First fit
  //----------
  double fitmin=0.8*mean;
  double fitmax=1.2*mean;
  if ( fabs(mean)<10 ) {
    fitmin= -50;
    fitmax= 50;
  }

  hist->Fit("gaus","RO","",fitmin,fitmax);
  mean = hist->GetFunction("gaus")->GetParameter(1);
  sigma = hist->GetFunction("gaus")->GetParameter(2);

  // Second fit
  //-----------
  fitmin = mean - X*sigma;
  fitmax = mean + X*sigma;
  TF1 *fitfunc = new TF1("fitfunc","gaus");
  fitfunc->SetRange(fitmin,fitmax);
  hist->Fit("fitfunc","RO");
  mean = hist->GetFunction("fitfunc")->GetParameter(1);
  sigma = hist->GetFunction("fitfunc")->GetParameter(2);

  // Third fit
  //----------
  fitmin = mean - X*sigma;
  fitmax = mean + X*sigma;
  fitfunc->SetRange(fitmin,fitmax);
  hist->Fit("fitfunc","RO");
}

#endif // VBFHZZLLBBFITFUNCTIONS_H
