{

    // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

/*  
  TH1D * histo = (TH1D*)jetParton_deltaRVSdeltaR->Clone();
  TProfile * profile = (TProfile*)jetParton_deltaRVSdeltaR_profile->Clone();
  profile->SetMarkerColor(kRed);
  profile->SetLineColor(kRed);

  double xmin = histo->GetXaxis()->GetXmin();
  double xmax = histo->GetXaxis()->GetXmax();

  TF1 *line = new TF1("line","[0]+[1]*x",xmin,xmax);
  line->SetParName(0,"constant");
  line->SetParameter(0,0);
  line->SetParName(1,"slope");
  line->SetParameter(1,1);

  profile->Fit(line,"","",0.8,2.);

  TLegend *legend = new TLegend(0.6,0.2,0.99,0.3);
  legend->SetFillColor(0);
  legend->SetTextFont(72);
  legend->SetTextSize(0.03);
  char nev[50];
  for ( int index = 0; index < line->GetNpar(); index++ ) {
    std::cout << line->GetParName(index) << ": " << line->GetParameter(index) << std::endl;
    sprintf(nev,"%s: %.2f #pm %.2f",line->GetParName(index),line->GetParameter(index),line->GetParError(index));
    legend->AddEntry(line->Clone(),nev,"");
  }

  TCanvas * G = new TCanvas ("G", "", 600, 600 );
  G->cd();
  histo->Draw();
  profile->Draw("same");
  legend->Draw();

  G->Print("plot_deltaR.jpg");
  G->Print("plot_deltaR.eps");
*/

//  TH1D * partonsDeltaR = (TH1D*)ZpartonsDeltaR->Clone();
  TH1D * partonsDeltaR = (TH1D*)TagPartonsDeltaR->Clone();
  partonsDeltaR->GetXaxis()->SetTitle("#DeltaR_{qq}");
  double ymax = partonsDeltaR->GetMaximum();
  double xmin = partonsDeltaR->GetXaxis()->GetXmin();
  double xmax = partonsDeltaR->GetXaxis()->GetXmax();
  int nbin = partonsDeltaR->GetNbinsX();
  double xbin = (xmax-xmin)/double(nbin);

  std::cout << "xbin: " << xbin << std::endl;
  int bin05 = int(0.5/xbin);
  std::cout << "integral: " << partonsDeltaR->Integral(0,bin05) << std::endl;

  TH1D * partonsDeltaR05 = new TH1D("partonsDeltaR05",partonsDeltaR->GetTitle(),bin05,xmin,0.5);
  for ( int ibin = 0; ibin <= bin05; ibin++ ) 
    partonsDeltaR05->SetBinContent(ibin,partonsDeltaR->GetBinContent(ibin));
  partonsDeltaR05->SetLineColor(kRed);
  partonsDeltaR05->SetFillColor(kRed);

  TH1D * TagPartonsDeltaR = (TH1D*)TagPartonsDeltaR->Clone();

//  TLegend *legend = new TLegend(0.4,0.7,0.85,0.8);
//  TLegend *legend = new TLegend(0.15,0.7,0.45,0.8);
  TLegend *legend = new TLegend(0.6,0.7,0.95,0.8);
  legend->SetFillColor(0);
  legend->SetTextFont(72);
  legend->SetTextSize(0.03);
  char nev[50];
  sprintf(nev,"entries: %d",partonsDeltaR->GetEntries());
  legend->AddEntry(partonsDeltaR,nev,"l");
  sprintf(nev,"fraction w/ #DeltaR<0.5: %.0f",double(partonsDeltaR05->Integral()/partonsDeltaR->GetEntries()));
  legend->AddEntry(partonsDeltaR05,nev,"f");


  TCanvas * G = new TCanvas ("G", "", 600, 600 );
  G->cd();
  partonsDeltaR->Draw();
  partonsDeltaR05->Draw("same");
  legend->Draw("same");

  G->Print("plot_deltaR.jpg");
  G->Print("plot_deltaR.eps");

}
