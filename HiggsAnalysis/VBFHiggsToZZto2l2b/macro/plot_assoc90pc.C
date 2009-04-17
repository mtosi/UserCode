{
  
  TH1D * histo = (TH1D*)jetParton_deltaR2_zoom->Clone();
  int nbins = histo->GetNbinsX();
  double xmin = histo->GetXaxis()->GetXmin();
  double xmax = histo->GetXaxis()->GetXmax();
  TString histoName = histo->GetName();
  TString histoTitle = histo->GetTitle();
  TH1D * histo90pc = new TH1D(histoName+"90pc","",nbins,xmin,xmax);

  histo90pc->SetFillColor(kRed);
  histo90pc->SetLineColor(kRed);

  int bin90pc = nbins;
  double totalIntegral = histo->Integral();
  std::cout << "totalIntegral : " << totalIntegral << std::endl;
  
  for ( int ibin = 1; ibin <= nbins ; ibin++ ) {
    std::cout << "nbins: " << nbins << " <--> ibin: " << ibin << std::endl;
    double iIntegral = histo->Integral(1,ibin);
    double frac = iIntegral/totalIntegral;
    std::cout << ibin << " frac: " << frac << std::endl;
    if ( frac <= 0.90 ) bin90pc = ibin;
    std::cout << "bin90pc: " << bin90pc << std::endl;
  }

  for ( int ibin = 1; ibin <= bin90pc; ibin++ ) {
    double icontent = histo->GetBinContent(ibin);
    histo90pc->SetBinContent(ibin,icontent);
    std::cout << "ibin: " << ibin << std::endl;
    std::cout << "histo: " << histo->GetBinContent(ibin) << " <--> histo90pc: " << histo90pc->GetBinContent(ibin) << std::endl;
  }

  TCanvas * G = new TCanvas ("G", "", 600, 600 );
  G->cd();
  histo->Draw();
  histo90pc->Draw("same");


//  G->Print("plot_EtRes.jpg");

}
