{

  
  TH1D * histo1 = (TH1D*)jetParton_004deltaEtRes->Clone();

  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(0.4); 

  TF1 * gaus = TF1("gaus","gaus");
  gaus->SetLineColor(kRed);

  histo1->Fit("gaus");
  
  TCanvas * G = new TCanvas ("G", "", 600, 600 );
  G->cd();
  histo1->Draw();

//  G->Print("plot_EtRes.jpg");

}
