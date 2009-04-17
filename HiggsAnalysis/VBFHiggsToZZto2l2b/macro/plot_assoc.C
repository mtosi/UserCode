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

  
  TH1F * A1 = new TH1F ("A1", "Ratio of matches for #DeltaR<0.5 and #DeltaE_{T}<10%", 41, 9.5, 50.5);
  TH1F * A2 = new TH1F ("A2", "Ratio of matches for #DeltaR<0.5 and #DeltaE_{T}<20%", 41, 9.5, 50.5);
  TH1F * A3 = new TH1F ("A3", "Ratio of matches for #DeltaR<0.5 and #DeltaE_{T}<30%", 41, 9.5, 50.5);
  TH1F * A4 = new TH1F ("A4", "Ratio of matches for #DeltaR<0.5", 41, 9.5, 50.5);

  TH1F * B1 = new TH1F ("B1", "Ratio of matches for #DeltaR<0.2 and #D#E_{T}<10%", 41, 9.5, 50.5);
  TH1F * B2 = new TH1F ("B2", "Ratio of matches for #DeltaR<0.2 and #D#E_{T}<20%", 41, 9.5, 50.5);
  TH1F * B3 = new TH1F ("B3", "Ratio of matches for #DeltaR<0.2 and #D#E_{T}<30%", 41, 9.5, 50.5);
  TH1F * B4 = new TH1F ("B4", "Ratio of matches for #DeltaR<0.2", 41, 9.5, 50.5);

  A1->Sumw2();
  A2->Sumw2();
  A3->Sumw2();
  A4->Sumw2();
  B1->Sumw2();
  B2->Sumw2();
  B3->Sumw2();
  B4->Sumw2();

  A1->Divide(Bothok_05_10pc,Allevents);
  A2->Divide(Bothok_05_20pc,Allevents);
  A3->Divide(Bothok_05_30pc,Allevents);
  A4->Divide(Bothok_05,Allevents);
  B1->Divide(Bothok_02_10pc,Allevents);
  B2->Divide(Bothok_02_20pc,Allevents);
  B3->Divide(Bothok_02_30pc,Allevents);
  B4->Divide(Bothok_02,Allevents);

  A4->SetMinimum(0);
  B4->SetMinimum(0);

  A1->SetLineColor(kRed);
  A2->SetLineColor(kBlue);
  A3->SetLineColor(kGreen);
  A4->SetLineColor(kBlack);
  A1->SetMarkerColor(kRed);
  A2->SetMarkerColor(kBlue);
  A3->SetMarkerColor(kGreen);
  A4->SetMarkerColor(kBlack);
  A1->SetMarkerStyle(20);
  A2->SetMarkerStyle(21);
  A3->SetMarkerStyle(24);
  A4->SetMarkerStyle(25);
  A1->SetMarkerSize(0.4); 
  A2->SetMarkerSize(0.4); 
  A3->SetMarkerSize(0.4); 
  A4->SetMarkerSize(0.4); 
  B1->SetLineColor(kRed);
  B2->SetLineColor(kBlue);
  B3->SetLineColor(kGreen);
  B4->SetLineColor(kBlack);
  B1->SetMarkerColor(kRed);
  B2->SetMarkerColor(kBlue);
  B3->SetMarkerColor(kGreen);
  B4->SetMarkerColor(kBlack);
  B1->SetMarkerStyle(20);
  B2->SetMarkerStyle(21);
  B3->SetMarkerStyle(24);
  B4->SetMarkerStyle(25);
  B1->SetMarkerSize(0.4); 
  B2->SetMarkerSize(0.4); 
  B3->SetMarkerSize(0.4); 
  B4->SetMarkerSize(0.4); 
  
 TLegend *legend05 = new TLegend(0.45,0.45,0.99,0.6);
 legend05->SetFillColor(0);
 legend05->SetTextFont(72);
 legend05->SetTextSize(0.04);
 legend05->AddEntry(A4->Clone(),"#DeltaR<0.5",                   "lfp");
 legend05->AddEntry(A3->Clone(),"#DeltaR<0.5 && #DeltaE_{T}<30%","lfp");
 legend05->AddEntry(A2->Clone(),"#DeltaR<0.5 && #DeltaE_{T}<20%","lfp");
 legend05->AddEntry(A1->Clone(),"#DeltaR<0.5 && #DeltaE_{T}<10%","lfp");

 TLegend *legend02 = new TLegend(0.45,0.45,0.99,0.6);
 legend02->SetFillColor(0);
 legend02->SetTextFont(72);
 legend02->SetTextSize(0.04);
 legend02->AddEntry(B4->Clone(),"#DeltaR<0.2",                   "lfp");
 legend02->AddEntry(B3->Clone(),"#DeltaR<0.2 && #DeltaE_{T}<30%","lfp");
 legend02->AddEntry(B2->Clone(),"#DeltaR<0.2 && #DeltaE_{T}<20%","lfp");
 legend02->AddEntry(B1->Clone(),"#DeltaR<0.2 && #DeltaE_{T}<10%","lfp");

  TCanvas * G = new TCanvas ("G", "Probability of good matches as a function of minimum jet Et", 600, 600 );
  G->Divide(2,1);

  G->cd(1);
  A4->GetXaxis()->SetTitle("E_{T} cut (GeV)");
  A4->Draw("PE");
  A3->Draw("PESAME");
  A2->Draw("PESAME");
  A1->Draw("PESAME");
  legend05->Draw();
  G->cd(2);
  B4->GetXaxis()->SetTitle("E_{T} cut (GeV)");
  B4->Draw("PE");
  B3->Draw("PESAME");
  B2->Draw("PESAME");
  B1->Draw("PESAME");
  legend02->Draw();

  G->Print("plot_assoc.jpg");
  G->Print("plot_assoc.eps");

}
